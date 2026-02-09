%% run_hfmf_5LG_single_point.m
% Reproduce the Python HFMF (single-point SCF) for 5-layer graphene 2-band effective model
% Units follow the original Python by default:
%   Energy: meV, Length: cm, ke = 1.44e-4 meV*cm, a = 2.46e-8 cm, beta=58 1/meV (T~0.2K)

clear; clc;

%% ===================== User knobs (match Python) =====================
N_layer  = 5;

% k-mesh (dimensionless, unit of 1/a)
dk       = 1.0e-3;
kmax     = 0.150;
c3_mesh  = true;      % True -> 60-degree mesh like Python
fix_c3   = true;     % optional (need C3op), default off

% Coulomb & geometry (cm, meV)
a_lat    = 2.46e-8;     % cm
ke       = 1.44e-4;     % meV*cm
er       = 40;
d_gate   = 3.69e-6;     % cm
d_lay    = 3.35e-8;     % cm  (interlayer distance)

% HFMF controls
alp      = 0.03;        % valley-interchange strength
beta     = 58;          % 1/meV  (T ~ 0.2K)
mix      = 0.30;
tol      = 3e-8;
max_iter = 200;

% external field (the "U" in sweep), in meV (consistent with the Python code)
U_ext    = 0.0;

% density ne in cm^-2 (Python uses e.g. 0.00e12)
ne       = 0.00e12;     % cm^-2

% LAFz seed strength (Python seed_val: +/-5). This is an energy scale in meV.
m_seed   = 5.0;

%% ===================== Build k-mesh & area =====================
[kx, ky, c3_x, c3_y, A] = make_kmesh(dk, kmax, c3_mesh, a_lat);
nk = size(kx,1);
fprintf('nk=%d, nk^2=%d, A=%.6e (cm^2)\n', nk, nk^2, A);

%% ===================== Build 8-band operators =====================
ops = build_ops_8band(); % eye8, Sx,Sz,Vx,Vz,Lx,Lz

%% ===================== SWM parameters (meV) =====================
SWM = get_SWM_meV(N_layer);

%% ===================== Build V_r (nk x nk x 8 x 8) =====================
V_r = build_Vr_python_units(kx, ky, ke, er, d_gate, d_lay, N_layer, a_lat);

%% ===================== Build H0(k) and HU(k) on the mesh =====================
SOC        = 0.0;
SOC_dir    = 0;
SOC_single = 1;

H0_k = build_H0_k(kx, ky, SWM, N_layer, SOC, SOC_dir, SOC_single, ops);
HU_k = build_HU_k(kx, ky, SWM, N_layer);   % diagonal only, as in Python sweep()

Hn_k = H0_k + U_ext * HU_k;                % non-interacting part used in SCF loop

%% ===================== Build LAFz seed_V =====================
% Python seed_val: K_up=+m, K_dn=-m, Kp_up=+m, Kp_dn=-m  -> seed = m * S_z * L_z
seed_V = m_seed * (ops.Sz * ops.Lz);   % (8x8)
%% ===================== Build ValleyHall seed_V =====================
% Python seed_val: K_up=+m, K_dn=-m, Kp_up=+m, Kp_dn=-m  -> seed = m * S_z * L_z
seed_V = m_seed * (ops.eye8 * ops.Lz);   % (8x8)

%% ===================== Initial diagonalization: H_seed = Hn + seed_V =====================
H_seed_k = Hn_k + reshape(seed_V,1,1,8,8);

[E_seed, Vec_seed] = diag_batch_hermitian(H_seed_k);

mu = solve_mu_bisect(E_seed, ne, A, beta);
occ = fermi_dirac_safe(E_seed, mu, beta);
rho_old = rho_from_eig(Vec_seed, occ);
rho_new = rho_old;

fprintf('Init: mu=%.6f meV\n', mu);

%% ===================== SCF loop =====================
for it = 1:max_iter
    rho_mix = mix * rho_new + (1-mix) * rho_old;

    % Build HF potential via FFT (Fock + VI + Hartree)  -> V_F(k)
    V_F_k = V_HF_builder_fft(rho_mix, V_r, A, alp, ke, er, d_gate, d_lay, N_layer);

    % Optional C3 symmetrization of V_F(k) (needs C3op; default off)
    if fix_c3
        C3op = eye(8); % placeholder; set to correct internal C3 operator if you use this
        V_F_k = C3_symmetrize_mat(nk, c3_x, c3_y, V_F_k, C3op);
    end

    % Effective Hamiltonian
    H_int_k = Hn_k + V_F_k;

    % Diagonalize
    [E_int, Vec_int] = diag_batch_hermitian(H_int_k);

    % Update mu and rho
    mu = solve_mu_bisect(E_int, ne, A, beta);

    occ = fermi_dirac_safe(E_int, mu, beta);
    rho_old = rho_new;
    rho_new = rho_from_eig(Vec_int, occ);

    % Energies (match Python)
    [E_tot, E_tot0] = total_energies(H_int_k, Hn_k, rho_new);

    % Error metric (match Python style)
    err = density_error(rho_new, rho_old);

    fprintf('Iter %3d: Err=%.3e, mu=%.6f meV, Etot=%.10e, Etot0=%.10e\n', ...
        it, err, mu, E_tot, E_tot0);

    if err < tol
        break;
    end
end

%% ===================== Diagnostics: LAF order parameter =====================
mLAF = order_parameter(rho_new, ops.Sz*ops.Lz);
fprintf('Order <Sz*Lz> per k-pt (trace averaged): %.6e\n', mLAF);

disp('Done.');
%%
figure()
hold on;
% band1=E_int(:,:,1);
% band2=E_int(:,:,5);
band1=E_seed(:,:,1);
band2=E_seed(:,:,5);
surf(kx,ky,band1,'EdgeColor','none')
surf(kx,ky,band2,'EdgeColor','none')
%% ========================================================================
%%                           Local functions
%% ========================================================================

function SWM = get_SWM_meV(N)
% SWM parameters in meV (exactly as Python snippet)
    switch N
        case 2
            SWM.gamma0 = 3160;  SWM.gamma1 = 500;  SWM.gamma2 =   0.00;
            SWM.gamma3 = -280;  SWM.gamma4 = -200; SWM.delta  =  -0.05;
        case 3
            SWM.gamma0 = 3160;  SWM.gamma1 = 460;  SWM.gamma2 = -17.0;
            SWM.gamma3 = -300;  SWM.gamma4 =  -86; SWM.delta  =  -1.1;
        case 4
            SWM.gamma0 = 3160;  SWM.gamma1 = 445;  SWM.gamma2 = -18.2;
            SWM.gamma3 = -319;  SWM.gamma4 =  -79; SWM.delta  =  -0.066;
        case 5
            SWM.gamma0 = 3160;  SWM.gamma1 = 435;  SWM.gamma2 = -18.5;
            SWM.gamma3 = -322;  SWM.gamma4 = -67.5;SWM.delta  =  -0.147;
        case 6
            SWM.gamma0 = 3160;  SWM.gamma1 = 430;  SWM.gamma2 = -18.5;
            SWM.gamma3 = -325;  SWM.gamma4 =  -73; SWM.delta  =  -1.6;
        otherwise
            error('N supported: 2-6');
    end
end

function ops = build_ops_8band()
% 8 band = spin(up,dn) x valley(K,K') x orbital(1A,NB)
    sig0 = eye(2);
    sigx = [0 1; 1 0];
    sigz = [1 0; 0 -1];

    ops.eye8 = kron(sig0, kron(sig0, sig0));
    ops.Sx   = kron(sigx, kron(sig0, sig0));
    ops.Sz   = kron(sigz, kron(sig0, sig0));
    ops.Vx   = kron(sig0, kron(sigx, sig0));
    ops.Vz   = kron(sig0, kron(sigz, sig0));
    ops.Lx   = kron(sig0, kron(sig0, sigx));
    ops.Lz   = kron(sig0, kron(sig0, sigz));
end

function [kx, ky, c3_x, c3_y, A] = make_kmesh(dk, kmax, c3_mesh, a_lat_cm)
% Build k mesh in dimensionless unit (1/a), and A in cm^2 as Python:
%   A = (2*pi*a/dk)^2; if c3_mesh -> A *= 2/sqrt(3)
    nk = round(2*kmax/dk) + 1;
    t = linspace(-kmax, kmax, nk);
    [tgx, tgy] = meshgrid(t, t);

    if c3_mesh
        kx = tgx + 0.5 * tgy;
        ky = (sqrt(3)/2) * tgy;
        A  = (2*pi*a_lat_cm/dk)^2 * (2/sqrt(3));

        % C3 partner indices (Python-style)
        % Build integer grid coords centered at 0 using 0-based logic then convert to 1-based.
        [I,J] = meshgrid(0:nk-1, 0:nk-1); % NOTE: MATLAB meshgrid gives X=cols, Y=rows; we'll follow Python mapping later
        tr_x0 = I - (nk-1)/2;
        tr_y0 = J - (nk-1)/2;

        tr_x1 = tr_x0 + 0.5 * tr_y0;
        tr_y1 = (sqrt(3)/2) * tr_y0;

        tr_x2 = tr_x1 * (-1/2)         + tr_y1 * (-sqrt(3)/2);
        tr_y2 = tr_x1 * ( sqrt(3)/2)   + tr_y1 * (-1/2);

        c3_x0 = round(tr_x2 - (1/sqrt(3)) * tr_y2 + (nk-1)/2);
        c3_y0 = round((2/sqrt(3)) * tr_y2 + (nk-1)/2);

        % convert to 1-based for MATLAB indexing
        c3_x = c3_x0 + 1;
        c3_y = c3_y0 + 1;
    else
        kx = tgx;
        ky = tgy;
        A  = (2*pi*a_lat_cm/dk)^2;

        c3_x = [];
        c3_y = [];
    end
end

function V_r = build_Vr_python_units(kx, ky, ke, er, d_gate, d_lay, N, a_lat)
% Match Python:
%   qq = sqrt(kx^2+ky^2)/a + eps   (cm^-1)
%   V_S_r = 2π ke/er * fft2( 1/qq * (cosh(2qq d_gate)-cosh(qq(N-1)d_lay))/sinh(2qq d_gate) )
%   V_D_r = 2π ke/er * fft2( 1/qq * (cosh(2qq d_gate-qq(N-1)d_lay)-1)/sinh(2qq d_gate) )
%   V_r   = kron(ones(4,4), [[V_S_r,V_D_r],[V_D_r,V_S_r]])  -> (nk,nk,8,8)
    q = sqrt(kx.^2 + ky.^2) ./ a_lat + 1e-15; % cm^-1

    kernel_S = (1 ./ q) .* (cosh(2*q*d_gate) - cosh(q*(N-1)*d_lay)) ./ sinh(2*q*d_gate);
    kernel_D = (1 ./ q) .* (cosh(2*q*d_gate - q*(N-1)*d_lay) - 1) ./ sinh(2*q*d_gate);

    V_S_r = 2*pi*ke/er * fft2(kernel_S);
    V_D_r = 2*pi*ke/er * fft2(kernel_D);

    [nk,~] = size(kx);
    V_r = zeros(nk, nk, 8, 8);

    % orbital block (2x2)
    V_orb = zeros(nk,nk,2,2);
    V_orb(:,:,1,1) = V_S_r;
    V_orb(:,:,2,2) = V_S_r;
    V_orb(:,:,1,2) = V_D_r;
    V_orb(:,:,2,1) = V_D_r;

    % expand to 4 flavor (spin×valley)
    Norb = 2; Nfl = 4;
    for f1 = 0:Nfl-1
        for f2 = 0:Nfl-1
            for o1 = 1:Norb
                for o2 = 1:Norb
                    i = f1*Norb + o1;
                    j = f2*Norb + o2;
                    V_r(:,:,i,j) = V_orb(:,:,o1,o2);
                end
            end
        end
    end
end

function H0_k = build_H0_k(kx, ky, SWM, N, SOC, SOC_dir, SOC_single, ops)
% H0(k) = kron(I_spin, blkdiag(H2B(K), H2B(K'))) + Ising_SOC
    [nk,~] = size(kx);
    H0_k = complex(zeros(nk,nk,8,8));

    for ix = 1:nk
        for iy = 1:nk
            H2K  = H2B_one_flavor(kx(ix,iy), ky(ix,iy), +1, SWM, N);
            H2Kp = H2B_one_flavor(kx(ix,iy), ky(ix,iy), -1, SWM, N);
            H_spin_valley = kron(eye(2), blkdiag(H2K, H2Kp));

            % Ising SOC (same structure as Python)
            Ising_SOC = SOC * (ops.Vz) * (SOC_dir * ops.Sx + (1 - SOC_dir) * ops.Sz) * ...
                ( (1 - SOC_single) * ops.Lz + SOC_single * (ops.eye8 - ops.Lz) / 2 );

            H0_k(ix,iy,:,:) = H_spin_valley + Ising_SOC;
        end
    end
end

function H = H2B_one_flavor(kx, ky, q, SWM, N)
% This matches the earlier Python H_2B (single flavor) structure.
    hv0a  = sqrt(3)/2 * SWM.gamma0;
    hv3a  = sqrt(3)/2 * SWM.gamma3;
    hv4a  = sqrt(3)/2 * SWM.gamma4;
    gm1   = SWM.gamma1;
    gm2   = SWM.gamma2;
    delta = SWM.delta;

    Pi = q*kx + 1i*ky;

    H = complex(zeros(2,2));

    x2 = abs(hv0a * Pi / gm1).^2;
    rN   = (x2.^N     - 1) ./ (x2 - 1);
    rNm1 = (x2.^(N-1) - 1) ./ (x2 - 1);

    Pi3_plus = Pi.^3 + conj(Pi).^3;

    % linear v3 correction in rN (as in the H_2B snippet)
    % rN = rN - hv0a^2 * hv3a / gm1^3 .* real(Pi3_plus) .* (1 + 2*x2 + 3*x2.^2);

    % H_ch
    H(1,2) = H(1,2) + (-gm1) * (hv0a * conj(Pi) / (-gm1))^N;
    H(2,1) = H(2,1) + (-gm1) * (hv0a * Pi        / (-gm1))^N;

    % H_s
    H(1,1) = H(1,1) + delta - 2 * abs(Pi)^2 * hv0a * hv4a / gm1 * rNm1;
    H(2,2) = H(2,2) + delta - 2 * abs(Pi)^2 * hv0a * hv4a / gm1 * rNm1;

    % H_tr
    pref_tr = ((N-2)*gm2/2 - (N-1)*hv0a*hv3a*abs(Pi)^2/gm1);
    H(1,2) = H(1,2) + pref_tr * (hv0a * conj(Pi) / (-gm1))^(N-3);
    H(2,1) = H(2,1) + pref_tr * (hv0a * Pi        / (-gm1))^(N-3);

    % second order in v3 (pentalayer note)
    % pref_v3_2 = (3*hv0a*hv3a^2*abs(Pi)^2/gm1^2 - hv3a*gm2/gm1);
    % H(1,2) = H(1,2) + pref_v3_2 * Pi;
    % H(2,1) = H(2,1) + pref_v3_2 * conj(Pi);

    % pref_diag = (hv0a*hv3a*hv4a/gm1^2) * real(Pi3_plus) * (1 + 3*x2 + 5*x2.^2);
    % H(1,1) = H(1,1) + pref_diag;
    % H(2,2) = H(2,2) + pref_diag;

    H = H ./ rN;
end

function HU_k = build_HU_k(kx, ky, SWM, N)
% Match your Python sweep() H_U construction.
% HU is diagonal; built for first 4 bands then copied to 5:8.
    hv0a = sqrt(3)/2 * SWM.gamma0;
    hv3a = sqrt(3)/2 * SWM.gamma3;
    gm1  = SWM.gamma1;

    [nk,~] = size(kx);
    HU_k = zeros(nk,nk,8,8);

    for ix = 1:nk
        for iy = 1:nk
            % spin-up block indices: 1..4 (valley K:1-2, valley K':3-4)
            for vv = 1:2
                if vv == 1
                    q  = +1;
                    i1 = 1; i2 = 2;
                else
                    q  = -1;
                    i1 = 3; i2 = 4;
                end

                Pi = q*kx(ix,iy) + 1i*ky(ix,iy);

                abs_term = abs(hv0a * conj(Pi) / gm1); % scalar

                % r_gap = sum_{i=0}^{N-2} (1 - 2 i/(N-1)) * abs_term^(2i)
                r_gap = 0;
                for ii = 0:(N-2)
                    r_gap = r_gap + (1 - 2*ii/(N-1)) * (abs_term^(2*ii));
                end

                % r_N (uses abs_term)
                rN = (abs_term^(2*N) - 1) / (abs_term^2 - 1);

                % pentalayer correction to rN (note factor -2 in your sweep snippet)
                % x2 = abs(hv0a * Pi / gm1)^2;
                % Pi3_plus = Pi^3 + conj(Pi)^3;
                % rN = rN - 2 * hv0a^2 * hv3a / gm1^3 * real(Pi3_plus) * (1 + 2*x2 + 3*x2^2);

                % r_gap2 (pentalayer only)
                % r_gap2 = hv3a * hv0a^2 / gm1^3 * real(Pi3_plus) * ( x2 + 3*x2^2 );

                r_gap2 =0;
                val = 0.5 * (r_gap + r_gap2) / rN;

                HU_k(ix,iy,i1,i1) =  val;
                HU_k(ix,iy,i2,i2) = -val;
            end

            % copy spin-up(1:4) -> spin-down(5:8)
            HU_k(ix,iy,5:8,5:8) = HU_k(ix,iy,1:4,1:4);
        end
    end
end

function V_F_k = V_HF_builder_fft(rho_k, V_r, A, alp, ke, er, d_gate, d_lay, N)
% Match Python V_HF_builder:
%   rho_r = fft2(rho_k)
%   V_F_r = - rho_r .* V_r / A
%   V_F   = ifftshift(ifft2(V_F_r))
% plus valley interchange (VI) and Hartree q->0 terms (k-independent)
    [nk,~,Nb,~] = size(rho_k);
    V_F_k = complex(zeros(nk,nk,Nb,Nb));

    % Fock via FFT (elementwise in band space)
    for a = 1:Nb
        for b = 1:Nb
            rho_r = fft2(rho_k(:,:,a,b));
            VFr   = - rho_r .* V_r(:,:,a,b) / A;
            V_F_k(:,:,a,b) = ifftshift(ifft2(VFr));  % IMPORTANT: match Python ifftshift
        end
    end

    % k-summed density matrix
    rho_k_sum = squeeze(sum(sum(rho_k,1),2));  % [Nb x Nb]

    % ---- valley interchange term (VI) ----
    V_VI = complex(zeros(Nb,Nb));
    if abs(alp) > 0
        C_VI = 2*pi*ke/er * (d_gate - ((N-1)*d_lay)^2/(4*d_gate));
        tV_VI = - alp * C_VI * rho_k_sum / A;

        % multiply by kron(ones(4,4), I2)
        tV_VI = tV_VI .* kron(ones(4,4), eye(2));

        N_BL  = Nb/4;  % =2
        rows1 = [1:N_BL, 2*N_BL+1:3*N_BL];      % blocks 0 & 2 (MATLAB 1-based)
        rows2 = [N_BL+1:2*N_BL, 3*N_BL+1:4*N_BL];% blocks 1 & 3

        V_VI(rows1,rows1) = tV_VI(rows2,rows2);
        V_VI(rows2,rows2) = tV_VI(rows1,rows1);
    end

    % ---- Hartree term (q~0 constants) ----
    diag_D = zeros(Nb,1);
    diag_D(1:2:end) = real(trace(rho_k_sum(2:2:end,2:2:end)));
    diag_D(2:2:end) = real(trace(rho_k_sum(1:2:end,1:2:end)));

    diag_S = zeros(Nb,1);
    diag_S(2:2:end) = diag_D(1:2:end);
    diag_S(1:2:end) = diag_D(2:2:end);

    V_D0 = 2*pi*ke/er*( d_gate - (N-1)*d_lay + ((N-1)*d_lay)^2/(4*d_gate) );
    V_S0 = 2*pi*ke/er*( d_gate               - ((N-1)*d_lay)^2/(4*d_gate) );

    V_H = (V_D0*diag(diag_D) + V_S0*diag(diag_S)) / A;

    % add k-independent VI + Hartree to every k
    for ix = 1:nk
        for iy = 1:nk
            V_F_k(ix,iy,:,:) = squeeze(V_F_k(ix,iy,:,:)) + V_VI + V_H;
        end
    end
end

function [E, VEC] = diag_batch_hermitian(Hk)
% Diagonalize Hk(ix,iy,:,:) for all k
    [nk,~,Nb,~] = size(Hk);
    E   = zeros(nk,nk,Nb);
    VEC = complex(zeros(nk,nk,Nb,Nb));
    for ix = 1:nk
        for iy = 1:nk
            H = squeeze(Hk(ix,iy,:,:));
            H = (H + H')/2; % enforce Hermitian numerically
            [V,D] = eig(H);
            e = real(diag(D));
            [e,idx] = sort(e,'ascend');
            V = V(:,idx);
            E(ix,iy,:) = e;
            VEC(ix,iy,:,:) = V;
        end
    end
end

function occ = fermi_dirac_safe(E, mu, beta)
% Safe FD to avoid overflow: E,mu in meV; beta in 1/meV
    x = beta*(E - mu);
    occ = zeros(size(E));
    occ(x < -40) = 1;
    occ(x >  40) = 0;
    m = (x >= -40) & (x <= 40);
    occ(m) = 1 ./ (1 + exp(x(m)));
end

function rho = rho_from_eig(VEC, occ)
% rho(k) = V diag(occ) V^\dagger
    [nk,~,Nb,~] = size(VEC);
    rho = complex(zeros(nk,nk,Nb,Nb));
    for ix = 1:nk
        for iy = 1:nk
            V = squeeze(VEC(ix,iy,:,:));
            f = squeeze(occ(ix,iy,:));
            rho(ix,iy,:,:) = V * (diag(f) * V');
        end
    end
end

function mu = solve_mu_bisect(E, ne, A, beta)
% Match Python get_mu:
% find_mu(mu) = sum f(E,mu) - nk^2*(Nb/2) - ne*A = 0
    [nk,~,Nb] = size(E);
    Nk_tot = nk*nk;

    Evmax = max(E(:,:,1:Nb/2),[],'all');
    Ecmin = min(E(:,:,Nb/2+1:Nb),[],'all');
    mu0 = 0.5*(Evmax + Ecmin);

    lo = mu0 - 100;  % meV
    hi = mu0 + 100;  % meV

    f_lo = find_mu(lo);
    f_hi = find_mu(hi);

    % expand bracket if needed
    tries = 0;
    while f_lo*f_hi > 0 && tries < 20
        lo = lo - 100;
        hi = hi + 100;
        f_lo = find_mu(lo);
        f_hi = find_mu(hi);
        tries = tries + 1;
    end
    if f_lo*f_hi > 0
        error('Failed to bracket chemical potential. Try larger bracket.');
    end

    for it = 1:200
        mid = 0.5*(lo+hi);
        f_mid = find_mu(mid);
        if f_lo*f_mid <= 0
            hi = mid; f_hi = f_mid;
        else
            lo = mid; f_lo = f_mid;
        end
        if abs(hi-lo) < 1e-10
            break;
        end
    end
    mu = 0.5*(lo+hi);

    function val = find_mu(mu_try)
        occ = fermi_dirac_safe(E, mu_try, beta);
        val = sum(occ,'all') - Nk_tot*(Nb/2) - ne*A;
    end
end

function [Etot, Etot0] = total_energies(H_int_k, H_n_k, rho)
% Match Python:
% Etot  = sum Tr[ 1/2*(H_int + H_n) rho ]
% Etot0 = sum Tr[ H_n rho ]
    [nk,~,Nb,~] = size(H_int_k);
    Etot = 0; Etot0 = 0;
    for ix = 1:nk
        for iy = 1:nk
            Hint = squeeze(H_int_k(ix,iy,:,:));
            Hn   = squeeze(H_n_k(ix,iy,:,:));
            r    = squeeze(rho(ix,iy,:,:));
            Etot  = Etot  + real(trace(0.5*(Hint+Hn)*r));
            Etot0 = Etot0 + real(trace(Hn*r));
        end
    end
end

function err = density_error(rho_new, rho_old)
% Match Python:
% error = sum ||rho_new - rho_old|| / (Nb*nk)^2
    [nk,~,Nb,~] = size(rho_new);
    s = 0;
    for ix = 1:nk
        for iy = 1:nk
            d = squeeze(rho_new(ix,iy,:,:) - rho_old(ix,iy,:,:));
            s = s + norm(d,'fro');
        end
    end
    err = s / ((Nb*nk)^2);
end

function op = order_parameter(rho, M)
% average over k of Tr[ rho(k) M ] / Nk
    [nk,~,~,~] = size(rho);
    s = 0;
    for ix = 1:nk
        for iy = 1:nk
            r = squeeze(rho(ix,iy,:,:));
            s = s + real(trace(r*M));
        end
    end
    op = s / (nk*nk);
end

function V_F = C3_symmetrize_mat(nk, c3_x, c3_y, V_F, C3op)
% Match Python C3_symmetrize
    V0 = V_F;
    for i = 1:nk
        for j = 1:nk
            c1i = c3_y(i,j); c1j = c3_x(i,j);
            if c1i>=1 && c1i<=nk && c1j>=1 && c1j<=nk
                c2i = c3_y(c1i,c1j); c2j = c3_x(c1i,c1j);
                V_F(i,j,:,:) = ( ...
                    squeeze(V0(i,j,:,:)) + ...
                    (C3op') * squeeze(V0(c1i,c1j,:,:)) * C3op + ...
                    C3op * squeeze(V0(c2i,c2j,:,:)) * (C3op') ) / 3;
            end
        end
    end
end
