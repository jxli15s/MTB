clc;
clear;
delete(gcp('nocreate'))
% Nworkers  = 6;   % 外层 parfor worker 数
% Nthreads  = 1;    % 每个 worker 内部 BLAS 线程数
% pctRunOnAll maxNumCompThreads(Nthreads);
parpool('local',6);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RhG Model      %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the structure
g = MTB.geometry("RhG");
a = 2.46;   % Ang
g.a =  [1/2,-sqrt(3)/2,0;...
        1/2,sqrt(3)/2,0;...
        0,0,1/a]*a;

g.b=inv(g.a')*2*pi;
Layer_N = 2;
pars = struct('v0',3.16, 'gamma1',0.46, 'gamma2',-0.017, ...   % eV*Ang (示例数)
              'gamma3',-0.30, 'gamma4',-0.086, ... % eV
              'uext',0.015, 'delta',-0.0011, 'xi',1, 'N',Layer_N);        % eV
% % 
% Layer_N = 3;
% pars = struct('v0',3.16, 'gamma1',0.435, 'gamma2',-0.0185, ...   % eV*Ang (示例数)
%               'gamma3',-0.322, 'gamma4',-0.0675, ... % eV
%               'uext',0.02, 'delta',-0.000147, 'xi',1, 'N',Layer_N);
% get the ham
%kpoint=[0.5,0.0,0.0];
%[h,hx,hy,hz] = MTB.ham.get_ham_kp_RhG(kpoint, pars);
model=@MTB.ham.get_ham_kp_RhG;
% model=@MTB.ham.get_ham_kp_RhG_valley;
g.dim_kp=Layer_N*2;
nbands=g.dim_kp;
knum=301;
tic;
efermi=get_ef(g,pars,model,nbands,knum);
toc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Get the 3D band in plane   %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=300;
kxline=[-0.03,0.03];
kyline=[-0.03,0.03];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
model=@MTB.ham.get_ham_kp_RhG;
nbands=g.dim_kp;
part=(kxline(2)-kxline(1))^2;
knum=size(Kx,1);
weights=1/knum^2/det(g.a)*part*pi*2;

tic;
% [Unk,Enk,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity(pars,model,nbands,Kx,Ky,Kz);
[Unk,Enk,ham,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity_withH(pars,model,nbands,Kx,Ky,Kz);
toc;

efermi=calculate_ef(Enk(:),0.5);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RhG Model      %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the structure
g = MTB.geometry("RhG");
a = 2.46;   % Ang
g.a =  [1/2,-sqrt(3)/2,0;...
        1/2,sqrt(3)/2,0;...
        0,0,1/a]*a;

g.b=inv(g.a')*2*pi;
Layer_N = 2;
pars = struct('v0',3.16, 'gamma1',0.46, 'gamma2',-0.017, ...   % eV*Ang (示例数)
              'gamma3',-0.30, 'gamma4',-0.086, ... % eV
              'uext',0.015, 'delta',-0.0011, 'xi',1, 'N',Layer_N);        % eV
% % 
% Layer_N = 3;
% pars = struct('v0',3.16, 'gamma1',0.435, 'gamma2',-0.0185, ...   % eV*Ang (示例数)
%               'gamma3',-0.322, 'gamma4',-0.0675, ... % eV
%               'uext',0.02, 'delta',-0.000147, 'xi',1, 'N',Layer_N);
% get the ham
%kpoint=[0.5,0.0,0.0];
%[h,hx,hy,hz] = MTB.ham.get_ham_kp_RhG(kpoint, pars);
model=@MTB.ham.get_ham_kp_RhG;
% model=@MTB.ham.get_ham_kp_RhG_valley;
g.dim_kp=Layer_N*2;
nbands=g.dim_kp;
knum=501;
tic;
efermi=get_ef(g,pars,model,nbands,knum);
toc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Get the 3D band in plane   %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=300;
kxline=[-0.03,0.03];
kyline=[-0.03,0.03];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
model=@MTB.ham.get_ham_kp_RhG;
nbands=g.dim_kp;
part=(kxline(2)-kxline(1))^2;
knum=size(Kx,1);
weights=1/knum^2/det(g.a)*part*pi*2;

tic;
% [Unk,Enk,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity(pars,model,nbands,Kx,Ky,Kz);
[Unk,Enk,ham,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity_withH(pars,model,nbands,Kx,Ky,Kz);
toc;

efermi=calculate_ef(Enk(:),0.5);
%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.1, 1000);   % hbar*omega in eV
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 0.0;       % eV (0 => step)
eta = 0.001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.band_list   = 1:4;     % only sum within these bands (speedup)
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;

% opts.tau = [0,0,0];

% [sigma_abc, out] = shift_current_plane_fd_energy_skew( Kx, Ky, Kz, pars, model, Eph_list, Ef, kT, eta, opts);
[sigma_abc,out] = shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, opts);
%%
figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(sigma_abc(i,j,k,:));
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%%
%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.1, 1000);   % hbar*omega in eV
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 0.0;       % eV (0 => step)
eta = 0.001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;
%%
[eta_abc, out] = mic_metric_plane_fromUE_mn(Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, opts);
%%
figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(eta_abc(i,j,k,:))*10^-12*10^6/10;
            % sig_xxy = squeeze(sigma_abc(i,j,k,:))*10^6/10;
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%%
function [eta_abc, out] = mic_metric_plane_fromUE_mn( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%MIC_METRIC_PLANE_FROMUE_MN
% =========================================================================
% Magnetic Injection Current (MIC) for LINEAR polarization (metric-type),
% multi-band form: SUM over all interband pairs (m,n), m≠n, on a 2D k-grid.
%
% Energy (eV) form with Lorentzian delta_eV (unit 1/eV):
%
%   eta^{abc}(Eph) = -(pi * g_s * e^2)/(2*hbar) *
%       ∫ d^2k/(2π)^2  Σ_{m≠n}  f_nm(k) * Δv^a_mn(k) * [2 g^{bc}_mn(k)] *
%       δ_eV( (E_m - E_n) - Eph )
%
% where (linear polarization => symmetric in b,c):
%   2 g^{bc}_{mn} = r^b_{mn} r^c_{nm} + r^c_{mn} r^b_{nm}
%   r^b_{mn} = i <u_m | ∂_{k_b} u_n>   (length unit follows k-unit)
%   v^a_n    = (1/ħ) ∂E_n/∂k_a  (SI: m/s if k in 1/m; here length unit follows k-unit)
%   Δv^a_mn  = v^a_m - v^a_n
%   f_nm     = f_n - f_m   (initial minus final; consistent with ΔE=E_m-E_n)
%
% Implementation notes:
% - All energies (E, Eph_list, Ef, kT, eta) are in eV.
% - Lorentzian broadening parameter 'eta' is in eV.
% - By default we keep only absorption-like transitions ΔE>0 using a mask
%   (branchless, no if inside heavy loops).
%
% Inputs:
%   Kx,Ky     : Nkx x Nky (Cartesian k-grid components, unit 1/Ang or 1/m ...)
%   U         : nb x nb_sel x Nkx x Nky  (eigenvectors; columns are bands)
%               or nb x nb_sel x (Nkx*Nky)
%   E         : Nkx x Nky x nb_sel (eV) or (Nkx*Nky) x nb_sel
%   Eph_list  : 1 x Nw photon energies (eV)
%   Ef,kT,eta : eV
%   opts fields (optional):
%       opts.periodicFD    (default false)
%       opts.trimBoundary  (default true)
%       opts.doGaugeFix    (default true)
%       opts.g_s           (default 1)
%       opts.verbose       (default true)
%       opts.positiveDE    (default true)   % keep only ΔE>0
%       opts.saveFullMN    (default true)   % store large mn-resolved arrays in out
%
% Outputs:
%   eta_abc : 2 x 2 x 2 x Nw  (a,b,c in {x,y})
%   out     : struct with k/band-resolved ingredients (optional)
% =========================================================================

    arguments
        Kx double
        Ky double
        U  {mustBeNumeric}
        E  double
        Eph_list double
        Ef double
        kT double
        eta double
        opts struct = struct()
    end

    % -------------------- options --------------------
    periodicFD   = get_opt(opts,'periodicFD',false);
    trimBoundary = get_opt(opts,'trimBoundary',true);
    doGaugeFix   = get_opt(opts,'doGaugeFix',true);
    g_s          = get_opt(opts,'g_s',1);
    verbose      = get_opt(opts,'verbose',true);
    positiveDE   = get_opt(opts,'positiveDE',true);
    saveFullMN   = get_opt(opts,'saveFullMN',true);

    % -------------------- constants (SI) --------------------
    e_charge = 1.602176634e-19;   % C
    hbar_Js  = 1.054571817e-34;   % J*s

    % MIC prefactor (as in your comment):
    %   -(pi g_s e^2)/(2 ħ)
    pref = -pi * g_s * (e_charge^2) / (2*hbar_Js);

    % velocity prefactor when E is in eV:
    %   v = (1/ħ) ∂E_J/∂k = (e/ħ) ∂E_eV/∂k
    vel_pref = e_charge / hbar_Js;

    % -------------------- reshape & checks --------------------
    [Nkx, Nky] = size(Kx);
    if ~isequal(size(Ky), [Nkx,Nky])
        error('Kx and Ky must have the same size Nkx x Nky.');
    end

    % reshape U to nb x nb_sel x Nkx x Nky
    szU = size(U);
    if numel(szU) == 3
        nb = szU(1); nb_sel = szU(2);
        if szU(3) ~= Nkx*Nky
            error('If U is 3D, its 3rd dim must be Nkx*Nky.');
        end
        U = reshape(U, [nb, nb_sel, Nkx, Nky]);
    elseif numel(szU) == 4
        nb = szU(1); nb_sel = szU(2);
        if szU(3)~=Nkx || szU(4)~=Nky
            error('U must be nb x nb_sel x Nkx x Nky.');
        end
    else
        error('U must be 4D (nb x nb_sel x Nkx x Nky) or 3D (nb x nb_sel x Nkx*Nky).');
    end

    % reshape E to Nkx x Nky x nb_sel
    if ~isequal(size(E), [Nkx,Nky,nb_sel])
        if isequal(size(E), [Nkx*Nky, nb_sel])
            E = reshape(E, [Nkx, Nky, nb_sel]);
        else
            error('E must be Nkx x Nky x nb_sel (or (Nkx*Nky) x nb_sel).');
        end
    end

    Eph_list = Eph_list(:).'; % 1 x Nw
    Nw = numel(Eph_list);

    % -------------------- skew grid Jacobian --------------------
    [dk1, dk2] = grid_step_vectors(Kx, Ky); % 2x1
    J = [dk1(:), dk2(:)];                   % 2x2
    invJT = inv(J.');                       % J^{-T}

    % k-space weight: d^2k/(2π)^2, cell area = |det(J)|
    d2k = abs(det(J));
    w_k = d2k / (2*pi)^2;

    if verbose
        fprintf('[MIC metric mn] Nkx=%d Nky=%d nb=%d nb_sel=%d\n', Nkx, Nky, nb, nb_sel);
        fprintf('[MIC metric mn] dk1=(%.3e,%.3e), dk2=(%.3e,%.3e), det(J)=%.3e\n', dk1(1),dk1(2),dk2(1),dk2(2), det(J));
        fprintf('[MIC metric mn] periodicFD=%d, trimBoundary=%d, doGaugeFix=%d, g_s=%g, positiveDE=%d\n', ...
            periodicFD, trimBoundary, doGaugeFix, g_s, positiveDE);
    end

    % -------------------- gauge smoothing --------------------
    if doGaugeFix
        Ug = gauge_fix_parallel_transport(U);
    else
        Ug = U;
    end

    % -------------------- k-mask (interior only if trimBoundary & nonperiodic) ----
    mask_k = ones(Nkx, Nky);
    if trimBoundary && ~periodicFD
        mask_k(:,:) = 0;
        mask_k(2:Nkx-1, 2:Nky-1) = 1;
    end
    mask_k4 = reshape(mask_k, [Nkx, Nky, 1, 1]); % broadcast to (Nkx,Nky,nb,nb)

    % -------------------- eigenvector derivatives ∂u/∂kx, ∂u/∂ky --------------------
    [du_x, du_y] = fd_du_skew(Ug, invJT, periodicFD);

    % -------------------- r_mn(k) for all m,n --------------------
    % r_mn: Nkx x Nky x nb_sel x nb_sel x 2
    r_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 'like', Ug);

    ix_list = 1:Nkx;
    iy_list = 1:Nky;
    if trimBoundary && ~periodicFD
        ix_list = 2:(Nkx-1);
        iy_list = 2:(Nky-1);
    end

    for ix = ix_list
        for iy = iy_list
            U0  = squeeze(Ug(:,:,ix,iy));    % nb x nb_sel
            dux = squeeze(du_x(:,:,ix,iy));  % nb x nb_sel
            duy = squeeze(du_y(:,:,ix,iy));  % nb x nb_sel

            Mx = U0' * dux; % <u_m|∂x u_n>
            My = U0' * duy; % <u_m|∂y u_n>

            r_mn(ix,iy,:,:,1) = 1i * Mx;
            r_mn(ix,iy,:,:,2) = 1i * My;
        end
    end

    % enforce off-diagonal only using a mask (branchless)
    off = ones(nb_sel) - eye(nb_sel);
    mask_off = reshape(off, [1,1,nb_sel,nb_sel]); % broadcast

    % -------------------- build r_nm via conjugate transpose in band indices ----
    % r_nm(k,n,m,*) = conj(r_mn(k,m,n,*))
    % r_nm = permute(conj(r_mn), [1 2 4 3 5]); % Nkx x Nky x nb_sel(n) x nb_sel(m) x 2
    r_mn_H = 0.5 * ( r_mn + permute(conj(r_mn), [1 2 4 3 5]) );
    % For convenience, we want r_nm with indices (m,n) again:
    % r_nm_mn(k,m,n,*) = r_nm(k,n,m,*) => permute back:
    % r_nm_mn = permute(r_nm, [1 2 4 3 5]);     % Nkx x Nky x nb_sel(m) x nb_sel(n) x 2
    r_nm_mn = permute(r_mn_H, [1 2 4 3 5]);   % Nkx x Nky x m x n x 2


    % -------------------- metric kernel: g2^{bc}_{mn} (gauge invariant) ----
    % g2: Nkx x Nky x nb_sel x nb_sel x 2(b) x 2(c)
    g2 = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 2, 'like', Ug);
    for b = 1:2
        for cc = 1:2
            % 2 g^{bc} = r^b_mn r^c_nm + r^c_mn r^b_nm
            % g2(:,:,:,:,b,cc) = real( ...
            %     r_mn(:,:,:,:,b) .* r_nm_mn(:,:,:,:,cc) + ...
            %     r_mn(:,:,:,:,cc) .* r_nm_mn(:,:,:,:,b) );
            g2(:,:,:,:,b,cc) = real( ...
                r_mn_H(:,:,:,:,b)  .* r_nm_mn(:,:,:,:,cc) + ...
                r_mn_H(:,:,:,:,cc) .* r_nm_mn(:,:,:,:,b)  );
        end
    end

    % -------------------- energies: build ΔE_mn = E_m - E_n --------------------
    Em = reshape(E, [Nkx, Nky, nb_sel, 1]); % ... m
    En = reshape(E, [Nkx, Nky, 1, nb_sel]); % ... n
    dE_mn = Em - En;                        % Nkx x Nky x m x n (eV)

    % absorption-only mask ΔE>0 (branchless)
    if positiveDE
        mask_pos = double(dE_mn > 0);
    else
        mask_pos = ones(size(dE_mn));
    end

    % -------------------- Fermi factor f_nm = f_n - f_m --------------------
    f = fermi_dirac(E, Ef, kT); % Nkx x Nky x nb_sel
    fm = reshape(f, [Nkx, Nky, nb_sel, 1]); % final m
    fn = reshape(f, [Nkx, Nky, 1, nb_sel]); % initial n
    f_nm = fn - fm;                         % Nkx x Nky x m x n

    % -------------------- velocities v_n = (e/ħ) ∂E/∂k (multi-band) --------
    [dE_dx, dE_dy] = fd_tensor_skew(E, invJT, periodicFD); % Nkx x Nky x nb_sel
    v_band = zeros(Nkx, Nky, nb_sel, 2);
    v_band(:,:,:,1) = vel_pref * dE_dx;  % vx
    v_band(:,:,:,2) = vel_pref * dE_dy;  % vy

    % Δv^a_mn = v_m^a - v_n^a  (for a=x,y)
    dv_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2);
    for a = 1:2
        Vm = reshape(v_band(:,:,:,a), [Nkx, Nky, nb_sel, 1]);
        Vn = reshape(v_band(:,:,:,a), [Nkx, Nky, 1, nb_sel]);
        dv_mn(:,:,:,:,a) = Vm - Vn;
    end

    % -------------------- total (k,m,n) mask --------------------
    % use: interior k (mask_k4), off-diagonal (mask_off), positive ΔE (mask_pos)
    mask_total = mask_k4 .* mask_off .* mask_pos; % Nkx x Nky x m x n

    % -------------------- integrate over photon energy --------------------
    eta_abc = zeros(2,2,2,Nw);

    % Precompute the (k,m,n)-resolved part except delta(Eph):
    % For each (a,b,c):
    %   Kernel^{abc}_{mn}(k) = f_nm * Δv^a_mn * g2^{bc}_mn * mask_total
    %
    % Then eta^{abc}(Eph) = pref * w_k * Σ Kernel * delta_w
    for a = 1:2
        dvA = dv_mn(:,:,:,:,a); % Nkx x Nky x m x n
        for b = 1:2
            for cc = 1:2
                gBC = g2(:,:,:,:,b,cc); % Nkx x Nky x m x n
                Kabc = (f_nm .* dvA .* gBC) .* mask_total;

                % loop over photon energies (avoid huge 5D allocation)
                for iw = 1:Nw
                    Eph = Eph_list(iw);
                    delta_w = (1/pi) * eta ./ ((dE_mn - Eph).^2 + eta^2); % Nkx x Nky x m x n, unit 1/eV
                    eta_abc(a,b,cc,iw) = pref * w_k * sum(Kabc .* delta_w, 'all');
                end
            end
        end
    end

    % -------------------- outputs --------------------
    out = struct();
    out.pref      = pref;
    out.vel_pref  = vel_pref;
    out.w_k       = w_k;
    out.J         = J;
    out.invJT     = invJT;
    out.dk1       = dk1;
    out.dk2       = dk2;
    out.mask_k    = mask_k;
    out.mask_off  = off;
    out.positiveDE = positiveDE;

    out.Ug        = Ug;

    if saveFullMN
        out.r_mn    = r_mn;     % length_unit
        out.g2      = g2;       % length_unit^2, gauge-invariant
        out.dE_mn   = dE_mn;    % eV
        out.f_nm    = f_nm;     % dimensionless
        out.v_band  = v_band;   % length_unit/s
        out.dv_mn   = dv_mn;    % length_unit/s
    end
end

% ======================================================================
% helpers
% ======================================================================

function val = get_opt(opts, name, default)
    if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = default;
    end
end

function [dk1, dk2] = grid_step_vectors(Kx, Ky)
% infer two step vectors in Cartesian k-space directly from Kx,Ky arrays
    [Nkx,Nky] = size(Kx);

    dKx1 = diff(Kx(:,1)); dKy1 = diff(Ky(:,1));
    idx1 = find((abs(dKx1)+abs(dKy1))>0, 1, 'first');
    if isempty(idx1) || Nkx < 2
        error('Cannot infer dk1 from K-grid.');
    end
    dk1 = [dKx1(idx1); dKy1(idx1)];

    dKx2 = diff(Kx(1,:)); dKy2 = diff(Ky(1,:));
    idx2 = find((abs(dKx2)+abs(dKy2))>0, 1, 'first');
    if isempty(idx2) || Nky < 2
        error('Cannot infer dk2 from K-grid.');
    end
    dk2 = [dKx2(idx2); dKy2(idx2)];
end

function Ug = gauge_fix_parallel_transport(U)
% band-wise phase smoothing (nondegenerate)
    Ug = U;
    [~, nb_sel, Nkx, Nky] = size(U);

    for iy = 1:Nky
        for ix = 2:Nkx
            Uprev = squeeze(Ug(:,:,ix-1,iy));
            Ucur  = squeeze(Ug(:,:,ix,iy));
            for n = 1:nb_sel
                ov = Uprev(:,n)' * Ucur(:,n);
                ph = ov / max(abs(ov), 1e-30);
                Ucur(:,n) = Ucur(:,n) / ph;
            end
            Ug(:,:,ix,iy) = Ucur;
        end
    end

    for ix = 1:Nkx
        for iy = 2:Nky
            Uprev = squeeze(Ug(:,:,ix,iy-1));
            Ucur  = squeeze(Ug(:,:,ix,iy));
            for n = 1:nb_sel
                ov = Uprev(:,n)' * Ucur(:,n);
                ph = ov / max(abs(ov), 1e-30);
                Ucur(:,n) = Ucur(:,n) / ph;
            end
            Ug(:,:,ix,iy) = Ucur;
        end
    end
end

function [du_x, du_y] = fd_du_skew(U, invJT, periodicFD)
% Cartesian derivatives of eigenvectors on a skew grid
    [nb, nb_sel, Nkx, Nky] = size(U);
    du_x = zeros(nb, nb_sel, Nkx, Nky, 'like', U);
    du_y = zeros(nb, nb_sel, Nkx, Nky, 'like', U);

    if periodicFD
        idxp1 = [2:Nkx, 1];
        idxm1 = [Nkx, 1:Nkx-1];
        idxp2 = [2:Nky, 1];
        idxm2 = [Nky, 1:Nky-1];

        du_k1 = (U(:,:,idxp1,:) - U(:,:,idxm1,:)) / 2;
        du_k2 = (U(:,:,:,idxp2) - U(:,:,:,idxm2)) / 2;
    else
        du_k1 = zeros(size(U), 'like', U);
        du_k2 = zeros(size(U), 'like', U);

        ix = 2:(Nkx-1);
        iy = 2:(Nky-1);

        du_k1(:,:,ix,:) = (U(:,:,ix+1,:) - U(:,:,ix-1,:)) / 2;
        du_k2(:,:,:,iy) = (U(:,:,:,iy+1) - U(:,:,:,iy-1)) / 2;
    end

    du_x = invJT(1,1)*du_k1 + invJT(1,2)*du_k2;
    du_y = invJT(2,1)*du_k1 + invJT(2,2)*du_k2;
end

function [dF_dx, dF_dy] = fd_tensor_skew(F, invJT, periodicFD)
%FD_TENSOR_SKEW
% Central differences for a tensor field F(ix,iy,band) on skew grid,
% then transform to Cartesian derivatives using inv(J').
%
% Input:
%   F: Nkx x Nky x Nb
% Output:
%   dF_dx, dF_dy: Nkx x Nky x Nb

    [Nkx,Nky,Nb] = size(F);
    dF_dx = zeros(Nkx,Nky,Nb);
    dF_dy = zeros(Nkx,Nky,Nb);

    if periodicFD
        idxp1 = [2:Nkx, 1];
        idxm1 = [Nkx, 1:Nkx-1];
        idxp2 = [2:Nky, 1];
        idxm2 = [Nky, 1:Nky-1];

        d_k1 = (F(idxp1,:,:) - F(idxm1,:,:)) / 2;
        d_k2 = (F(:,idxp2,:) - F(:,idxm2,:)) / 2;
    else
        d_k1 = zeros(Nkx,Nky,Nb);
        d_k2 = zeros(Nkx,Nky,Nb);

        ix = 2:(Nkx-1);
        iy = 2:(Nky-1);

        d_k1(ix,:,:) = (F(ix+1,:,:) - F(ix-1,:,:)) / 2;
        d_k2(:,iy,:) = (F(:,iy+1,:) - F(:,iy-1,:)) / 2;
    end

    dF_dx = invJT(1,1)*d_k1 + invJT(1,2)*d_k2;
    dF_dy = invJT(2,1)*d_k1 + invJT(2,2)*d_k2;
end

function f = fermi_dirac(E, Ef, kT)
% energies in eV
    if kT <= 0
        f = double(E < Ef);
    else
        f = 1 ./ (1 + exp((E - Ef)/kT));
    end
end


function ef=get_ef(g,pars,model,nbands,knum)
    kxline=[-0.1,0.1];
    kyline=[-0.1,0.1];
    u=0.5;
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_kp(pars,model,nbands,Kx,Ky,Kz);
    ef=calculate_ef(Enk(:), u);
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end
