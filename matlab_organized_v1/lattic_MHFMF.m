% clear;
% clear all;
g = MTB.geometry("Rgra_5s");
g = MTB.read_poscar(g,"data/Graphene/5s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/5s/wannier90_hr_p1.dat','data/Graphene/5s/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%%
% 假设 g 已经在工作区（来自 Wannier TB）
opts = struct;
opts.Nk = [101 101 1];
opts.nspin = 1;
opts.kT = 0;
opts.vmodel = '2D';     % 或 'keldysh'
opts.eps    = 20.0;
opts.r0     = 10.0;     % 仅 keldysh 用
opts.max_iter = 60;
opts.tol      = 1e-6;
opts.mixing   = 0.7;
opts.use_wpos_cart = true;      % 若 g.wpos 是分数坐标，把它设 false
opts.use_tau_phase = false;      % S_ab(q) 形状因子
opts.fock_drop_q0  = true;
opts.hartree_enable= true;

% 可选：设置参考电荷（例如中性/初态层密度）
% opts.nref = your_reference_vector(:);

out = hf_run_wannier(g, opts);
%%
% 查看迭代收敛
semilogy(out.history,'-o'), grid on
xlabel('SCF iter'), ylabel('||\Delta\rho||/||\rho||')
%%
figure()
for j=1:101
for i=1:20
    hold on;
plot(squeeze(out.E_k(j,:,1,i)))
end
end


%%
function out = hf_run_wannier(g, opts)
% HF self-consistent calculation in Wannier/orbital basis (uniform phase).
% Atomic gauge by default (no intra-cell phase in density vertex).
%
% INPUT:
%   g: struct with fields
%      .a      (3x3) real-space lattice vectors (Angstrom)
%      .b      (3x3) reciprocal lattice vectors (1/Angstrom), b·a = 2*pi*I
%      .hopr   (Lx3) hopping lattice vectors in fractional coords
%      .ham    (m x m x L) hopping matrices H(R)
%      .wpos   (m x 3) orbital positions; set opts.use_wpos_cart accordingly
%   opts: struct (key fields; defaults below)
%      .Nk           [Nk1 Nk2 Nk3]  k-mesh
%      .Efilling     total electrons per cell (including spin)
%      .nspin        1 if spin is explicit in g.ham; 2 if spin-degenerate bands
%      .kT           temperature in eV (Fermi smearing)
%      .max_iter     SCF max iterations
%      .tol          SCF tolerance on ||Δρ||/||ρ||
%      .mixing       linear mixing (0..1)
%      .vmodel       '2D' or 'keldysh'
%      .eps          dielectric constant (for '2D')
%      .e2           e^2 (in eV·Angstrom), default 14.3996
%      .r0           Keldysh length (Angstrom), if vmodel='keldysh'
%      .use_wpos_cart   (true) wpos is Cartesian; false if fractional
%      .use_tau_phase   (false) atomic gauge => no exp(-iq·(tau_b - tau_a))
%      .hartree_enable  (true) include Hartree (uniform, layer-resolved)
%      .verbose         (true) print SCF diagnostics each iteration
%      .nref         (m x 1) reference charges per orbital (neutral background)
%
% OUTPUT:
%   out: struct with fields
%      .rho_k, .E_k, .U_k, .SigmaF_k, .SigmaH_k, .H0_k, .HMF_k, .mu, .history, .log

% -------------------------- defaults --------------------------
if ~isfield(opts, 'Nk'),              opts.Nk = [24 24 1];      end
if ~isfield(opts, 'nspin'),           opts.nspin = 2;           end
if ~isfield(opts, 'kT'),              opts.kT = 1e-4;           end  % eV
if ~isfield(opts, 'max_iter'),        opts.max_iter = 100;      end
if ~isfield(opts, 'tol'),             opts.tol = 1e-7;          end
if ~isfield(opts, 'mixing'),          opts.mixing = 0.7;        end
if ~isfield(opts, 'vmodel'),          opts.vmodel = '2D';       end
if ~isfield(opts, 'eps'),             opts.eps = 4.0;           end
if ~isfield(opts, 'e2'),              opts.e2 = 14.3996;        end % eV·Å
if ~isfield(opts, 'r0'),              opts.r0 = 10.0;           end
if ~isfield(opts, 'use_wpos_cart'),   opts.use_wpos_cart = true;end
if ~isfield(opts, 'use_tau_phase'),   opts.use_tau_phase = false;end % atomic gauge
if ~isfield(opts, 'hartree_enable'),  opts.hartree_enable = true;end
if ~isfield(opts, 'verbose'),         opts.verbose = true;      end

% -------------------------- sizes -----------------------------
m     = size(g.ham,1);                 % # orbitals
Nk1   = opts.Nk(1);  Nk2 = opts.Nk(2); Nk3 = opts.Nk(3);
NkTot = Nk1*Nk2*Nk3;

% -------------------- orbital positions ----------------------
tau = g.wpos;                          % Å if use_wpos_cart
if ~opts.use_wpos_cart
    tau = (g.a * tau.').';             % fractional -> Cartesian (Å)
end

% -------------------------- k-mesh ----------------------------
[kgrid, kfrac, kind] = make_kmesh(g.b, opts.Nk); % kgrid: Nk1xNk2xNk3x3 (cart)

% -------------------------- H0(k) -----------------------------
H0_k = zeros(Nk1,Nk2,Nk3,m,m);
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    kvec = squeeze(kgrid(i1,i2,i3,:)).';
    H0_k(i1,i2,i3,:,:) = H0k_from_wannier(g, kvec);
end

% -------- target electrons per cell (including spin) ----------
if ~isfield(opts,'Efilling') || isempty(opts.Efilling)
    % default: half-filling per spin => total ~ m * nspin / 1? Adjust to your case.
    opts.Efilling = m * opts.nspin / 2; % change if you need a different default
end
Efilling = opts.Efilling;

% --------------- initial diagonalization, rho -----------------
E_k   = zeros(Nk1,Nk2,Nk3,m);
U_k   = cell(NkTot,1);
rho_k = zeros(Nk1,Nk2,Nk3,m,m);
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    Hk = squeeze(H0_k(i1,i2,i3,:,:));
    [V,D] = eig((Hk+Hk')/2);
    Ek = real(diag(D));
    E_k(i1,i2,i3,:) = Ek;
    U_k{ik} = V;
end
% chemical potential (non-interacting)
mu = find_mu(E_k, Efilling, opts.kT, opts.nspin);
% build rho from non-interacting bands
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    Ek = squeeze(E_k(i1,i2,i3,:)).';
    fk = fermi_dirac(Ek, mu, opts.kT);   % per spin
    Uk = U_k{ik};
    rho_spin = Uk * diag(fk) * Uk';      % per spin density
    rho_k(i1,i2,i3,:,:) = opts.nspin * rho_spin;   % total density (incl spin)
end

% ---------------- precompute V(q) for Fock --------------------
Vq = build_Vq(g, tau, opts, kfrac);     % (Nk1 x Nk2 x Nk3 x m x m)

% -------------------------- SCF loop --------------------------
history = zeros(opts.max_iter,1);
logmat  = nan(opts.max_iter, 5); % [iter, nr, max|dρ|, μ, Etot]
mu_prev   = mu;
Etot_prev = NaN;

for it = 1:opts.max_iter
    % Fock: Σ^F(k) via FFT-based convolution
    SigmaF_k = fock_from_fft(Vq, rho_k, opts);

    % Hartree (uniform, layer-resolved, diagonal, k-independent)
    if opts.hartree_enable
        SigmaH_k = hartree_uniform_smallq(g, tau, rho_k, opts);
    else
        SigmaH_k = zeros(Nk1,Nk2,Nk3,m,m);
    end

    % Mean-field Hamiltonian
    HMF_k = H0_k + SigmaH_k + SigmaF_k;

    % Diagonalize H_MF, update rho_new
    E_k_new = zeros(Nk1,Nk2,Nk3,m);
    U_k_new = cell(NkTot,1);
    for ik = 1:NkTot
        [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
        Hk = squeeze(HMF_k(i1,i2,i3,:,:));
        [V,D] = eig((Hk+Hk')/2);
        Ek = real(diag(D));
        E_k_new(i1,i2,i3,:) = Ek;
        U_k_new{ik} = V;
    end
    % Update chemical potential for interacting bands
    mu = find_mu(E_k_new, Efilling, opts.kT, opts.nspin);

    rho_new = zeros(Nk1,Nk2,Nk3,m,m);
    for ik = 1:NkTot
        [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
        Ek = squeeze(E_k_new(i1,i2,i3,:)).';
        fk = fermi_dirac(Ek, mu, opts.kT);  % per spin
        Uk = U_k_new{ik};
        rho_spin = Uk * diag(fk) * Uk';
        rho_new(i1,i2,i3,:,:) = opts.nspin * rho_spin;
    end

    % -------- diagnostics BEFORE mixing ----------
    drho   = rho_new - rho_k;
    nr     = norm(drho(:)) / max(1, norm(rho_k(:)));  % relative Frobenius
    nr_inf = max(abs(drho(:)));                       % infinity-norm
    Etot   = hf_total_energy(H0_k, SigmaH_k, SigmaF_k, rho_new);

    if opts.verbose
        if isnan(Etot_prev)
            fprintf('SCF %3d | nr=%.3e | max|dρ|=%.3e | Δμ=% .3e eV | E=% .10f eV\n',...
                it, nr, nr_inf, (mu-mu_prev), Etot);
        else
            fprintf('SCF %3d | nr=%.3e | max|dρ|=%.3e | Δμ=% .3e eV | E=% .10f eV (ΔE=% .3e)\n',...
                it, nr, nr_inf, (mu-mu_prev), Etot, (Etot-Etot_prev));
        end
    end
    logmat(it,:) = [it, nr, nr_inf, mu, Etot];
    mu_prev   = mu;
    Etot_prev = Etot;

    % ---------------- mixing & convergence ----------------
    rho_k = (1-opts.mixing)*rho_k + opts.mixing*rho_new;
    E_k   = E_k_new;  U_k = U_k_new;

    history(it) = nr;
    if nr < opts.tol
        history = history(1:it);
        logmat  = logmat(1:it,:);
        if opts.verbose
            fprintf('HF converged at iter %d: nr=%.3e\n', it, nr);
        end
        break;
    end
    if it == opts.max_iter && opts.verbose
        fprintf('HF stopped at max_iter %d: nr=%.3e\n', it, nr);
        logmat = logmat(1:it,:);
    end
end

% --------------------------- output ---------------------------
out.rho_k    = rho_k;
out.E_k      = E_k;
out.U_k      = U_k;
out.SigmaF_k = SigmaF_k;
out.SigmaH_k = SigmaH_k;
out.H0_k     = H0_k;
out.HMF_k    = HMF_k;
out.mu       = mu;
out.history  = history;
out.log      = logmat;
out.kgrid    = kgrid;   % Cartesian k (Å^-1)
out.kfrac    = kfrac;   % fractional k in [0,1)
out.kind     = kind;

end

%%
function out = hf_run_wannier_old(g, opts)
% HF self-consistent calculation in Wannier/orbital basis (uniform phase).
%
% INPUT:
%   g: struct with fields
%      .a      (3x3) real-space lattice vectors (Angstrom)
%      .b      (3x3) reciprocal lattice vectors (1/Angstrom), b·a = 2*pi*I
%      .hopr   (Lx3) hopping lattice vectors in fractional coords
%      .ham    (m x m x L) hopping matrices H(R)
%      .wpos   (m x 3) orbital positions; set opts.use_wpos_cart accordingly
%   opts: struct (key fields listed below; all have sensible defaults)
%      .Nk         [Nk1 Nk2 Nk3]      k-mesh
%      .Efilling   scalar OR vector  target electrons per cell (total, incl spin)
%      .nspin      1 or 2 (collinear, no-SOC)
%      .kT         temperature in eV (Fermi smearing)
%      .max_iter   SCF max iterations
%      .tol        SCF tolerance on density matrix
%      .mixing     linear mixing for rho (0~1)
%      .vmodel     '2D' or 'keldysh'
%      .eps        dielectric constant (for '2D')
%      .e2         e^2 (in eV·Angstrom), default 14.3996 (1/(4pi eps0))
%      .r0         Keldysh length (Angstrom), if vmodel='keldysh'
%      .use_wpos_cart  (true) wpos is Cartesian; false if fractional
%      .use_tau_phase  (true) include S_ab(q)=exp(-iq·(tau_b-tau_a)) in Vq
%      .fock_drop_q0   (true) set V(q=0)=0 in Fock to avoid self-interaction
%      .hartree_enable (true) include Hartree (uniform)
%      .nref       (m x 1) reference charges per orbital (neutral background)
%
% OUTPUT:
%   out: struct with fields
%      .rho_k      (Nk1 x Nk2 x Nk3 x m x m) density matrix in k
%      .E_k        (Nk1 x Nk2 x Nk3 x m)     band energies (eV)
%      .U_k        cell{NkTot} of (m x m) eigenvectors (column unitary)
%      .SigmaF_k   (Nk1 x Nk2 x Nk3 x m x m) Fock self-energy
%      .SigmaH_k   (Nk1 x Nk2 x Nk3 x m x m) Hartree self-energy (diag, k-indep)
%      .H0_k       (Nk1 x Nk2 x Nk3 x m x m) TB Hamiltonian
%      .HMF_k      (Nk1 x Nk2 x Nk3 x m x m) mean-field Hamiltonian
%      .mu         final chemical potential
%      .history    SCF history (norms)
%
% -------------------------------------------------------------------------

% ---------- defaults ----------
if ~isfield(opts, 'Nk'),              opts.Nk = [24 24 1];      end
if ~isfield(opts, 'nspin'),           opts.nspin = 2;           end
if ~isfield(opts, 'kT'),              opts.kT = 1e-4;           end  % eV
if ~isfield(opts, 'max_iter'),        opts.max_iter = 100;      end
if ~isfield(opts, 'tol'),             opts.tol = 1e-7;          end
if ~isfield(opts, 'mixing'),          opts.mixing = 0.7;        end
if ~isfield(opts, 'vmodel'),          opts.vmodel = '2D';       end
if ~isfield(opts, 'eps'),             opts.eps = 4.0;           end
if ~isfield(opts, 'e2'),              opts.e2 = 14.3996;        end % eV*Ang
if ~isfield(opts, 'r0'),              opts.r0 = 10.0;           end
if ~isfield(opts, 'use_wpos_cart'),   opts.use_wpos_cart = true;end
if ~isfield(opts, 'use_tau_phase'),   opts.use_tau_phase = true;end
if ~isfield(opts, 'fock_drop_q0'),    opts.fock_drop_q0 = true; end
if ~isfield(opts, 'hartree_enable'),  opts.hartree_enable = true;end

% ---------- sizes ----------
m    = size(g.ham,1);                 % number of orbitals
Nk1  = opts.Nk(1);  Nk2 = opts.Nk(2); Nk3 = opts.Nk(3);
NkTot = Nk1*Nk2*Nk3;

% ---------- orbital positions (Cartesian) ----------
tau = g.wpos;
if ~opts.use_wpos_cart
    tau = (g.a * tau.').';    % fractional -> Cartesian
end

% ---------- k-mesh ----------
[kgrid, kfrac, kind] = make_kmesh(g.b, opts.Nk); % kgrid: Nk1xNk2xNk3x3 (cart)

% ---------- build H0(k) ----------
H0_k = zeros(Nk1,Nk2,Nk3,m,m);
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    kvec = squeeze(kgrid(i1,i2,i3,:)).';
    H0_k(i1,i2,i3,:,:) = H0k_from_wannier(g, kvec);
end

% ---------- chemical potential target electrons ----------
if ~isfield(opts,'Efilling') || isempty(opts.Efilling)
    % 默认填充 = m 电子/自旋 * nspin（半充满可改）
    Efilling = m * opts.nspin;  % 每胞电子数（总）
else
    Efilling = opts.Efilling;
end

% ---------- initial diagonalization, density ----------
E_k  = zeros(Nk1,Nk2,Nk3,m);
U_k  = cell(NkTot,1);
rho_k = zeros(Nk1,Nk2,Nk3,m,m);
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    Hk = squeeze(H0_k(i1,i2,i3,:,:));
    [V,D] = eig((Hk+Hk')/2);    % Hermitian safety
    Ek = real(diag(D));
    E_k(i1,i2,i3,:) = Ek;
    U_k{ik} = V;
end
% mu for non-interacting filling
mu = find_mu(E_k, Efilling, opts.kT, opts.nspin);
% build rho
for ik = 1:NkTot
    [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
    Ek = squeeze(E_k(i1,i2,i3,:)).';
    fk = fermi_dirac(Ek, mu, opts.kT);   % per spin
    % 每自旋占据矩阵 UFU^\dag
    Uk = U_k{ik};
    rho_spin = Uk * diag(fk) * Uk';
    rho_k(i1,i2,i3,:,:) = opts.nspin * rho_spin;   % Hartree用总密度
end

% ---------- precompute V(q) for Fock ----------
Vq = build_Vq(g, tau, opts, kfrac);   % (Nk1 x Nk2 x Nk3 x m x m)

% ---------- SCF loop ----------
history = zeros(opts.max_iter,1);
for it = 1:opts.max_iter
    % Fock self-energy: Σ^F(k) = - FFT^{-1}[ V(r) .* ρ(r) ]
    SigmaF_k = fock_from_fft(Vq, rho_k, opts);
    
    % Hartree (uniform, layer-resolved, diag)
    if opts.hartree_enable
        SigmaH_k = hartree_uniform_smallq(g, tau, rho_k, opts);
    else
        SigmaH_k = zeros(Nk1,Nk2,Nk3,m,m);
    end
    
    % H_MF = H0 + ΣH + ΣF
    HMF_k = H0_k + SigmaH_k + SigmaF_k;
    
    % diagonalize H_MF, update rho
    E_k_new  = zeros(Nk1,Nk2,Nk3,m);
    U_k_new  = cell(NkTot,1);
    rho_new  = zeros(Nk1,Nk2,Nk3,m,m);
    for ik = 1:NkTot
        [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
        Hk = squeeze(HMF_k(i1,i2,i3,:,:));
        [V,D] = eig((Hk+Hk')/2);
        Ek = real(diag(D));
        E_k_new(i1,i2,i3,:) = Ek;
        U_k_new{ik} = V;
    end
    % 更新化学势以满足总电子数
    mu = find_mu(E_k_new, Efilling, opts.kT, opts.nspin);
    for ik = 1:NkTot
        [i1,i2,i3] = ind2sub([Nk1 Nk2 Nk3], ik);
        Ek = squeeze(E_k_new(i1,i2,i3,:)).';
        fk = fermi_dirac(Ek, mu, opts.kT);
        Uk = U_k_new{ik};
        rho_spin = Uk * diag(fk) * Uk';      % per spin
        rho_new(i1,i2,i3,:,:) = opts.nspin * rho_spin;
    end
    
    % mixing & check convergence
    drho = rho_new - rho_k;
    nr = norm(drho(:)) / max(1, norm(rho_k(:)));
    history(it) = nr;
    rho_k = (1-opts.mixing)*rho_k + opts.mixing*rho_new;
    E_k  = E_k_new; U_k = U_k_new;
    if nr < opts.tol
        history = history(1:it);
        break;
    end
end

% 输出
out.rho_k    = rho_k;
out.E_k      = E_k;
out.U_k      = U_k;
out.SigmaF_k = SigmaF_k;
out.SigmaH_k = SigmaH_k;
out.H0_k     = H0_k;
out.HMF_k    = HMF_k;
out.mu       = mu;
out.history  = history;
end

function Hk = H0k_from_wannier(g, kvec)
% H(k) = sum_R H(R) e^{i k·R} ; g.hopr in fractional coords.

m = size(g.ham,1);
Hk = zeros(m,m);
L  = size(g.hopr,1);
for l = 1:L
    Rfrac = g.hopr(l,:);         % fractional
    Rvec  = (g.a * Rfrac.').';   % Cartesian (Angstrom)
    phase = exp(1i * dot(kvec, Rvec));
    Hk = Hk + phase * g.ham(:,:,l);
end
% Hermitian safety
Hk = (Hk + Hk')/2;
end


function [kgrid, kfrac, kind] = make_kmesh(B, Nk)
% B: (3x3) reciprocal lattice basis (1/Angstrom), row-vectors or columns
% Nk: [Nk1 Nk2 Nk3]
% Output:
%   kgrid (Nk1 x Nk2 x Nk3 x 3): Cartesian k
%   kfrac (Nk1 x Nk2 x Nk3 x 3): fractional in [0,1)
%   kind  linear indices

Nk1 = Nk(1); Nk2 = Nk(2); Nk3 = Nk(3);
[k1, k2, k3] = ndgrid(0:Nk1-1, 0:Nk2-1, 0:Nk3-1);
kfrac = zeros(Nk1,Nk2,Nk3,3);
kfrac(:,:,:,1) = k1 / Nk1;
kfrac(:,:,:,2) = k2 / Nk2;
kfrac(:,:,:,3) = k3 / Nk3;

kgrid = zeros(Nk1,Nk2,Nk3,3);
for i=1:Nk1
  for j=1:Nk2
    for k=1:Nk3
      kgrid(i,j,k,:) = (B * squeeze(kfrac(i,j,k,:))).';
    end
  end
end
kind = reshape(1:(Nk1*Nk2*Nk3), [Nk1 Nk2 Nk3]);
end


function f = fermi_dirac(E, mu, kT)
% E: 1xM energies (eV)
% mu, kT: eV
if kT<=0
    f = double(E<mu);
else
    x = (E - mu)/kT;
    % 防溢出
    x = max(min(x, 40), -40);
    f = 1 ./ (1 + exp(x));
end
end

function mu = find_mu(E_k, Nelec, kT, nspin)
% total electrons per cell (including spin) = Nelec
% E_k: (Nk1 x Nk2 x Nk3 x m)

E = sort(E_k(:));
Emin = E(1) - 5*max(kT,1e-4);
Emax = E(end) + 5*max(kT,1e-4);
NkTot = numel(E_k)/numel(E);

m = size(E_k,4);

% bisection
for it=1:80
    mu  = 0.5*(Emin+Emax);
    fsum = 0;
    % integrate over BZ: sum_k sum_n f(Ekn)
    fsum = sum( fermi_dirac(E_k(:), mu, kT) );
    % 这是"每自旋"的求和；实际每 k 点/每能带只有一次计数，
    % 所以总电子数 = nspin * fsum / NkTot
    Ne_mu = nspin * fsum / NkTot;
    if Ne_mu > Nelec
        Emax = mu;
    else
        Emin = mu;
    end
    if abs(Ne_mu - Nelec) < 1e-10*max(1,Nelec)
        break;
    end
end
end


function Vq = build_Vq(g, tau, opts, kfrac)
% Vq: (Nk1 x Nk2 x Nk3 x m x m)

[Nk1,Nk2,Nk3,~] = size(kfrac);
m = size(tau,1);
Vq = zeros(Nk1,Nk2,Nk3,m,m);

% pre z-positions for layer screening
% 若层法向不是 z，请改成 tau*nhat
z = tau(:,3);

for i1 = 1:Nk1
  for i2 = 1:Nk2
    for i3 = 1:Nk3
      qfrac = squeeze(kfrac(i1,i2,i3,:)).';  % fractional
      qvec  = (g.b * qfrac.').';             % Cartesian (1/Ang)
      qmag  = norm(qvec);
      if qmag < 1e-14
          % q=0 点在 Fock 里通常置 0 防止自作用；Hartree 另行处理
          Vscalar = 0.0;
      else
          switch lower(opts.vmodel)
              case '2d'
                  Vscalar = 2*pi*opts.e2/(opts.eps*qmag);   % eV
              case 'keldysh'
                  % 常用形式，e^2/[(eps1+eps2)2eps0] 的整体系数
                  % 以 eV*Angstrom 记法近似为 opts.e2；实现上仍保持 ~1/q(1+r0 q)
                  Vscalar = opts.e2 / (qmag*(1+opts.r0*qmag));
              otherwise
                  error('Unknown vmodel');
          end
      end
      
      % orbital-by-orbital kernel
      for a = 1:m
        for b = 1:m
            layer_scr = exp(-qmag * abs(z(a)-z(b)));  % 2D electrons + interlayer distance
            if opts.use_tau_phase
                phase = exp(-1i * dot(qvec, tau(b,:) - tau(a,:)));
            else
                phase = 1.0;
            end
            Vq(i1,i2,i3,a,b) = Vscalar * layer_scr * phase;
        end
      end
    end
  end
end
end


function SigmaF_k = fock_from_fft(Vq, rho_k, opts)
% Vq: (Nk1 x Nk2 x Nk3 x m x m)
% rho_k: same shape
% OUTPUT: SigmaF_k (Nk1 x Nk2 x Nk3 x m x m)

[Nk1,Nk2,Nk3,m,~] = size(Vq);
SigmaF_k = zeros(Nk1,Nk2,Nk3,m,m);

% FFT convention: MATLAB fftn/ifftn with implicit 1/N factors (ifftn has 1/N)
% 我们用 ifftn(Vq).*ifftn(rho) 再 fftn 回去

for a = 1:m
  for b = 1:m
    Vq_ab   = squeeze(Vq(:,:,:,a,b));
    rho_ab  = squeeze(rho_k(:,:,:,a,b));
    Vr_ab   = ifftn(Vq_ab);
    rhor_ab = ifftn(rho_ab);
    Sig_ab  = - fftn( Vr_ab .* rhor_ab );
    SigmaF_k(:,:,:,a,b) = Sig_ab;
  end
end

% 可选：去掉 q=0 自作用（已在 build_Vq 里把 q=0=0 了，这里无需再处理）
% 若你在 build_Vq 没有处理 q=0，这里可以手动减掉平均值等。
end


function SigmaH_k = hartree_uniform_smallq(g, tau, rho_k, opts)
% Returns diagonal Hartree self-energy (k-independent), replicated on k-grid.
% SigmaH_{aa} = sum_gamma V_{aγ}(q->0) * ( nbar_γ - nref_γ )
% remove common-mode average.

[Nk1,Nk2,Nk3,m,~] = size(rho_k);
NkTot = Nk1*Nk2*Nk3;

% average occupation per orbital (per cell, total incl spin)
nbar = zeros(m,1);
for a=1:m
    nbar(a) = (1/NkTot) * sum( reshape(real(rho_k(:,:,:,a,a)), [],1) );
end

% reference (neutralizing background)
if isfield(opts,'nref') && numel(opts.nref)==m
    nref = opts.nref(:);
else
    nref = ones(m,1) * (sum(nbar)/m);
end
dn = nbar - nref;

% approximate V_ab(q->0) using a few smallest non-zero q
qsel = [ 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1 ];
Vsmall = zeros(m,m, size(qsel,1));

% layer direction: take z-component
z = tau(:,3);
for iq=1:size(qsel,1)
    qfrac = qsel(iq,:)./ [Nk1 Nk2 Nk3];
    qvec  = (g.b * qfrac.').';     % 1/Ang
    qmag  = norm(qvec);
    if qmag < 1e-14, continue; end
    switch lower(opts.vmodel)
        case '2d'
            Vscalar = 2*pi*opts.e2/(opts.eps*qmag);
        case 'keldysh'
            Vscalar = opts.e2 /( qmag*(1+opts.r0*qmag) );
        otherwise
            error('Unknown vmodel');
    end
    for a=1:m
      for b=1:m
        Vsmall(a,b,iq) = Vscalar * exp(-qmag*abs(z(a)-z(b)));
      end
    end
end
Vab0 = mean(Vsmall,3,'omitnan');   % m x m

% Hartree diagonal (k-independent)
SigH_diag = Vab0 * dn;             % m x 1

% remove common-mode (gauge)
SigH_diag = SigH_diag - mean(SigH_diag);

SigmaH_k = zeros(Nk1,Nk2,Nk3,m,m);
for a=1:m
    SigmaH_k(:,:,:,a,a) = SigH_diag(a);
end
end


function Etot = hf_total_energy(H0_k, SigmaH_k, SigmaF_k, rho_k)
% E_HF = (1/Nk) sum_k Tr[ H0(k) ρ(k) + 1/2 (ΣH(k)+ΣF(k)) ρ(k) ]

[Nk1,Nk2,Nk3,m,~] = size(rho_k);
NkTot = Nk1*Nk2*Nk3;
Etot = 0.0;

for i1=1:Nk1
  for i2=1:Nk2
    for i3=1:Nk3
      H0  = squeeze(H0_k(i1,i2,i3,:,:));
      SH  = squeeze(SigmaH_k(i1,i2,i3,:,:));
      SF  = squeeze(SigmaF_k(i1,i2,i3,:,:));
      rho = squeeze(rho_k(i1,i2,i3,:,:));
      % 保证数值实性
      H0  = (H0+H0')/2; SH = (SH+SH')/2; SF=(SF+SF')/2; rho=(rho+rho')/2;
      Etot = Etot + real(trace( H0*rho + 0.5*(SH+SF)*rho ));
    end
  end
end

Etot = Etot / NkTot;
end
