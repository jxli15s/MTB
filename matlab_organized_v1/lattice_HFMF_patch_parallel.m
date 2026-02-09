%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Construct the g.ham                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("Rgra_5s");
g = MTB.read_poscar(g,"data/Graphene/5s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/5s/wannier90_hr_p1.dat','data/Graphene/5s/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.wpos=g.atoms*g.a;
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
%%
s3=[1  0
    0  -1];
laf=-diag([0.002,0.002,0.001,0.001,0.0,0.0,-0.001,-0.001,-0.002,-0.002]);
Zeeman=kron(s3,laf);
g.add_zeeman(Zeeman)
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.9,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
hkpoints={[1/3,2/3,0.0]*0.9,...
          [1/3,2/3,0.0],...
          [1/3,2/3,0.0]+([0.5,0.5,0.0]-[1/3,2/3,0.0])*0.15,...
          };% hkpoints-high symmetry k points
nk=251;
efermi=-0.0;

Electric_field_in_evpA=0.00; %0.08-0.12 V/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
%%

kpoint=[1/3,2/3,0.0];
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
%%
figure()
plot(Energy)
%%
abs(Psik(:,10).*Psik(:,10))
%%
%%
% === User-prepared g ===
% g.a, g.b, g.hopr, g.ham, g.wpos ...
% 单位：e2/(4π ε0) = 14.399645 eV·Å
opts.e2     = 14.399645;     
opts.eps    = 27;             % 相对介电常数
opts.alpha  = 2*pi*opts.e2/opts.eps;   % 统一：V(q)=alpha/q
opts.r0     = 10;
opts.vmodel = '2d';          % 2d 或 'keldysh' 并设置 opts.r0
opts.qmin   = 1e-6;
opts.layer_z= g.wpos(:,3);
opts.nspin = 1;
opts.Nelec = 10; 
opts.parfor = true;
opts.mix    = 0.5;
opts.maxiter= 200;
opts.tol    = 1e-9;
opts.fock_drop_q0=true;

%%
% patch 网格
opts.Ncoarse = 5; opts.Ndense = 30; opts.R_dense = 0.08;
out=main_hfmf_patch(g, opts);
%%
figure()
for i=9:12
    hold on;
    plot(out.E_k(:,i))
end
%%
q_cart=[0.3,0.3,0.0]*g.b;
alpha=opts.alpha
layer_z=g.wpos(:,3)
V_ab=Vq_2d_layered_core(q_cart, alpha, layer_z, opts)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Function                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = main_hfmf_patch(g, opts)
% MAIN_HFMF_PATCH  Self-consistent HF (Wannier TB, non-uniform k-patch)
% Row-basis convention: rows of g.b/g.a are b_i/a_i (Cartesian).
%
% Inputs:
%   g.a        : [d x d], rows are real-space primitive vectors (a1; a2; [a3])
%   g.b        : [d x d], rows are reciprocal primitives (b1; b2; [b3])
%   g.hopr     : [L x d], hopping ΔR in fractional coords (lattice gauge)
%   g.ham      : [m x m x L], hopping matrices t(ΔR)
%   g.wpos     : [m x 3], orbital Cartesian positions (only z used by default)
%
%   opts.Ncoarse, opts.Ndense, opts.R_dense : patch mesh parameters (hex BZ)
%   opts.kT       : temperature (energy units)
%   opts.nspin    : spin degeneracy
%   opts.Nelec    : target electrons per cell (incl. spin)
%   opts.mix      : density mixing (0~0.8)
%   opts.maxiter  : SCF max iterations
%   opts.tol      : convergence tolerance on rho_orb
%   opts.parfor   : logical, parallelize over k in eig and Fock
%   opts.fock_drop_q0 : logical, drop q=0 in Fock (self-exchange regularization)
%
%   --- Unified V(q) model (set by make_Vq_model) ---
%   opts.vmodel   : '2d' | 'keldysh' | 'lattice'
%   opts.e2, opts.eps, opts.r0, opts.alpha, opts.layer_z, ...
%   (see make_Vq_model.m for details)
%
% Outputs (out):
%   .k_cart, .k_frac, .w_k, .Wnorm
%   .H0_k, .Hmf_k, .SigmaF_k, .SigmaH_k     : [m x m x Nk]
%   .E_k                                     : [Nk x m]
%   .U3d                                     : [m x m x Nk]
%   .mu
%   .rho_k                                   : [m x m x Nk] (per-cell)
%   .rho_orb                                 : [m x m]      (k-avg per-cell)
%   .history  : [iter, mu, err]
%
% Notes:
%   * This routine assumes g.b/g.a are row-based. If your g.a is column-based,
%     change R_cart = R_frac * g.a into R_cart = (g.a * R_frac.').' in build_H0k_from_wannier.

    % ---------- defaults ----------
    if ~isfield(opts,'kT'),       opts.kT = 0; end
    if ~isfield(opts,'nspin'),    opts.nspin = 2; end
    if ~isfield(opts,'mix'),      opts.mix = 0.5; end
    if ~isfield(opts,'maxiter'),  opts.maxiter = 50; end
    if ~isfield(opts,'tol'),      opts.tol = 1e-8; end
    if ~isfield(opts,'parfor'),   opts.parfor = false; end
    if ~isfield(opts,'fock_drop_q0'), opts.fock_drop_q0 = true; end

    % ---------- k-mesh (hex, your function expects column b1,b2) ----------
    b1_col = g.b(1,1:2).';                 % convert row to column (2×1)
    b2_col = g.b(2,1:2).';
    [k_cart, k_frac, w_k, ~] = make_kmesh_hex_patch(b1_col, b2_col, ...
        opts.Ncoarse, opts.Ndense, opts.R_dense);
    Nk    = size(k_cart,1);
    k_frac = [k_frac, zeros(Nk,1)]; % 2D->3D compatible
    k_cart = k_frac*g.b;
    Wnorm = sum(w_k);                     % = Ω_BZ


    % ---------- Build H0(k) from Wannier TB (row-basis aware) ----------
    [H0_k, E0_k, U0_3d] = build_H0k_from_wannier(g, k_frac, opts.parfor);

    % ---------- Unified V(q) model for BOTH Fock & Hartree ----------
    opts = make_Vq_model(opts, g);        % sets opts.Vq_fun / Vq_vec_fun consistently

    % ---------- Initial density from non-interacting bands ----------
    m   = size(g.ham,1);
    kT  = opts.kT;
    nsp = opts.nspin;

    mu  = find_mu_patch(E0_k, w_k, Wnorm, opts.Nelec, kT, nsp);
    f0  = fermi_dirac(E0_k, mu, kT);      % [Nk x m]
    rho_k = rho_from_UF(U0_3d, f0);       % [m x m x Nk], per-cell
    rho_orb = kavg_rho(rho_k, w_k, Wnorm);

    history = [];

    % ---------- SCF loop ----------
    for it = 1:opts.maxiter
        % Fock Σ_F(k): direct sum on patch mesh with q folding (row-basis)
        SigmaF_k = fock_direct_patch(g, k_frac, w_k, Wnorm, rho_k, opts);

        % Hartree Σ_H(k): k-independent, diagonal; uses the SAME V(q) model at q→0
        SigmaH_k = hartree_uniform_smallq_patch(g, g.wpos, rho_k, w_k, Wnorm, opts);

        % Total mean-field Hamiltonian
        Hmf_k = H0_k + SigmaF_k + SigmaH_k;

        % Diagonalize (optionally parfor)
        [E_k, U3d] = diag_all_k(Hmf_k, opts.parfor);

        % New chemical potential and occupations
        mu = find_mu_patch(E_k, w_k, Wnorm, opts.Nelec, kT, nsp);
        f  = fermi_dirac(E_k, mu, kT);

        % New density
        rho_k_new  = rho_from_UF(U3d, f);
        rho_orb_new= kavg_rho(rho_k_new, w_k, Wnorm);

        % Density mixing (linear)
        rho_k   = (1-opts.mix)*rho_k   + opts.mix*rho_k_new;
        rho_orb = (1-opts.mix)*rho_orb + opts.mix*rho_orb_new;

        % Convergence monitor
        err = norm(rho_orb_new - rho_orb, 'fro') / max(1, norm(rho_orb,'fro'));
        history = [history; it, mu, err]; %#ok<AGROW>
        fprintf('SCF %3d: mu = %.6f, err = %.3e\n', it, mu, err);

        if err < opts.tol, break; end
    end

    % ---------- pack outputs ----------
    out.k_cart   = k_cart;
    out.k_frac   = k_frac;
    out.w_k      = w_k;
    out.Wnorm    = Wnorm;

    out.H0_k     = H0_k;
    out.Hmf_k    = Hmf_k;
    out.SigmaF_k = SigmaF_k;
    out.SigmaH_k = SigmaH_k;

    out.E_k      = E_k;
    out.U3d      = U3d;
    out.mu       = mu;

    out.rho_k    = rho_k;
    out.rho_orb  = rho_orb;
    out.history  = history;
end

function rho_orb = kavg_rho(rho_k, w_k, Wnorm)
% K-AVERAGE of rho(k) to per-cell density matrix.
    [m,~,Nk] = size(rho_k);
    acc = zeros(m,m);
    for ik=1:Nk
        acc = acc + w_k(ik) * rho_k(:,:,ik);
    end
    rho_orb = (1/Wnorm) * acc;
end

function rho_k = rho_from_UF(U3d, f)
% Build rho(k) = U(k) diag(f_k) U(k)^\dagger for all k (3D arrays).
% Inputs:
%   U3d : [m x m x Nk], columns are eigenvectors at each k
%   f   : [Nk x m], occupations at each k
% Output:
%   rho_k : [m x m x Nk]
    [m,~,Nk] = size(U3d);
    F3 = zeros(m,m,Nk);
    for ik=1:Nk
        F3(:,:,ik) = diag(f(ik,:));
    end
    % f=diag([0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0]);
    %rho = U * F * U^\dagger   (page-wise)
    rho_k = pagemtimes( pagemtimes(U3d, F3), permute(conj(U3d), [2 1 3]) );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------Fock term------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SigmaF_k = fock_direct_patch(g, klist_frac, w_k, Wnorm, rho_k, opts)
% Σ_F_ab(k_i) = -(1/Ω_BZ) ∑_{k'} w(k') V_ab(q) ρ_ab(k'), q = fold(k_i-k')
% 行基：q_cart = q_frac * g.b

    Brows = g.b;                               % 行基
    Nk  = size(klist_frac,1);
    m   = size(rho_k,1);

    if ~isfield(opts,'fock_drop_q0'), opts.fock_drop_q0 = false; end
    if ~isfield(opts,'parfor'),       opts.parfor       = false; end

    use_vecV = isfield(opts,'Vq_vec_fun') && ~isempty(opts.Vq_vec_fun);
    SigmaF_k = zeros(m,m,Nk, 'like', rho_k);

    if opts.parfor
        % 直接在 parfor 里调用局部函数，不用函数句柄变量
        parfor i = 1:Nk
            SigmaF_k(:,:,i) = local_one(i, klist_frac, Brows, w_k, Wnorm, rho_k, opts, use_vecV);
        end
    else
        for i = 1:Nk
            SigmaF_k(:,:,i) = local_one(i, klist_frac, Brows, w_k, Wnorm, rho_k, opts, use_vecV);
        end
    end
end

function Sigma_i = local_one(i, kfrac, Brows, w_k, Wnorm, rho_k, opts, use_vecV)
    Nk = size(kfrac,1);
    m  = size(rho_k,1);
    Sigma_i = zeros(m,m, 'like', rho_k);

    ki_frac = kfrac(i,:);                  % 1×d
    q_frac  = ki_frac - kfrac;             % Nk×d
    % q_frac=kfrac-ki_frac;
    q_frac  = q_frac - round(q_frac);      % 折回 (-0.5,0.5]
    q_cart  = q_frac * Brows;              % Nk×d   行基

    if use_vecV
        Vpages = opts.Vq_vec_fun(q_cart, opts);    % m×m×Nk
        if opts.fock_drop_q0
            iz = vecnorm(q_cart,2,2) < 1e-14;
            if any(iz), Vpages(:,:,iz) = 0; end
        end
        for j = 1:Nk
            Sigma_i = Sigma_i - (w_k(j)/Wnorm) * ( Vpages(:,:,j) .* rho_k(:,:,j) );
        end
    else
        for j = 1:Nk
            if opts.fock_drop_q0 && all(abs(q_cart(j,:))<1e-14), continue; end
            V_ab = opts.Vq_fun(q_cart(j,:), opts); % m×m
            Sigma_i = Sigma_i - (w_k(j)/Wnorm) * ( V_ab .* rho_k(:,:,j) );
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------Hartree term--------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SigmaH_k = hartree_uniform_smallq_patch(g, tau, rho_k, w_k, Wnorm, opts)
% HARTREE_UNIFORM_SMALLQ_PATCH
% k-independent, diagonal Hartree shift with layer resolution; remove common mode.
% Uses the SAME V(q) model as Fock by sampling q→0 along each reciprocal basis direction.
%
% Inputs:
%   g.b    : [d x d] row-basis
%   tau    : [m x 3] orbital Cartesian positions (only z used here)
%   rho_k  : [m x m x Nk]  (per-cell)
%   w_k    : [Nk x 1], sum = Ω_BZ
%   Wnorm  : scalar, Ω_BZ
%   opts.Vq_fun : @(q_cart, opts) m×m matrix, unified with Fock
%   opts.nref   : optional [m x 1] reference density (per-cell)
%
% Output:
%   SigmaH_k : [m x m x Nk] (each page identical; only diagonal entries nonzero)

    [m,~,Nk] = size(rho_k);
    w_k = w_k(:);

    % k-averaged orbital occupations per cell: nbar(a) = (1/Ω) ∑_k w(k) ρ_aa(k)
    rho_flat = reshape(real(rho_k), m*m, Nk); % (m*m)×Nk
    idx_diag = 1:(m+1):m*m;
    rho_diag = rho_flat(idx_diag, :).';       % Nk×m
    nbar = (1/Wnorm) * (w_k.' * rho_diag).';  % m×1

    % reference background (e.g., neutrality per layer)
    if isfield(opts,'nref') && numel(opts.nref)==m
        nref = opts.nref(:);
    else
        nref = mean(nbar) * ones(m,1);
    end
    dn = nbar-nref;                        % m×1

    % V_ab(q→0): sample tiny q along each reciprocal direction (row-basis)
    d = size(g.b,1)-1;
    Vsmall = zeros(m,m,d);
    eps_q = 1e-4;
    for iq=1:d
        qf = zeros(1,3); qf(iq) = eps_q;      % tiny fractional step
        q_cart = qf * g.b;                    % 1×d
        Vsmall(:,:,iq) = opts.Vq_fun(q_cart, opts);  % SAME model as Fock
    end
    Vab0 = mean(Vsmall(:,:,1:d), 3, 'omitnan');      % m×m
     % Vab0 = opts.Vq_fun([0,0,0], opts); % m×m

    % Diagonal Hartree shift; remove common mode to avoid trivial rigid shift
    SigH_diag = Vab0 * dn;                    % m×1
    SigH_diag = SigH_diag - mean(SigH_diag);

    D = diag(SigH_diag);                      % m×m
    SigmaH_k = repmat(D, [1 1 Nk]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------Get_H0-----------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H0_k, E_k, U3d] = build_H0k_from_wannier(g, k_frac, use_parfor)
% BUILD_H0K_FROM_WANNIER
% Lattice-gauge TB: H0(k) = sum_{ΔR} t(ΔR) e^{ i k·ΔR }
% Row-basis convention: k_cart = k_frac * g.b,  ΔR_cart = ΔR_frac * g.a
%
% Inputs:
%   g.a      : [d x d], rows are a_i
%   g.b      : [d x d], rows are b_i
%   g.hopr   : [L x d], ΔR in fractional coords
%   g.ham    : [m x m x L], hopping t(ΔR)
%   k_frac   : [Nk x d], fractional k
%   use_parfor : logical
%
% Outputs:
%   H0_k     : [m x m x Nk]
%   E_k      : [Nk x m]
%   U3d      : [m x m x Nk]

    L  = size(g.hopr,1);
    m  = size(g.ham,1);
    Nk = size(k_frac,1);
    H0_k = zeros(m,m,Nk);

    for ik = 1:Nk
        kf = k_frac(ik,:);                % 1×d
        k_cart = kf * g.b;                % 1×d   (row-basis)

        Hk = zeros(m,m);
        for ell = 1:L
            dR_frac = g.hopr(ell,:);      % 1×d
            % If g.a is column-based in your data, replace next line by:
            dR_cart = dR_frac * g.a;      % 1×d   (row-basis for a)
            phase = exp(1i * dot(k_cart, dR_cart));
            Hk = Hk + g.ham(:,:,ell) * phase;
        end
        H0_k(:,:,ik) = Hk;
    end

    [E_k, U3d] = diag_all_k(H0_k, use_parfor);
end

function [E_k, U3d] = diag_all_k(Hk, use_parfor)
% DIAG_ALL_K  Diagonalize H(k) for all k, optionally with PARFOR.
% Inputs:
%   Hk         : [m x m x Nk]
%   use_parfor : logical (default: false)
% Outputs:
%   E_k        : [Nk x m]  energies sorted ascending
%   U3d        : [m x m x Nk]  eigenvectors per k (columns are eigenstates)

    if nargin < 2, use_parfor = false; end

    [m,~,Nk] = size(Hk);
    E_k = zeros(Nk, m);
    U3d = zeros(m,m,Nk);

    if use_parfor
        parfor ik = 1:Nk
            H = (Hk(:,:,ik) + Hk(:,:,ik)')/2;           % Hermitize slightly
            [U, d] = eig(H, 'vector');                  % d: eigenvalues
            [e, idx] = sort(real(d), 'ascend');
            U = U(:, idx);
            U3d(:,:,ik) = U;
            E_k(ik,:)   = e.';
        end
    else
        for ik = 1:Nk
            H = (Hk(:,:,ik) + Hk(:,:,ik)')/2;
            [U, d] = eig(H, 'vector');
            [e, idx] = sort(real(d), 'ascend');
            U = U(:, idx);
            U3d(:,:,ik) = U;
            E_k(ik,:)   = e.';
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------Chemical Potential & Fermi level mu----------------- -------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fermi_dirac(E, mu, kT)
% FERMI_DIRAC  Occupations per k and band.
% Inputs:
%   E  : [Nk x m]
%   mu : scalar
%   kT : scalar
% Output:
%   f  : [Nk x m]
    if kT <= 0
        f = double(E <= mu);
    else
        x = (E - mu)/kT;
        f = 1 ./ (1 + exp(x));
    end
end

function mu = find_mu_patch(E_k, w_k, Wnorm, Nelec, kT, nspin)
% FIND_MU_PATCH  Bisection on μ s.t. nspin*(1/Ω)∑_k w(k)∑_n f_{kn} = Nelec.
    emin = min(E_k(:)) - 5*max(kT,1e-4);
    emax = max(E_k(:)) + 5*max(kT,1e-4);
    for it=1:80
        mu = 0.5*(emin+emax);
        f  = fermi_dirac(E_k, mu, kT);                 % Nk x m
        Ne = nspin * (1/Wnorm) * sum( w_k(:) .* sum(f,2) );
        if Ne > Nelec, emax = mu; else, emin = mu; end
        if abs(Ne-Nelec) < 1e-12*max(1,Nelec), break; end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------2D Coulomb & Layer screening-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opts = make_Vq_model(opts, g)
% MAKE_VQ_MODEL  Build unified V(q) for both Fock and Hartree.
% Sets opts.Vq_fun (and Vq_vec_fun) so both terms use identical interaction.
%
% Inputs (typical):
%   opts.vmodel  : '2d' | 'keldysh' | 'lattice'
%   opts.e2      : e^2/(4π ε0) in your units (e.g., 14.3996 eV·Å)
%   opts.eps     : ε_r (relative dielectric, for '2d' / 'keldysh')
%   opts.r0      : Keldysh screening length (for 'keldysh'), length units
%   opts.alpha   : if empty, set alpha = 2π*e2/eps
%   opts.layer_z : [m x 1] layer z per orbital (default: g.wpos(:,3))
%
%   For 'lattice' (discrete Fourier):
%     opts.VR_R_cart : [NR x d] real-space vectors (Cartesian)
%     opts.VR_ab     : cell(NR,1), each [m x m] matrix V_ab(R)
%
% Outputs:
%   opts.Vq_fun     : @(q_cart, opts) -> [m x m]
%   opts.Vq_vec_fun : @(q_cart_all, opts) -> [m x m x Nk]  (optional, for speed)

    if ~isfield(opts,'layer_z') || isempty(opts.layer_z)
        opts.layer_z = g.wpos(:,3);
    end
    if ~isfield(opts,'vmodel'), opts.vmodel = '2d'; end

    switch lower(opts.vmodel)
        case '2d'
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;         % V(q)=alpha/q
            end
            opts.Vq_fun     = @(q_cart, opts2) Vq_2d_layered_core(q_cart, opts2.alpha, opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all, opts2) Vq_2d_layered_vec_core(q_cart_all, opts2.alpha, opts2.layer_z, opts2);

        case 'keldysh'
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;         % prefactor unified
            end
            opts.Vq_fun     = @(q_cart, opts2) Vq_keldysh_layered_core(q_cart, opts2.alpha, opts2.r0, opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all, opts2) Vq_keldysh_layered_vec_core(q_cart_all, opts2.alpha, opts2.r0, opts2.layer_z, opts2);

        case 'lattice'
            assert(isfield(opts,'VR_R_cart') && isfield(opts,'VR_ab'), ...
                'For vmodel=lattice, provide opts.VR_R_cart and opts.VR_ab.');
            opts.Vq_fun     = @(q_cart, opts2) Vq_from_VR_core(q_cart, opts2.VR_R_cart, opts2.VR_ab);
            opts.Vq_vec_fun = @(q_cart_all, opts2) Vq_from_VR_vec_core(q_cart_all, opts2.VR_R_cart, opts2.VR_ab);

        otherwise
            error('Unknown vmodel: %s', opts.vmodel);
    end
end

% ---- cores (shared) ----
function V_ab = Vq_2d_layered_core(q_cart, alpha, layer_z, opts)
    q = max(norm(q_cart), get_qmin(opts));
    DZ = abs(layer_z - layer_z.');
    V_ab = (alpha / q) * exp(-q * DZ);
end

function V_pages = Vq_2d_layered_vec_core(q_cart_all, alpha, layer_z, opts)
    Nk   = size(q_cart_all,1);
    qmag = vecnorm(q_cart_all,2,2);
    qmag = max(qmag, get_qmin(opts));
    DZ   = abs(layer_z - layer_z.');
    m    = numel(layer_z);
    V_pages = zeros(m,m,Nk);
    for j=1:Nk
        V_pages(:,:,j) = (alpha / qmag(j)) * exp(-qmag(j) * DZ);
    end
end

function V_ab = Vq_keldysh_layered_core(q_cart, alpha, r0, layer_z, opts)
    q = max(norm(q_cart), get_qmin(opts));
    DZ = abs(layer_z - layer_z.');
    V_ab = (alpha / (q*(1 + r0*q))) * exp(-q * DZ);
end

function V_pages = Vq_keldysh_layered_vec_core(q_cart_all, alpha, r0, layer_z, opts)
    Nk   = size(q_cart_all,1);
    qmag = vecnorm(q_cart_all,2,2);
    qmag = max(qmag, get_qmin(opts));
    DZ   = abs(layer_z - layer_z.');
    m    = numel(layer_z);
    V_pages = zeros(m,m,Nk);
    for j=1:Nk
        V_pages(:,:,j) = (alpha / (qmag(j)*(1 + r0*qmag(j)))) * exp(-qmag(j) * DZ);
    end
end

function V_ab = Vq_from_VR_core(q_cart, R_cart, V_ab_of_R)
% Fourier sum: V_ab(q) = ∑_R V_ab(R) e^{-i q·R}
    NR = size(R_cart,1);
    V_ab = 0;
    for ir=1:NR
        V_ab = V_ab + exp(-1i * dot(q_cart, R_cart(ir,:))) * V_ab_of_R{ir};
    end
    V_ab = real(V_ab);  % enforce Hermiticity (numerical)
end

function V_pages = Vq_from_VR_vec_core(q_cart_all, R_cart, V_ab_of_R)
    Nk  = size(q_cart_all,1);
    m   = size(V_ab_of_R{1},1);
    V_pages = zeros(m,m,Nk);
    for j=1:Nk
        V_pages(:,:,j) = Vq_from_VR_core(q_cart_all(j,:), R_cart, V_ab_of_R);
    end
end

function qmin = get_qmin(opts)
    if isfield(opts,'qmin') && ~isempty(opts.qmin), qmin = opts.qmin; else, qmin = 1e-6; end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1%%%%%%%%%%
%-----------------Cretate pached BZ and w_k at (0,1)----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense)
%MAKE_KMESH_HEX_PATCH_01  k-mesh on hexagonal BZ (frac coords in [0,1))
%   Coarse mesh on full BZ + dense patches around K1=(1/3,2/3), K2=(2/3,1/3).
%
% INPUT:
%   b1, b2   : [2x1] reciprocal lattice vectors, k = k1*b1 + k2*b2
%   Ncoarse  : coarse grid density per direction (fractional coords in [0,1))
%   Ndense   : dense patch density per direction
%   R_dense  : half-width of dense patch (in fractional coords)
%
% OUTPUT:
%   k_cart   : [Nk x 2] cartesian k points (kx, ky)
%   k_frac   : [Nk x 2] fractional coords (k1, k2) in [0,1)
%   weight   : [Nk x 1] integration weights, sum(weight) = area_BZ
%   region_id: [Nk x 1], 1 = coarse region, 2 = K1 patch, 3 = K2 patch

    % ---- 0. 输入检查 ----
    if nargin < 5
        error('Usage: make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense)');
    end

    b1 = b1(:);  b2 = b2(:);
    if numel(b1) ~= 2 || numel(b2) ~= 2
        error('b1 and b2 must be 2x1 column vectors.');
    end

    B = [b1, b2];                 % 2x2
    area_BZ = abs(det(B));        % BZ 平行四边形面积

    % ---- 1. K1 / K2 的 frac 坐标（[0,1)）----
    % 按你说的约定：Gamma=(0,0), K1=(1/3,2/3), K2=(2/3,1/3)
    K1_frac = [1/3, 2/3];
    K2_frac = [2/3, 1/3];

    % ---- 2. whole-BZ coarse grid: frac in [0,1) x [0,1) ----
    d_coarse = 1 / Ncoarse;
    k1c = 0 : d_coarse : 1 - d_coarse;
    k2c = 0 : d_coarse : 1 - d_coarse;
    [K1c, K2c] = meshgrid(k1c, k2c);
    kfrac_coarse = [K1c(:), K2c(:)];   % Nc x 2, in [0,1)

    % ---- 3. 在 K1 / K2 周围构造 dense patch (fractional) ----
    d_dense = (2 * R_dense) / Ndense;
    d1 = -R_dense : d_dense : R_dense;
    d2 = -R_dense : d_dense : R_dense;
    [D1, D2] = meshgrid(d1, d2);

    kfrac_dense_K1 = [D1(:) + K1_frac(1), D2(:) + K1_frac(2)];
    kfrac_dense_K2 = [D1(:) + K2_frac(1), D2(:) + K2_frac(2)];

    % 折回 [0,1)（mod 1）
    kfrac_dense_K1 = kfrac_dense_K1 - floor(kfrac_dense_K1);
    kfrac_dense_K2 = kfrac_dense_K2 - floor(kfrac_dense_K2);

    % ---- 4. 从 coarse 网格中移除落在两个 dense patch 里的点 ----
    % patch 条件（直接用 frac 空间的矩形窗口）
    deltaK1 = abs(kfrac_coarse - K1_frac);
    inK1    = (deltaK1(:,1) <= R_dense) & (deltaK1(:,2) <= R_dense);

    deltaK2 = abs(kfrac_coarse - K2_frac);
    inK2    = (deltaK2(:,1) <= R_dense) & (deltaK2(:,2) <= R_dense);

    in_patch = inK1 | inK2;
    kfrac_coarse_out = kfrac_coarse(~in_patch, :);

    % ---- 5. 合并 & 去重 ----
    kfrac_all = [kfrac_coarse_out;
                 kfrac_dense_K1;
                 kfrac_dense_K2];

    key = round(kfrac_all * 1e8);
    [~, ia] = unique(key, 'rows', 'stable');
    k_frac = kfrac_all(ia, :);          % 所有 frac 点 (k1,k2) in [0,1)

    Nk = size(k_frac, 1);

    % ---- 6. region_id: 1 coarse, 2 K1, 3 K2 ----
    region_id = ones(Nk, 1);

    deltaK1_all = abs(k_frac - K1_frac);
    inK1_all    = (deltaK1_all(:,1) <= R_dense) & (deltaK1_all(:,2) <= R_dense);

    deltaK2_all = abs(k_frac - K2_frac);
    inK2_all    = (deltaK2_all(:,1) <= R_dense) & (deltaK2_all(:,2) <= R_dense);

    region_id(inK1_all) = 2;
    region_id(inK2_all) = 3;

    % ---- 7. 积分权重（先粗略按区域分配，再归一化）----
    area_patch_frac  = (2 * R_dense)^2;   % 每个 patch 的 frac 面积
    area_coarse_frac = 1 - 2 * area_patch_frac;

    idx_coarse = (region_id == 1);
    idx_K1     = (region_id == 2);
    idx_K2     = (region_id == 3);

    N_coarse = nnz(idx_coarse);
    N_K1     = nnz(idx_K1);
    N_K2     = nnz(idx_K2);

    weight_frac = zeros(Nk, 1);

    if N_coarse > 0 && area_coarse_frac > 0
        weight_frac(idx_coarse) = area_coarse_frac / N_coarse;
    else
        weight_frac(idx_coarse) = 1;   % 极端情况占位
    end
    if N_K1 > 0
        weight_frac(idx_K1) = area_patch_frac / N_K1;
    end
    if N_K2 > 0
        weight_frac(idx_K2) = area_patch_frac / N_K2;
    end

    % 全局归一化：sum(weight) = area_BZ
    S_frac = sum(weight_frac);
    if S_frac == 0
        weight_frac(:) = 1 / Nk;
        S_frac = 1;
    end
    weight_frac = weight_frac / S_frac;
    weight = area_BZ * weight_frac;

    % ---- 8. 变成真实 k 向量 ----
    k_cart = (B * k_frac.').';     % Nk x 2, (kx,ky)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1%%%%%%%%%%
%-----------------Cretate pached BZ and w_k at (-0.5,0.5)-----------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense)
%MAKE_KMESH_HEX_PATCH  Build a k-mesh on a hexagonal lattice BZ
%   with coarse sampling everywhere and dense patches around K and K'.
%
% INPUT:
%   b1, b2   : [2x1] real-space column vectors of reciprocal lattice
%              (units: 1/length). They define k = k1*b1 + k2*b2.
%   Ncoarse  : scalar, coarse grid points per direction in fractional
%              coordinates (k1,k2) in [-0.5,0.5).
%   Ndense   : scalar, dense patch points per direction around each valley.
%   R_dense  : half-width of the dense patch (in fractional coordinates).
%              Patch in frac-space: |k1 - K1| <= R_dense, |k2 - K2| <= R_dense.
%
% OUTPUT:
%   k_cart   : [Nk x 2] array, cartesian (kx, ky) of all k-points.
%   k_frac   : [Nk x 2] array, fractional (k1, k2) such that k = k1*b1 + k2*b2.
%   weight   : [Nk x 1] integration weights (in k-space area units).
%              Sum_j weight(j) = area of the BZ parallelogram = |det([b1 b2])|.
%   region_id: [Nk x 1] integer flag:
%              1 = coarse region
%              2 = dense patch around K
%              3 = dense patch around K'
%
% NOTES:
%   - The BZ here is taken as the parallelogram in fractional coords
%     (k1, k2) in [-0.5, 0.5) x [-0.5, 0.5).
%   - K, K' are placed at (k1,k2) = (2/3, 1/3) and (1/3, 2/3),
%     and then folded into (-0.5, 0.5] via frac = frac - round(frac).
%   - R_dense should be small enough that the two patches do not overlap.
%
% EXAMPLE:
%   % Suppose real-space primitive vectors are:
%   a = 1.0;
%   a1 = a * [1; 0];
%   a2 = a * [0.5; sqrt(3)/2];
%   % Reciprocal lattice:
%   A  = [a1, a2];
%   B  = 2*pi * inv(A).';   % columns b1, b2
%   b1 = B(:,1); b2 = B(:,2);
%
%   [k_cart, k_frac, w, region_id] = make_kmesh_hex_patch(b1, b2, ...
%                                     20, 40, 0.08);
%
%   % k_cart, w 就可以直接用在 HFMF 的 k 积分里:
%   % sum_j w(j) * f(k_cart(j,:))

    % -------- 0. 基本检查 --------
    if nargin < 5
        error('Usage: make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense)');
    end

    b1 = b1(:);
    b2 = b2(:);
    if numel(b1) ~= 2 || numel(b2) ~= 2
        error('b1 and b2 must be 2x1 column vectors.');
    end

    % Reciprocal lattice matrix and its determinant (BZ area)
    B = [b1, b2];                 % 2x2
    area_BZ = abs(det(B));        % area of parallelogram BZ in k-space

    % -------- 1. 定义 K / K' 的 fractional 坐标并折回 --------
    % 原始 frac 坐标
    K_frac  = [2/3, 1/3];
    Kp_frac = [1/3, 2/3];

    % 折回到 (-0.5, 0.5] 区间（mod 1）
    K_frac  = K_frac  - round(K_frac);
    Kp_frac = Kp_frac - round(Kp_frac);

    % -------- 2. 构造 whole-BZ coarse grid (fractional coords) --------
    d_coarse = 1 / Ncoarse;  % step in frac space
    k1c = -0.5 : d_coarse : 0.5 - d_coarse;  % [-0.5, 0.5)
    k2c = -0.5 : d_coarse : 0.5 - d_coarse;
    [K1c, K2c] = meshgrid(k1c, k2c);
    kfrac_coarse = [K1c(:), K2c(:)];   % Nc x 2

    % -------- 3. 在 K / K' 周围构造 dense patch (fractional coords) --------
    % patch 是中心在 K / K'，边长 2*R_dense 的方形区域
    d_dense = (2*R_dense) / Ndense;
    d1 = -R_dense : d_dense : R_dense;
    d2 = -R_dense : d_dense : R_dense;
    [D1, D2] = meshgrid(d1, d2);

    kfrac_dense_K  = [D1(:) + K_frac(1),  D2(:) + K_frac(2)];
    kfrac_dense_Kp = [D1(:) + Kp_frac(1), D2(:) + Kp_frac(2)];

    % 再次折回到 (-0.5, 0.5] 以保证都在同一个 fundamental domain
    kfrac_dense_K  = kfrac_dense_K  - round(kfrac_dense_K);
    kfrac_dense_Kp = kfrac_dense_Kp - round(kfrac_dense_Kp);

    % -------- 4. 从 coarse 网格中移除落在两个 dense patch 里的点 --------
    % 使用 frac 空间定义 patch: |k1 - K1| <= R_dense 且 |k2 - K2| <= R_dense
    % 注意这里已经在 (-0.5,0.5] 里，因此不考虑额外的周期性
    deltaK  = abs(kfrac_coarse - K_frac);   % Nc x 2
    inK     = (deltaK(:,1) <= R_dense) & (deltaK(:,2) <= R_dense);

    deltaKp = abs(kfrac_coarse - Kp_frac);
    inKp    = (deltaKp(:,1) <= R_dense) & (deltaKp(:,2) <= R_dense);

    in_patch = inK | inKp;

    kfrac_coarse_out = kfrac_coarse(~in_patch, :);   % coarse outside patches

    % -------- 5. 合并 coarse_out + dense_K + dense_Kp 并去重 --------
    kfrac_all = [kfrac_coarse_out;
                 kfrac_dense_K;
                 kfrac_dense_Kp];

    % 为了 unique 稳定一点，先做缩放+round 再 unique
    key = round(kfrac_all * 1e8);
    [~, idx_unique] = unique(key, 'rows', 'stable');
    kfrac = kfrac_all(idx_unique, :);
    % 
    % % -------- 6. 标记 region_id: coarse(1), K_patch(2), Kp_patch(3) --------
    % % 我们再次用 frac 空间的窗口来检查
    % Nk = size(kfrac, 1);
    % region_id = zeros(Nk, 1);
    % 
    % % coarse region: 先设为 1
    % region_id(:) = 1;
    % 
    % % K patch
    % deltaK_all  = abs(kfrac - K_frac);
    % inK_all     = (deltaK_all(:,1) <= R_dense) & (deltaK_all(:,2) <= R_dense);
    % 
    % % K' patch
    % deltaKp_all = abs(kfrac - Kp_frac);
    % inKp_all    = (deltaKp_all(:,1) <= R_dense) & (deltaKp_all(:,2) <= R_dense);
    % 
    % region_id(inK_all)  = 2;
    % region_id(inKp_all) = 3;

    % -------- 6. 标记 region_id: coarse(1), K_patch(2), Kp_patch(3) --------
Nk = size(kfrac, 1);
region_id = ones(Nk, 1);   % 先全部设为 coarse=1

% 给边界留一点浮点余量
tol = 10 * eps(R_dense);   % 或者 tol = 1e-10; 根据你 R_dense 大小调整

% K patch
deltaK_all  = abs(kfrac - K_frac);
inK_all     = (deltaK_all(:,1) <= R_dense + tol) & ...
              (deltaK_all(:,2) <= R_dense + tol);

% K' patch
deltaKp_all = abs(kfrac - Kp_frac);
inKp_all    = (deltaKp_all(:,1) <= R_dense + tol) & ...
              (deltaKp_all(:,2) <= R_dense + tol);

region_id(inK_all)  = 2;
region_id(inKp_all) = 3;


    % -------- 7. 计算积分权重 --------
    % frac-space 中 BZ 面积是 1.0 (平行四边形 [0,1)x[0,1) 或 [-0.5,0.5)x[-0.5,0.5))
    % 我们有两类区域：
    %   - coarse outside patches
    %   - two dense patches around K / K'
    %
    % coarse 区域的面积（frac-space）:
    %   area_total_frac = 1
    %   area_patch_K_frac  = (2*R_dense)^2
    %   area_patch_Kp_frac = (2*R_dense)^2
    %   area_coarse_frac   = 1 - area_patch_K_frac - area_patch_Kp_frac

    area_patch_K_frac  = (2*R_dense)^2;
    area_patch_Kp_frac = (2*R_dense)^2;
    area_coarse_frac   = 1 - area_patch_K_frac - area_patch_Kp_frac;

    if area_coarse_frac < 0
        warning('R_dense too large: patches overlap or exceed BZ area in frac space.');
    end

    % 每个区域的点数
    idx_coarse = find(region_id == 1);
    idx_K      = find(region_id == 2);
    idx_Kp     = find(region_id == 3);

    N_coarse = numel(idx_coarse);
    N_K      = numel(idx_K);
    N_Kp     = numel(idx_Kp);

    weight_frac = zeros(Nk, 1);

    % 均匀分配每个区域的 frac-area
    if N_coarse > 0
        w_coarse = area_coarse_frac / N_coarse;
        weight_frac(idx_coarse) = w_coarse;
    end
    if N_K > 0
        w_K = area_patch_K_frac / N_K;
        weight_frac(idx_K) = w_K;
    end
    if N_Kp > 0
        w_Kp = area_patch_Kp_frac / N_Kp;
        weight_frac(idx_Kp) = w_Kp;
    end

    % 从 frac-space 面积转换到真实 k-space 面积（乘以 |det(B)|）
    weight = area_BZ * weight_frac;

    % -------- 8. 变成真实 k 向量 --------
    k_cart = (B * kfrac.').';   % Nk x 2, each row (kx, ky)
    k_frac=kfrac;

end




