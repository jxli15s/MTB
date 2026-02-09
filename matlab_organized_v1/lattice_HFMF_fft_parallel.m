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
kpoint=[1/3,2/3,0.0];
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
%%
s3=[1  0
    0  -1];
laf=-diag([0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,-0.001,-0.001]);
Zeeman=kron(s3,laf);
g.add_zeeman(Zeeman)
%%
Nk1 = 100; Nk2 = 100;
opts.kT      = 0;
opts.nspin   = 1;
opts.mix     = 0.5;
opts.maxiter = 100;
opts.tol     = 1e-6;
opts.parfor  = true;
opts.Nelec = 10; % Default half-filling for electrons
opts.vmodel  = 'keldysh';       % or 'keldysh' / 'lattice'
opts.r0      =  80;
opts.e2      = 14.3996;
opts.eps     = 100.0;
opts.fock_drop_q0 = true;

% 可选：patch 半宽 (fractional)
opts.R_patch = 0.1;
opts.qmin    = 1/Nk1*opts.R_patch*norm(g.b(1,:))*0.2;

% 电子数（默认半填）:
m = size(g.ham,1); opts.Nelec = m*opts.nspin/2;

out = main_hfmf_fft_patchK(g, opts, Nk1, Nk2);
%%
%%
ek=reshape(out.E_k,Nk1,Nk2,[]);
figure()
hold on;
surf(ek(:,:,10))
surf(ek(:,:,11))
%%

%%
figure()
for i=9:12
    hold on;
    plot(out.E_k(:,i));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Hartree-Fock on a uniform K-valley patch using FFT Fock          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = main_hfmf_fft_patchK(g, opts, Nk1, Nk2)
% MAIN_HFMF_FFT_PATCHK
%   Self-consistent Hartree–Fock on a *uniform rectangular k-mesh*
%   in a small patch around the K valley, using FFT-based Fock term.
%
%   Patch is centered at K_frac = (1/3, 2/3) in fractional coords of g.b.
%   k-grid: Nk1 x Nk2 uniform mesh in that small rectangle.
%
%   Row-basis convention:
%       g.a, g.b are [d x d], each row is a_i / b_i (Cartesian).
%
% INPUT:
%   g.a      [d x d]  real-space primitive vectors (rows)
%   g.b      [d x d]  reciprocal primitive vectors (rows)
%   g.hopr   [L x d]  hopping ΔR in fractional coords
%   g.ham    [m x m x L]  hopping matrices t(ΔR)
%   g.wpos   [m x 3]  orbital positions (Cartesian, z 包含层信息)
%
%   opts.kT       : temperature (energy units)
%   opts.nspin    : spin degeneracy
%   opts.Nelec    : target electrons per cell (incl. spin) [可选, 默认半填]
%   opts.mix      : density mixing (0~0.8)
%   opts.maxiter  : SCF max iterations
%   opts.tol      : convergence tolerance on rho_orb
%   opts.parfor   : logical, parallelize diagonalization
%
%   Patch 参数:
%     opts.R_patch : patch 半宽 (fractional units along b1,b2)，
%                    patch 区域为 [1/3±R_patch]×[2/3±R_patch]，默认 0.08.
%
%   Interaction model (统一给 Fock + Hartree):
%     opts.vmodel  : '2d' | 'keldysh' | 'lattice'
%     opts.e2      : e^2/(4π ε0) in your units (e.g., 14.3996 eV·Å)
%     opts.eps     : ε_r (relative dielectric)
%     opts.r0      : Keldysh screening length (for 'keldysh'), length units
%     opts.layer_z : optional, [m x 1] layer z; 默认用 g.wpos(:,3)
%     opts.qmin    : small-q regulator (默认 1e-2 in |q|)
%
%   Fock 数值设置:
%     opts.fock_drop_q0 : logical, 是否在 Fock 中丢掉 q=0 自交换 (默认 true)
%     opts.qeps         : 判断"q≈0"的阈值 (默认 1e-8)
%
% OUTPUT (out):
%   .k_cart, .k_frac, .w_k, .Wnorm
%   .H0_k, .Hmf_k, .SigmaF_k, .SigmaH_k   [m x m x Nk]
%   .E_k                                  [Nk x m]
%   .U3d                                  [m x m x Nk]
%   .mu
%   .rho_k                                [m x m x Nk]
%   .rho_orb                              [m x m]
%   .history   [iter, mu, err]
%   .Nk1, .Nk2

    % ---------- defaults ----------
    if ~isfield(opts,'kT'),            opts.kT       = 0;        end
    if ~isfield(opts,'nspin'),         opts.nspin    = 2;        end
    if ~isfield(opts,'mix'),           opts.mix      = 0.5;      end
    if ~isfield(opts,'maxiter'),       opts.maxiter  = 50;       end
    if ~isfield(opts,'tol'),           opts.tol      = 1e-8;     end
    if ~isfield(opts,'parfor'),        opts.parfor   = false;    end
    if ~isfield(opts,'fock_drop_q0'),  opts.fock_drop_q0 = true; end
    if ~isfield(opts,'vmodel'),        opts.vmodel   = '2d';     end
    if ~isfield(opts,'qmin') || isempty(opts.qmin),  opts.qmin = 1e-2; end
    if ~isfield(opts,'R_patch') || isempty(opts.R_patch), opts.R_patch = 0.08; end

    R_frac = opts.R_patch;

    % ---------- k-mesh: uniform patch around K = (1/3,2/3) ----------
    [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = ...
        make_kmesh_valley_patch_fft(g, Nk1, Nk2, R_frac);
    Nk = size(k_cart,1);

    % ---------- H0(k) from Wannier TB ----------
    [H0_k, E0_k, U0_3d] = build_H0k_from_wannier(g, k_frac, opts.parfor);
    m   = size(g.ham,1);
    kT  = opts.kT;
    nsp = opts.nspin;

    % ---------- electrons per cell ----------
    if ~isfield(opts,'Nelec') || isempty(opts.Nelec)
        opts.Nelec = m * nsp / 2;    % 默认半填 (单 valley)
    end

    % ---------- Unified V(q) model for BOTH Fock & Hartree ----------
    opts = make_Vq_model(opts, g);   % 使用单-valley g 设置 V(q)

    % ---------- initial density from non-interacting bands ----------
    mu  = find_mu_patch(E0_k, w_k, Wnorm, opts.Nelec, kT, nsp);
    f0  = fermi_dirac(E0_k, mu, kT);      % [Nk x m]
    rho_k   = rho_from_UF(U0_3d, f0);     % [m x m x Nk]
    rho_orb = kavg_rho(rho_k, w_k, Wnorm);

    % ---------- Precompute V(q) -> Vr_Fock (FFT kernel) ----------
    % 注意：这里对于 Fock 的 q-grid，我们把 g.b 的前两行缩放 2*R_frac，
    % 等价于把 patch 视为一个"有效 BZ"，q 步长 = (2R_frac/Nk1, 2R_frac/Nk2)。
    B_patch        = g.b;
    B_patch(1,:)   = 2*R_frac * B_patch(1,:);
    B_patch(2,:)   = 2*R_frac * B_patch(2,:);
    Vr_Fock = precompute_Vr_Fock_uniform_fromB(B_patch, m, Nk1, Nk2, opts);

    history = [];
    SigmaF_k = zeros(m,m,Nk);
    SigmaH_k = zeros(m,m,Nk);
    Hmf_k    = H0_k;
    E_k      = E0_k;
    U3d      = U0_3d;

    % ---------- SCF loop ----------
    for it = 1:opts.maxiter
        % 1) Fock Σ_F(k): FFT-based convolution on uniform K-patch
        SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2);

        % 2) Hartree Σ_H(k): k-independent diagonal term,
        %    使用统一 V(q) 模型在 q→0 极限
        SigmaH_k = hartree_uniform_smallq_patch(g, g.wpos, rho_k, w_k, Wnorm, opts);

        % 3) Total mean-field Hamiltonian
        Hmf_k = H0_k + SigmaF_k + SigmaH_k;

        % 4) Diagonalize Hmf(k)
        [E_k, U3d] = diag_all_k(Hmf_k, opts.parfor);

        % 5) Update μ and occupations
        mu = find_mu_patch(E_k, w_k, Wnorm, opts.Nelec, kT, nsp);
        f  = fermi_dirac(E_k, mu, kT);

        % 6) New density
        rho_k_new   = rho_from_UF(U3d, f);
        rho_orb_new = kavg_rho(rho_k_new, w_k, Wnorm);

        % 7) Mixing
        rho_k   = (1-opts.mix)*rho_k   + opts.mix*rho_k_new;
        rho_orb = (1-opts.mix)*rho_orb + opts.mix*rho_orb_new;

        % 8) Convergence monitor
        err = norm(rho_orb_new - rho_orb, 'fro') / max(1, norm(rho_orb,'fro'));
        history = [history; it, mu, err]; %#ok<AGROW>
        fprintf('SCF(FFT-K) %3d : mu = %.6f, err = %.3e\n', it, mu, err);

        if err < opts.tol
            break;
        end
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
    out.Nk1      = Nk1;
    out.Nk2      = Nk2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- k-average & density builder ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho_orb = kavg_rho(rho_k, w_k, Wnorm)
% K-AVERAGE of rho(k) to per-cell density matrix:
%   ρ_orb = (1/Ω) ∑_k w(k) ρ(k)
    [m,~,Nk] = size(rho_k);
    acc = zeros(m,m);
    for ik=1:Nk
        acc = acc + w_k(ik) * rho_k(:,:,ik);
    end
    rho_orb = (1/Wnorm) * acc;
end

function rho_k = rho_from_UF(U3d, f)
% rho(k) = U(k) diag(f_k) U(k)^\dagger  (page-wise)
% U3d : [m x m x Nk]
% f   : [Nk x m]
    [m,~,Nk] = size(U3d);
    F3 = zeros(m,m,Nk);
    for ik=1:Nk
        F3(:,:,ik) = diag(f(ik,:));
    end
    rho_k = pagemtimes( pagemtimes(U3d, F3), permute(conj(U3d), [2 1 3]) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ Uniform k-mesh on a patch around K = (1/3,2/3) in frac ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = ...
    make_kmesh_valley_patch_fft(g, Nk1, Nk2, R_frac)
% MAKE_KMESH_VALLEY_PATCH_FFT
%   Uniform rectangular k-mesh on a *small patch* of the BZ,
%   centered at K = (1/3,2/3) in fractional coords.
%
%   Fractional coords region:
%       k1 ∈ [1/3 - R_frac, 1/3 + R_frac)
%       k2 ∈ [2/3 - R_frac, 2/3 + R_frac)
%
%   g.b : [3 x 3], rows are b1,b2,b3.
%
%   输出:
%     k_cart : [Nk x 3], Cartesian k (Nk = Nk1*Nk2)
%     k_frac : [Nk x 3], fractional k (k1,k2,k3=0)
%     w_k    : [Nk x 1], sum(w_k) = area_patch
%     Wnorm  : scalar, = area_patch

    % ---- 1. 整个 2D BZ 面积 ----
    area_BZ_2d = norm(cross(g.b(1,:), g.b(2,:)));

    % ---- 2. valley center in fractional coords ----
    K_frac = [1/3, 2/3];

    % ---- 3. local coordinates in [-R_frac, +R_frac) ----
    dk1 = ((0:Nk1-1)/Nk1 - 0.5) * 2*R_frac;
    dk2 = ((0:Nk2-1)/Nk2 - 0.5) * 2*R_frac;
    [D1,D2] = ndgrid(dk1, dk2);

    k1_patch = K_frac(1) + D1;
    k2_patch = K_frac(2) + D2;

    Nk = Nk1 * Nk2;
    k_frac        = zeros(Nk,3);
    k_frac(:,1)   = k1_patch(:);
    k_frac(:,2)   = k2_patch(:);
    k_frac(:,3)   = 0;

    k_cart = k_frac * g.b;

    % patch 面积 (frac-space) = (2R)^2
    area_patch_frac = (2*R_frac)^2;
    area_patch      = area_patch_frac * area_BZ_2d;

    % uniform weights
    w_k   = ones(Nk,1) * (area_patch / Nk);
    Wnorm = area_patch;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Precompute Vr_Fock(q) --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vr_Fock = precompute_Vr_Fock_uniform_fromB(Brows, m, Nk1, Nk2, opts)
% PRECOMPUTE_VR_FOCK_UNIFORM_FROMB
%   使用给定的 row-basis reciprocal 基底 Brows (通常是 g.b 或 patch 缩放后的 B_patch)，
%   在 (Nk1, Nk2) 的 q-grid 上构造 V_ab(q)，再做 FFT2 得到 Vr_Fock。
%
% 输入:
%   Brows : [d x d]，行是 b1, b2, ...（只用前两行）
%   m     : 轨道数 = size(g.wpos,1)
%   Nk1,Nk2 : q-grid 尺寸
%   opts  : 内含 opts.Vq_fun / opts.Vq_vec_fun, fock_drop_q0, qmin 等
%
% 输出:
%   Vr_Fock : [Nk1 x Nk2 x m x m] = FFT2( V_ab(q) )

    Nk = Nk1 * Nk2;

    % 构造 q_frac 网格: [0,1)×[0,1)，然后 wrap 到 [-0.5,0.5)
    q1v = (0:Nk1-1)/Nk1;
    q2v = (0:Nk2-1)/Nk2;
    [Q1,Q2] = ndgrid(q1v, q2v);         % (Nk1,Nk2)

    q_frac = zeros(Nk,3);
    q_frac(:,1) = Q1(:);
    q_frac(:,2) = Q2(:);
    q_frac(:,3) = 0;

    % wrap 到 [-0.5,0.5)
    q_frac = wrap_frac_pm(q_frac);

    % q_cart = q_frac * Brows
    q_cart = q_frac * Brows;              % (Nk x 3)

    % 用统一的 V(q) 模型计算 V_ab(q)
    if isfield(opts,'Vq_vec_fun') && ~isempty(opts.Vq_vec_fun)
        Vpages = opts.Vq_vec_fun(q_cart, opts);  % (m x m x Nk)
    else
        Vpages = zeros(m,m,Nk);
        for j = 1:Nk
            Vpages(:,:,j) = opts.Vq_fun(q_cart(j,:), opts);
        end
    end

    % 如果需要，drop q = 0：把 qmag 很小的那一点 V_ab(q) = 0
    if isfield(opts,'fock_drop_q0') && opts.fock_drop_q0
        qmag = vecnorm(q_cart,2,2);     % (Nk x 1)
        if isfield(opts,'qeps') && ~isempty(opts.qeps)
            qeps = opts.qeps;
        else
            qeps = 1e-8;
        end
        iz = (qmag < qeps);
        if any(iz)
            Vpages(:,:,iz) = 0;
        end
    end

    % reshape -> (Nk1,Nk2,m,m)，再 FFT2
    Vq4 = reshape(permute(Vpages, [3 1 2]), [Nk1, Nk2, m, m]);   % (Nk1,Nk2,m,m)
    Vr_Fock = fft2(Vq4);                                         % FFT(V(q))

    % 数值上保证每个 (a,b) matrix 是 Hermitian
    Vr_Fock = 0.5 * (Vr_Fock + permute(conj(Vr_Fock), [1 2 4 3]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Fock via FFT convolution -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2)
% FOCK_FFT_UNIFORM
%   FFT-based Fock self-energy on a uniform (patch) k-mesh.
%
%   输入:
%     rho_k    : [m x m x Nk]
%     w_k      : [Nk x 1] 均匀权重 (sum = Wnorm)
%     Wnorm    : = sum(w_k) = area_patch
%     Vr_Fock  : [Nk1 x Nk2 x m x m]  = FFT2(V(q))
%     Nk1,Nk2  : k-grid 尺寸, Nk = Nk1*Nk2
%
%   输出:
%     SigmaF_k : [m x m x Nk]

    [m, ~, Nk] = size(rho_k);
    assert(Nk == Nk1*Nk2, 'Nk mismatch with Nk1*Nk2');

    % reshape rho_k: (m,m,Nk) -> (Nk1,Nk2,m,m)
    rho4 = reshape(permute(rho_k, [3 1 2]), [Nk1, Nk2, m, m]);

    % 2D FFT over k1,k2 for each (a,b) channel
    rho_r = fft2(rho4);             % (Nk1,Nk2,m,m)

    % 权重因子:
    %   Σ_F(k) = -(1/Ω) ∑_{k'} w(k') V(k-k') ρ(k')
    %   对于均匀 patch: w_k/Wnorm = 1/Nk
    w_factor = w_k(1) / Wnorm;      % ideally = 1/Nk

    % 在 r-space 做逐点乘法：
    %   Σ_F(r) = - w_factor * V(r) .* ρ(r)，其中 V(r) = Vr_Fock = FFT(V(q))
    SigmaF_r = - w_factor * Vr_Fock .* rho_r;    % (Nk1,Nk2,m,m)

    % 回到 k 空间：ifft2(Σ_F(r)) 就是卷积
    SigmaF4 = ifft2(SigmaF_r);      % (Nk1,Nk2,m,m)

    % 数值上保证每页 Hermitian
    SigmaF4 = 0.5 * (SigmaF4 + permute(conj(SigmaF4), [1 2 4 3]));

    % reshape 回 (m,m,Nk)
    tmp       = reshape(SigmaF4, [Nk1*Nk2, m, m]);  % (Nk,m,m)
    SigmaF_k  = permute(tmp, [2 3 1]);              % (m,m,Nk)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------- Hartree term ------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SigmaH_k = hartree_uniform_smallq_patch(g, tau, rho_k, w_k, Wnorm, opts)
% HARTREE_UNIFORM_SMALLQ_PATCH
%   k-independent, diagonal Hartree shift with layer resolution.
%
%   使用 *完全相同* 的 V(q) 模型 (opts.Vq_fun)，在 q_cart → 0 的位置求 V_ab(q→0)。
%
%   公式:
%     nbar_a  = (1/Ω) ∑_k w(k) ρ_aa(k)
%     δn      = nbar - n_ref
%     Σ_H^a   = ∑_b V_ab(q→0) δn_b  (然后去掉整体平均避免 rigid shift)

    [m,~,Nk] = size(rho_k);
    w_k = w_k(:);

    % ---- 1) k-averaged orbital occupation per cell ----
    rho_flat = reshape(real(rho_k), m*m, Nk); % (m*m)×Nk
    idx_diag = 1:(m+1):m*m;
    rho_diag = rho_flat(idx_diag, :).';       % Nk×m
    nbar = (1/Wnorm) * (w_k.' * rho_diag).';  % m×1

    % ---- 2) reference background density n_ref ----
    if isfield(opts,'nref') && numel(opts.nref)==m
        nref = opts.nref(:);
    else
        nref = mean(nbar) * ones(m,1);        % 默认：整体平均
    end
    dn = nbar - nref;                         % m×1

    % ---- 3) V_ab(q→0)：统一通过 Vq_fun(q_cart≈0,opts) 计算 ----
    d_full = size(g.b,1);
    d = min(2,d_full);           % 2D 体系: 只取前两个方向
    Vsmall = zeros(m,m,d);
    eps_q = 1e-4;
    for iq=1:d
        qf = zeros(1,d_full);
        qf(iq) = eps_q;                    % tiny fractional step
        q_cart = qf * g.b;                % 1×d_full
        Vsmall(:,:,iq) = opts.Vq_fun(q_cart, opts);
    end
    Vab0 = mean(Vsmall(:,:,1:d), 3, 'omitnan');      % m×m

    % ---- 4) Diagonal Hartree shift Σ_H = V(q→0) * dn ----
    SigH_diag = Vab0 * dn;                    % m×1

    % 去掉整体平均，避免 rigid shift（只留层/轨道间相对势能差）
    SigH_diag = SigH_diag - mean(SigH_diag);

    D = diag(SigH_diag);                      % m×m
    SigmaH_k = repmat(D, [1 1 Nk]);          % 每个 k 上相同
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ H0(k) ------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H0_k, E_k, U3d] = build_H0k_from_wannier(g, k_frac, use_parfor)
% BUILD_H0K_FROM_WANNIER
%   Lattice-gauge TB:
%       H0(k) = ∑_{ΔR} t(ΔR) e^{ i k·ΔR }.
%
%   Row-basis convention:
%       k_cart = k_frac * g.b,  ΔR_cart = ΔR_frac * g.a
%
%   INPUT:
%     g.a      : [d x d], rows are a_i
%     g.b      : [d x d], rows are b_i
%     g.hopr   : [L x d], ΔR in fractional coords
%     g.ham    : [m x m x L], hopping t(ΔR)
%     k_frac   : [Nk x d], fractional k
%     use_parfor : logical
%
%   OUTPUT:
%     H0_k     : [m x m x Nk]
%     E_k      : [Nk x m]
%     U3d      : [m x m x Nk]

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
%   Hk: [m x m x Nk]
%   返回:
%     E_k : [Nk x m]    energies sorted ascending
%     U3d : [m x m x Nk] eigenvectors per k (columns are eigenstates)

    if nargin < 2, use_parfor = false; end

    [m,~,Nk] = size(Hk);
    E_k = zeros(Nk, m);
    U3d = zeros(m,m,Nk);

    if use_parfor
        parfor ik = 1:Nk
            H = (Hk(:,:,ik) + Hk(:,:,ik)')/2;     % Hermitize
            [U, d] = eig(H, 'vector');
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
%--------------------- μ & Fermi-Dirac occupations ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = fermi_dirac(E, mu, kT)
% FERMI_DIRAC  Occupations per k and band.
    if kT <= 0
        f = double(E <= mu);
    else
        x = (E - mu)/kT;
        f = 1 ./ (1 + exp(x));
    end
end

function mu = find_mu_patch(E_k, w_k, Wnorm, Nelec, kT, nspin)
% FIND_MU_PATCH
%   Bisection on μ so that:
%       nspin*(1/Ω)∑_k w(k)∑_n f_{kn} = Nelec.

    emin = min(E_k(:)) - 5*max(kT,1e-4);
    emax = max(E_k(:)) + 5*max(kT,1e-4);
    for it=1:80
        mu = 0.5*(emin+emax);
        f  = fermi_dirac(E_k, mu, kT);
        Ne = nspin * (1/Wnorm) * sum( w_k(:) .* sum(f,2) );
        if Ne > Nelec
            emax = mu;
        else
            emin = mu;
        end
        if abs(Ne-Nelec) < 1e-12*max(1,Nelec), break; end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Unified V(q) model (重要) ------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opts = make_Vq_model(opts, g)
% MAKE_VQ_MODEL
%   Build unified V(q) for both Fock and Hartree.
%
%   设置:
%     opts.Vq_fun(q_cart,opts)     -> [m x m]
%     opts.Vq_vec_fun(q_cart_all,opts) -> [m x m x Nk] (可选)
%
%   支持:
%     vmodel = '2d'        : V(q) = α/q * exp(-q |z_a - z_b|)
%     vmodel = 'keldysh'   : V(q) = α/[q(1+r0 q)] * exp(-q |z_a - z_b|)
%     vmodel = 'lattice'   : user-provided V_ab(R) real-space Fourier sum

    if ~isfield(opts,'layer_z') || isempty(opts.layer_z)
        opts.layer_z = g.wpos(:,3);
    end
    if ~isfield(opts,'vmodel'), opts.vmodel = '2d'; end

    switch lower(opts.vmodel)
        case '2d'
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;         % V(q)=alpha/q
            end
            opts.Vq_fun     = @(q_cart, opts2) ...
                Vq_2d_layered_core(q_cart, opts2.alpha, opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all, opts2) ...
                Vq_2d_layered_vec_core(q_cart_all, opts2.alpha, opts2.layer_z, opts2);

        case 'keldysh'
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;
            end
            if ~isfield(opts,'r0') || isempty(opts.r0)
                error('For vmodel=''keldysh'', please set opts.r0.');
            end
            opts.Vq_fun     = @(q_cart, opts2) ...
                Vq_keldysh_layered_core(q_cart, opts2.alpha, opts2.r0, ...
                                        opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all, opts2) ...
                Vq_keldysh_layered_vec_core(q_cart_all, opts2.alpha, opts2.r0, ...
                                            opts2.layer_z, opts2);

        case 'lattice'
            assert(isfield(opts,'VR_R_cart') && isfield(opts,'VR_ab'), ...
                'For vmodel=''lattice'', provide opts.VR_R_cart and opts.VR_ab.');
            opts.Vq_fun     = @(q_cart, opts2) ...
                Vq_from_VR_core(q_cart, opts2.VR_R_cart, opts2.VR_ab);
            opts.Vq_vec_fun = @(q_cart_all, opts2) ...
                Vq_from_VR_vec_core(q_cart_all, opts2.VR_R_cart, opts2.VR_ab);

        otherwise
            error('Unknown vmodel: %s', opts.vmodel);
    end
end

% ---- 2D Coulomb, layered (core) ----------------------------------------
function V_ab = Vq_2d_layered_core(q_cart, alpha, layer_z, opts)
    q = max(norm(q_cart), get_qmin(opts));
    DZ = abs(layer_z - layer_z.');
    V_ab = (alpha / q) * exp(-q * DZ);
end

% function V_ab = Vq_2d_layered_core(q_cart, alpha, layer_z, opts)
%     q = norm(q_cart);
%     if q < 1e-12
%         q = 1e-12;  % 数值防守，不作为物理 cutoff
%     end
%     DZ = abs(layer_z - layer_z.');
%     V_ab = (alpha / q) * exp(-q * DZ);
% end


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

% ---- Keldysh, layered (core) -------------------------------------------
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

% ---- Lattice Fourier V(R) (optional) -----------------------------------
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

% ---- small-q regulator --------------------------------------------------
function qmin = get_qmin(opts)
    if isfield(opts,'qmin') && ~isempty(opts.qmin)
        qmin = opts.qmin;
    else
        qmin = 1e-2;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ wrap_frac_pm -----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xw = wrap_frac_pm(x)
% wrap each component to [-1/2, 1/2)
    xw = x - round(x);
end
