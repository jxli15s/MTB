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

% clc;
% clear;
% g = MTB.geometry("Rgra_5s");
% g = MTB.read_poscar(g,"/Volumes/T9/work/tb/matlab/data/Graphene/3s/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('/Volumes/T9/work/tb/matlab/data/Graphene/3s/wannier90_hr_p1.dat','/Volumes/T9/work/tb/matlab/data/Graphene/3s/wannier90_hr_p2.dat');
% g.wpos=g.atoms*g.a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.wpos=g.atoms*g.a;
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
%%
% s3=[1  0
%     0  -1];
% laf=-diag([0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,-0.001,-0.001]);
% Zeeman=kron(s3,laf);
% g.add_zeeman(Zeeman)
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.9,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
hkpoints={[1/3,2/3,0.0]*0.9,...
          [1/3,2/3,0.0],...
          [1/3,2/3,0.0]+([0.5,0.5,0.0]-[1/3,2/3,0.0])*0.1,...
          };% hkpoints-high symmetry k points
nk=251;
efermi=-0.0;

Electric_field_in_evpA=0.00; %0.08-0.12 V/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
%%
%% ===================== 2. HF options =====================
opts = struct();

opts.kT     = 0;       % eV, effectively T≈0
opts.nspin  = 1;          % spin degeneracy
opts.Nelec  = nbands * opts.nspin / 2;   % half-filling, 2 band spinful

opts.mix    = 0.5;        % density mixing
opts.maxiter= 300;
opts.tol    = 1e-7;
opts.parfor = true;       % if you have Parallel Toolbox

% --- Coulomb interaction model ---
opts.vmodel = '2d';       % 2D Coulomb, V(q) = (alpha/q)*exp(-q|z_a-z_b|)
opts.e2     = 14.3996;    % e^2/(4π ε0) in eV*Å
opts.eps    = 27.0;        % relative dielectric constant
opts.r0     = 10.0;       % only used for 'keldysh'
opts.qmin=1e-4;
% layer z 取自 g.wpos
opts.layer_z = g.wpos(:,3);

% Fock 自能是否去掉 q~0 self-exchange
opts.fock_drop_q0 = true;
opts.qeps         = 1e-8;

% 可选：Hartree 的参考密度 (比如 charge neutrality)
% opts.nref = [1; 1];     % 每个轨道的 reference occupation (per cell)

% ===================== 3. Run HF on uniform BZ with FFT Fock =====================
Nk1 = 50;        % b1 方向 k 点数
Nk2 = 50;        % b2 方向 k 点数

out = main_hfmf_fft(g, opts, Nk1, Nk2);
%%
ek=reshape(out.E_k,Nk1,Nk2,[]);
figure()
surf(squeeze(ek(:,:,11)),'EdgeColor','none')
hold on;
surf(squeeze(ek(:,:,10)),'EdgeColor','none')
%%
ek=reshape(out.E_k,Nk1,Nk2,[]);
figure()
for i=9:12
    hold on;
    plot(out.E_k(:,i));
end
%%
enk=reshape(out.U3d,20,20,Nk1,Nk2);
kmesh=reshape(out.k_frac(:,1),Nk2,[]);
% enn=enk(:,:,63,40);

figure()
hold on
for i=9:12
    hold on;
surf(squeeze(ek(:,:,i)),'EdgeColor','none')
end
%%
function out = main_hfmf_fft(g, opts, Nk1, Nk2)
% MAIN_HFMF_FFT  Self-consistent HF on a UNIFORM k-mesh with FFT Fock.
% Row-basis convention: rows of g.b/g.a are b_i/a_i (Cartesian).
%
% Inputs:
%   g         : TB/Wannier model struct
%               g.a    [d x d], rows are a_i
%               g.b    [d x d], rows are b_i
%               g.hopr [L x d], hopping ΔR in fractional coords
%               g.ham  [m x m x L], hopping matrices t(ΔR)
%               g.wpos [m x 3], orbital positions (Cartesian)
%
%   opts      : same options as main_hfmf_patch
%               opts.kT, opts.nspin, opts.Nelec, opts.mix, opts.maxiter,
%               opts.tol, opts.parfor, opts.vmodel, opts.e2, opts.eps,
%               opts.r0, opts.qmin, etc.  (see make_Vq_model)
%
%   Nk1,Nk2   : uniform k-grid points along b1, b2 directions
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

    % ---------- defaults ----------
    if ~isfield(opts,'kT'),       opts.kT       = 0;        end
    if ~isfield(opts,'nspin'),    opts.nspin    = 2;        end
    if ~isfield(opts,'mix'),      opts.mix      = 0.5;      end
    if ~isfield(opts,'maxiter'),  opts.maxiter  = 50;       end
    if ~isfield(opts,'tol'),      opts.tol      = 1e-8;     end
    if ~isfield(opts,'parfor'),   opts.parfor   = false;    end
    if ~isfield(opts,'fock_drop_q0'), opts.fock_drop_q0 = false; end

    % ---------- UNIFORM k-mesh on full BZ (rectangular in frac coords) ----------
    R_frac = 0.15;       % patch 只有 BZ 的 5% 宽度（非常小）
    [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = make_kmesh_valley_patch_fft(g, Nk1, Nk2, R_frac);

    % [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = make_kmesh_uniform_rect_fft(g, Nk1, Nk2);
    Nk = size(k_cart,1);   % Nk = Nk1 * Nk2

    % ---------- H0(k) from Wannier TB (row-basis aware,与你现有的一致) ----------
    [H0_k, E0_k, U0_3d] = build_H0k_from_wannier(g, k_frac, opts.parfor);
    m   = size(g.ham,1);
    kT  = opts.kT;
    nsp = opts.nspin;

    % ---------- Unified V(q) model for BOTH Fock & Hartree ----------
    opts = make_Vq_model(opts, g);   % sets opts.Vq_fun / Vq_vec_fun

    % ---------- Initial density from non-interacting bands ----------
    if ~isfield(opts,'Nelec') || isempty(opts.Nelec)
        opts.Nelec = m * nsp / 2;    % 默认半填
    end
    mu  = find_mu_patch(E0_k, w_k, Wnorm, opts.Nelec, kT, nsp);
    f0  = fermi_dirac(E0_k, mu, kT);      % [Nk x m]
    rho_k   = rho_from_UF(U0_3d, f0);     % [m x m x Nk], per-cell
    rho_orb = kavg_rho(rho_k, w_k, Wnorm);

    % ---------- Precompute Vr_Fock(q) for FFT Fock ----------
    Vr_Fock = precompute_Vr_Fock_uniform(g, Nk1, Nk2, opts);

    history = [];
    % 为了 debug，先初始化一下（避免没收敛时 out 变量没定义）
    SigmaF_k = zeros(m,m,Nk);
    SigmaH_k = zeros(m,m,Nk);
    Hmf_k    = H0_k;
    E_k      = E0_k;
    U3d      = U0_3d;

    % ---------- SCF loop ----------
    for it = 1:opts.maxiter
        % 1) Fock Σ_F(k): FFT-based convolution on uniform BZ
        SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2);

        % 2) Hartree Σ_H(k): k-independent, diagonal; SAME V(q) model at q→0
        SigmaH_k = hartree_uniform_smallq_patch(g, g.wpos, rho_k, w_k, Wnorm, opts);

        % 3) Total mean-field Hamiltonian
        Hmf_k = H0_k + SigmaF_k + SigmaH_k;

        % 4) Diagonalize (optionally parfor)
        [E_k, U3d] = diag_all_k(Hmf_k, opts.parfor);

        % 5) New chemical potential and occupations
        mu = find_mu_patch(E_k, w_k, Wnorm, opts.Nelec, kT, nsp);
        f  = fermi_dirac(E_k, mu, kT);

        % 6) New density
        rho_k_new   = rho_from_UF(U3d, f);
        rho_orb_new = kavg_rho(rho_k_new, w_k, Wnorm);

        % 7) Density mixing (linear)
        rho_k   = (1-opts.mix)*rho_k   + opts.mix*rho_k_new;
        rho_orb = (1-opts.mix)*rho_orb + opts.mix*rho_orb_new;

        % 8) Convergence monitor
        err = norm(rho_orb_new - rho_orb, 'fro') / max(1, norm(rho_orb,'fro'));
        history = [history; it, mu, err]; %#ok<AGROW>
        fprintf('SCF(FFT) %3d : mu = %.6f, err = %.3e\n', it, mu, err);

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
    % rho = U * F * U^\dagger   (page-wise)
    rho_k = pagemtimes( pagemtimes(U3d, F3), permute(conj(U3d), [2 1 3]) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------- Uniform k-mesh ------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = make_kmesh_uniform_rect_fft(g, Nk1, Nk2)
% MAKE_KMESH_UNIFORM_RECT_FFT
% Uniform rectangular k-mesh over the full BZ parallelogram.
% Fractional coords in [0,1) x [0,1), k_cart = k_frac * g.b (row-basis).
%
% g.b: [3 x 3], each row is a reciprocal basis vector (b1; b2; b3).

    % 仅用前两行计算 2D BZ 面积
    scale=0.4;
    area_BZ_2d = norm(cross(g.b(1,:),g.b(2,:)));      % BZ 平行四边形面积

    % uniform fractional grid in [0,1)
    k1v = (0:Nk1-1) / Nk1;
    k2v = (0:Nk2-1) / Nk2;
    [K1,K2] = ndgrid(k1v, k2v);     % (Nk1,Nk2)

    Nk  = Nk1*Nk2;
    k_frac        = zeros(Nk,3);
    k_frac(:,1)   = K1(:);
    k_frac(:,2)   = K2(:);
    k_frac(:,3)   = 0;              % kz = 0 for 2D
    k_frac=k_frac.*scale-[0.5,0.5,0]+[1/3,2/3,0.0];

    % cartesian k vectors: (Nk×3)*(3×3) = (Nk×3)
    k_cart = k_frac * g.b;          % row-basis

    % uniform weights: sum(w_k) = area_BZ
    w_k   = ones(Nk,1) * (area_BZ_2d / Nk);
    Wnorm = area_BZ_2d;               % = area_BZ_2d
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ Uniform k-mesh on a small patch around (1/3,2/3) in frac ---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = ...
    make_kmesh_valley_patch_fft(g, Nk1, Nk2, R_frac)
% MAKE_KMESH_VALLEY_PATCH_FFT
%   Uniform rectangular k-mesh on a *small patch* of the BZ,
%   centered at K = (1/3,2/3) in fractional coords.
%
%   Fractional coords region:
%       k1 \in [1/3 - R_frac, 1/3 + R_frac)
%       k2 \in [2/3 - R_frac, 2/3 + R_frac)
%   Then folded back into [0,1) via mod 1.
%
% Inputs:
%   g.b   : [3 x 3], rows are b1,b2,b3 (row-basis, Cartesian)
%   Nk1   : number of k-points along "b1" direction in the patch
%   Nk2   : number of k-points along "b2" direction in the patch
%   R_frac: half-width of the patch in fractional coords, e.g. 0.05
%
% Outputs:
%   k_cart : [Nk x 3], Cartesian k (Nk = Nk1*Nk2)
%   k_frac : [Nk x 3], fractional k (k1,k2,k3=0) in [0,1)
%   w_k    : [Nk x 1], integration weights, sum(w_k) = area_patch
%   Wnorm  : scalar, = area_patch
%   Nk1,Nk2: returned as given
%
% NOTE:
%   - If you want to treat this patch as an *effective BZ* for valley theory,
%     use Wnorm = area_patch, and the Coulomb V(q) 应该也用这个 patch 的 q。

    % ---- 1. BZ 基本信息：整 BZ 的面积 (仅用 b1,b2 的 2D 部分) ----
    % b1 = g.b(1,1:2).';
    % b2 = g.b(2,1:2).';
    % B2 = [b1, b2];                    % 2x2
    % area_BZ_2d = abs(det(B2));        % 整个 BZ 的平行四边形面积
    area_BZ_2d = norm(cross(g.b(1,:),g.b(2,:))); 


    % ---- 2. valley 中心 (K 点) 的分数坐标 ----
    K_frac = [1/3, 2/3];              % 你要的 (1/3,2/3)
    % 也可以改成 K' = [2/3,1/3] 等

    % ---- 3. 在 [-R_frac, +R_frac) 上均匀取 Nk1,Nk2 个点 ----
    % 这里先构造一个"局部坐标" delta_k1, delta_k2
    dk1 = ( (0:Nk1-1)/Nk1 - 0.5 ) * 2*R_frac;   % in [-R_frac, +R_frac)
    dk2 = ( (0:Nk2-1)/Nk2 - 0.5 ) * 2*R_frac;   % in [-R_frac, +R_frac)
    [D1,D2] = ndgrid(dk1, dk2);                 % (Nk1,Nk2)

    % ---- 4. shift 到以 K_frac 为中心的 patch，并折回 [0,1) ----
    k1_patch = K_frac(1) + D1;   % frac coords around K
    k2_patch = K_frac(2) + D2;

    % % 折回 [0,1) (mod 1)
    % k1_patch = k1_patch - floor(k1_patch);
    % k2_patch = k2_patch - floor(k2_patch);

    % ---- 5. 组装成 k_frac (Nk x 3) ----
    Nk = Nk1*Nk2;
    k_frac        = zeros(Nk,3);
    k_frac(:,1)   = k1_patch(:);
    k_frac(:,2)   = k2_patch(:);
    k_frac(:,3)   = 0;         % 2D -> k3 = 0

    % ---- 6. Cartesian k: k_cart = k_frac * g.b ----
    k_cart = k_frac * g.b;     % (Nk x 3) * (3 x 3)

    % ---- 7. patch 的面积 & 积分权重 ----
    % 整个 BZ 对应 frac 区域的面积是 1.0
    % 我们现在取的是 frac 区域:
    %   [K1-R_frac, K1+R_frac) x [K2-R_frac, K2+R_frac)
    % 面积是 (2R_frac)^2 (在 frac 坐标里)
    area_patch_frac = (2*R_frac)^2;          % frac-space area
    area_patch      = area_patch_frac * area_BZ_2d;   % 真实 k-space 面积

    % 均匀权重
    w_k   = ones(Nk,1) * (area_patch / Nk);
    Wnorm = area_patch;
    % Wnorm = area_BZ_2d;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Precompute Vr_Fock(q) --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vr_Fock = precompute_Vr_Fock_uniform(g, Nk1, Nk2, opts)
% PRECOMPUTE_VR_FOCK_UNIFORM
% Build V_ab(q) on uniform (Nk1,Nk2) BZ grid, then FFT2→Vr for Fock.
%
% 返回:
%   Vr_Fock : (Nk1, Nk2, m, m)   频域核 (already FFT(V(q))),
%             供 fock_fft_uniform 使用

    m  = size(g.wpos,1);
    Nk = Nk1 * Nk2;

    % 构造 q_frac 网格: [0,1)×[0,1)，再 wrap 到 [-0.5,0.5) 保持和之前 q 处理一致
    q1v = (0:Nk1-1)/Nk1;
    q2v = (0:Nk2-1)/Nk2;
    [Q1,Q2] = ndgrid(q1v, q2v);         % (Nk1,Nk2)

    q_frac = zeros(Nk,3);
    q_frac(:,1) = Q1(:);
    q_frac(:,2) = Q2(:);
    q_frac(:,3) = 0;

    q_frac = wrap_frac_pm(q_frac);      % 每个分量到 [-0.5,0.5)

    % q_cart = q_frac * g.b (row-basis)
    q_cart = q_frac * g.b;              % (Nk x 3)

    % 用统一的 Vq_vec_fun 或 Vq_fun 计算 V_ab(q)
    if isfield(opts,'Vq_vec_fun') && ~isempty(opts.Vq_vec_fun)
        Vpages = opts.Vq_vec_fun(q_cart, opts);  % (m x m x Nk)
    else
        Vpages = zeros(m,m,Nk);
        for j = 1:Nk
            Vpages(:,:,j) = opts.Vq_fun(q_cart(j,:), opts);
        end
    end

    % 如果需要，去掉 q ~ 0 的 self-exchange
    if isfield(opts,'fock_drop_q0') && opts.fock_drop_q0
        qmag = vecnorm(q_cart,2,2);     % (Nk x 1)
        % 阈值：取 max(qmin, qeps)，避免数值问题
        if isfield(opts,'qeps') && ~isempty(opts.qeps)
            qeps = opts.qeps;
        else
            qeps = 1e-6;
        end
        z_th = max(qeps, 1e-8);
        iz = (qmag < z_th);
        if any(iz)
            Vpages(:,:,iz) = 0;
        end
    end

    % reshape 到 (Nk1,Nk2,m,m)，再做一次 FFT2
    Vq4 = reshape(permute(Vpages, [3 1 2]), [Nk1, Nk2, m, m]);   % (Nk1,Nk2,m,m)
    Vr_Fock = fft2(Vq4);                                         % F(V(q))

    % 可选：保证每个 q 对应的矩阵 Hermitian（数值上）
    Vr_Fock = 0.5 * (Vr_Fock + permute(conj(Vr_Fock), [1 2 4 3]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Fock via FFT convolution -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2)
% FOCK_FFT_UNIFORM
% FFT-based Fock self-energy on a uniform BZ.
%
% 输入:
%   rho_k    : [m x m x Nk]
%   w_k      : [Nk x 1] 均匀权重 (sum = Wnorm)
%   Wnorm    : = sum(w_k) = area_BZ
%   Vr_Fock  : [Nk1 x Nk2 x m x m]  = fft2(V(q))
%   Nk1,Nk2  : k-grid 尺寸
%
% 输出:
%   SigmaF_k : [m x m x Nk]

    [m, ~, Nk] = size(rho_k);
    assert(Nk == Nk1*Nk2, 'Nk mismatch with Nk1*Nk2');

    % reshape rho_k: (m,m,Nk) -> (Nk1,Nk2,m,m)
    rho4 = reshape(permute(rho_k, [3 1 2]), [Nk1, Nk2, m, m]);

    % 2D FFT over k1,k2 for each (a,b) channel
    rho_r = fft2(rho4);             % (Nk1,Nk2,m,m)

    % 权重因子: w_k/Wnorm = 1/(Nk1*Nk2) (均匀 BZ 情况)
    w_factor = w_k(1) / Wnorm;      % 理论上 = 1/Nk

    % Σ_F(r) = - (w_factor) * V(r) .* ρ(r)，其中 V(r) = Vr_Fock = FFT(V(q))
    % SigmaF_r = - w_factor * Vr_Fock .* rho_r;
    SigmaF_r = - w_factor*Vr_Fock .* rho_r;

    % 回到 k 空间：ifft2 给出圆卷积 conv_k
    SigmaF4 = ifft2(SigmaF_r);      % (Nk1,Nk2,m,m)

    % 强制 Hermitian（数值对称化）
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
    dn = nbar - nref;                         % m×1

    % V_ab(q→0): sample tiny q along each reciprocal direction (row-basis)
    d = size(g.b,1);
    Vsmall = zeros(m,m,d);
    eps_q = 1e-4;
    for iq=1:d
        qf = zeros(1,d); qf(iq) = eps_q;      % tiny fractional step
        q_cart = qf * g.b;                    % 1×d
        Vsmall(:,:,iq) = opts.Vq_fun(q_cart, opts);  % SAME model as Fock
    end
    Vab0 = mean(Vsmall(:,:,1:d), 3, 'omitnan');      % m×m

    % Diagonal Hartree shift; remove common mode to avoid trivial rigid shift
    SigH_diag = Vab0 * dn;                    % m×1
    SigH_diag = SigH_diag - mean(SigH_diag);

    D = diag(SigH_diag);                      % m×m
    SigmaH_k = repmat(D, [1 1 Nk]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ H0(k) ------------------------------------%
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
            % dR_cart = (g.a * dR_frac.').';
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
%--------------------- μ & Fermi-Dirac occupations ----------------------%
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
%---------------------- Unified V(q) model (你的版本) --------------------%
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
    if isfield(opts,'qmin') && ~isempty(opts.qmin)
        qmin = opts.qmin;
    else
        qmin = 1e-6;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ wrap_frac_pm -----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xw = wrap_frac_pm(x)
% wrap each component to [-1/2, 1/2)
    xw = x - round(x);
end
