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

% s3=[1  0
%     0  -1];
% laf=-diag([0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,-0.001,-0.001])*100;
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
kpoint=[1/3,2/3,0.0];
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
% kpoint=[-1.2/3,-2/3,0.0];
% [Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- 设定高对称线 (Γ -> K) 的 frac 坐标 ---------
Nk_line = 200;

Gamma_frac = [0,   0];
K_frac     = [2/3, 1/3];

t = linspace(0, 1, Nk_line).';   % 参数 0 → 1
k_line_frac = (1 - t).*Gamma_frac + t.*K_frac;   % [Nk_line x 2]

% --------- 调用 valley_block_line ---------
out_line = valley_block_line(g, k_line_frac);

% out_line.E_valley: [Nk_line x 2m] 的本征值，可以用来画能带
E_valley = out_line.E_valley;

% --------- 示例：画 valley-block 的能带 ---------
figure;
plot(t, E_valley, 'k-');  % 简单全部画出来
xlabel('t along \Gamma \rightarrow K');
ylabel('Energy (eV)');
title('Valley-block band structure along \Gamma-K: diag(H(k), H(-k))');
grid on;
%%
% 1. 定义三个拐点（你这段 hkpoints）
hkpoints = {
    [1/3, 2/3, 0.0] * 0.9, ...
    [1/3, 2/3, 0.0], ...
    [1/3, 2/3, 0.0] + ([0.5,0.5,0.0] - [1/3,2/3,0.0]) * 0.1 ...
};

% 2. 沿每一段做线性插值，拼成一条完整的 k_line_frac
Nk_seg = 50;  % 每一小段取 50 个点
k_line_frac = [];

for s = 1:(numel(hkpoints)-1)
    k0 = hkpoints{s};
    k1 = hkpoints{s+1};

    t  = linspace(0, 1, Nk_seg+1)';  % 0→1
    seg = (1-t).*k0 + t.*k1;         % 线性插值

    if s > 1
        seg = seg(2:end,:);          % 去掉与上一段重复的端点
    end
    k_line_frac = [k_line_frac; seg]; %#ok<AGROW>
end

% 3. 调用 valley_block_line，得到 [H(k),0;0,H(-k)]
out_line = valley_block_line(g, k_line_frac);

% 4. 画 valley-block 能带
E_valley = out_line.E_valley;  % [Nk_total x 2m]
figure;
plot(E_valley, 'k-');
xlabel('k-point index along line');
ylabel('Energy (eV)');
title('Valley-block band structure around K');
grid on;
hold on;
plot(E_valley(:,17:20),'r--')
%%
%-------------------------------------------
% 0. 准备 Wannier TB 结构 g
%   (你已经有了, 比如 load('g_TaIrTe4.mat') 得到 g)
%-------------------------------------------

Nk1 = 40;    % k-grid along b1 in K-patch
Nk2 = 40;    % k-grid along b2 in K-patch

opts = struct();

% ----- HF 参数 -----
opts.kT      = 0;
opts.nspin   = 1;
opts.mix     = 0.5;
opts.maxiter = 80;
opts.tol     = 1e-6;
opts.parfor  = true;

% ----- 电子数: double-valley half filling -----
m0 = size(g.ham,1);   % single valley
m  = 2*m0;            % double valley
opts.Nelec = m * opts.nspin / 2;

% ----- 相互作用模型: 比如 double-gate -----
opts.vmodel = 'doublegate';
opts.e2     = 14.3996;    % eV·Å
opts.eps    = 10.0;        % 有效介电常数
opts.Dgate  = 300.0;       % Å

% small-q & Fock q=0
opts.qmin        = 1e-4;   % regulator in |q|
opts.fock_drop_q0 = true;  % 丢掉 Fock 的 q=0

% (可选) Hartree 参考密度: 比如 neutrality
% opts.nref = something;  % [2*m0 x 1]，默认用平均值

%-------------------------------------------
% 1. 跑 double-valley HFMF + FFT Fock
%-------------------------------------------
out = main_hfmf_fft(g, opts, Nk1, Nk2);

% out 里包含:
%   out.E_k      [Nk x (2*m0)]  HF 本征能
%   out.U3d      [2*m0 x 2*m0 x Nk]  本征矢
%   out.mu       标量化学势
%   out.SigmaF_k, out.SigmaH_k
%   out.rho_k, out.rho_orb
%   out.k_frac, out.k_cart, out.w_k, etc.

%-------------------------------------------
% 2. 如何区分 valley
%-------------------------------------------
% m0 = size(g.ham,1);
% 对任意 k：
%   - valley K 的轨道是 1:m0
%   - valley K' 的轨道是 m0+1:2*m0

% 例如: 取某个 k 点 ik 看两个 valley 的能谱
ik = 1;
E_k_valleyK  = out.E_k(ik, 1:m0);
E_k_valleyKp = out.E_k(ik, m0+1:2*m0);


% ==== 4. 简单看一下结果 ====
figure;
plot(out.history(:,1), out.history(:,3), '-o');
xlabel('SCF iter'); ylabel('err'); title('SCF convergence');

fprintf('Final mu = %.6f eV, final err = %.3e\n', out.mu, out.history(end,3));


%%
ek=reshape(out.E_k,Nk1,Nk2,[]);
figure()
for i=17:24
    hold on;
    plot(out.E_k(:,i));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- k-average & density builder ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = main_hfmf_fft(g, opts, Nk1, Nk2)
% MAIN_HFMF_FFT (double-valley version)
%   Self-consistent Hartree–Fock on a UNIFORM k-mesh (valley patch)
%   using FFT-based Fock term, with an explicit double valley:
%
%       H_dv(k) = blkdiag( H(k), H(-k) )
%
%   Row-basis convention:
%       g.a, g.b are [d x d], each row is a_i / b_i (Cartesian).
%
%   INPUT:
%     g.a      [d x d]  real-space primitive vectors (rows)
%     g.b      [d x d]  reciprocal primitive vectors (rows)
%     g.hopr   [L x d]  hopping ΔR in fractional coords
%     g.ham    [m0 x m0 x L]  hopping matrices t(ΔR)  (single valley)
%     g.wpos   [m0 x 3]  orbital positions (Cartesian, z 包含层信息)
%
%     opts.kT       : temperature (energy units)
%     opts.nspin    : spin degeneracy
%     opts.Nelec    : target electrons per cell (incl. spin & both valleys)
%     opts.mix      : density mixing (0~0.8)
%     opts.maxiter  : SCF max iterations
%     opts.tol      : convergence tolerance on rho_orb
%     opts.parfor   : logical, parallelize diagonalization
%
%     Interaction model (统一给 Fock + Hartree):
%       opts.vmodel  : '2d' | 'keldysh' | 'doublegate' | 'lattice'
%       opts.e2      : e^2/(4π ε0) in your units (e.g., 14.3996 eV·Å)
%       opts.eps     : ε_r (relative dielectric)
%       opts.r0      : Keldysh screening length (for 'keldysh'), length units
%       opts.Dgate   : gate separation parameter D (for 'doublegate'), same units as z
%       opts.layer_z : optional, [m x 1] layer z; 这里会被覆盖成 double-valley 的
%
%       可选数值参数:
%         opts.qmin        : small-q regulator (默认 1e-2 in |q|)
%         opts.fock_drop_q0: logical, 是否在 Fock 中丢掉 q=0 自交换 (默认 false)
%         opts.nref        : [m x 1] reference density for Hartree 背景
%
%   Nk1, Nk2 : k-grid 点数 (valley patch 内部的矩形网格)
%
%   NOTE:
%     - k-mesh 仍然是 K=(1/3,2/3) 附近的小 patch：
%           make_kmesh_valley_patch_fft(...)
%
%   OUTPUT (out):
%     .k_cart, .k_frac, .w_k, .Wnorm
%     .H0_k, .Hmf_k, .SigmaF_k, .SigmaH_k   [m x m x Nk], m = 2*m0
%     .E_k                                  [Nk x m]
%     .U3d                                  [m x m x Nk]
%     .mu
%     .rho_k                                [m x m x Nk]
%     .rho_orb                              [m x m]
%     .history   [iter, mu, err]
%     .Nk1, .Nk2

    % ---------- defaults ----------
    if ~isfield(opts,'kT'),            opts.kT       = 0;        end
    if ~isfield(opts,'nspin'),         opts.nspin    = 2;        end
    if ~isfield(opts,'mix'),           opts.mix      = 0.5;      end
    if ~isfield(opts,'maxiter'),       opts.maxiter  = 50;       end
    if ~isfield(opts,'tol'),           opts.tol      = 1e-8;     end
    if ~isfield(opts,'parfor'),        opts.parfor   = false;    end
    if ~isfield(opts,'fock_drop_q0'),  opts.fock_drop_q0 = false;end
    if ~isfield(opts,'vmodel'),        opts.vmodel   = '2d';     end

    % ---------- k-mesh: valley patch around K = (1/3,2/3) ----------
    R_frac = 0.08;   % patch 半宽 (fractional)，即 [K1±R, K2±R]
    [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = ...
        make_kmesh_valley_patch_fft(g, Nk1, Nk2, R_frac);
    Nk = size(k_cart,1);

    % ---------- single-valley H0(k) and H0(-k) ----------
    % H_plus(k)  = H(k)
    % H_minus(k) = H(-k)  (注意 frac 坐标要 wrap 一下)
    [H0_k_plus,  ~,  ~] = build_H0k_from_wannier(g, k_frac, opts.parfor);

    k_frac_minus = -k_frac;
    k_frac_minus = wrap_frac_pm(k_frac_minus);  % frac ∈ [-1/2,1/2) 再进 H(k)
    [H0_k_minus, ~,  ~] = build_H0k_from_wannier(g, k_frac_minus, opts.parfor);

    % ---------- build double-valley H0_dv(k) = blkdiag(H(k), H(-k)) ----------
    m0 = size(g.ham,1);      % single-valley orbital count
    m  = 2*m0;               % double-valley orbital count

    H0_k = zeros(m,m,Nk);
    for ik = 1:Nk
        H0_k(:,:,ik) = blkdiag(H0_k_plus(:,:,ik), H0_k_minus(:,:,ik));
    end

    % 对 double-valley H0(k) 对角化，得到 E0_k, U0_3d
    [E0_k, U0_3d] = diag_all_k(H0_k, opts.parfor);

    kT  = opts.kT;
    nsp = opts.nspin;

    % ---------- electrons per cell (for double valley) ----------
    if ~isfield(opts,'Nelec') || isempty(opts.Nelec)
        % 默认：double-valley + spin 的 half-filling
        opts.Nelec = m * nsp / 2;
    end

    % ---------- small-q regulator qmin ----------
    if ~isfield(opts,'qmin') || isempty(opts.qmin)
        opts.qmin = 1e-2;   % in |q| units
    end

    % ---------- construct double-valley geometry for interactions ----------
    % 只为了相互作用：把 wpos 按 valley 复制一份，层信息从 z 读出。
    wpos_dv = [g.wpos; g.wpos];   % [2*m0 x 3]
    g_dv = struct();
    g_dv.a    = g.a;
    g_dv.b    = g.b;
    g_dv.wpos = wpos_dv;

    % 强制用 double-valley 的 z 作为 layer_z
    opts.layer_z = wpos_dv(:,3);

    % ---------- Unified V(q) model for BOTH Fock & Hartree (double valley) ----------
    opts = make_Vq_model(opts, g_dv);   % 设置 opts.Vq_fun / Vq_vec_fun 等，基于 2*m0

    % ---------- initial density from non-interacting double-valley bands ----------
    mu  = find_mu_patch(E0_k, w_k, Wnorm, opts.Nelec, kT, nsp);
    f0  = fermi_dirac(E0_k, mu, kT);      % [Nk x m]
    rho_k   = rho_from_UF(U0_3d, f0);     % [m x m x Nk], double valley
    rho_orb = kavg_rho(rho_k, w_k, Wnorm);

    % ---------- Precompute V(q) -> Vr_Fock (FFT kernel, double valley) ----------
    Vr_Fock = precompute_Vr_Fock_uniform(g_dv, Nk1, Nk2, opts);

    history = [];
    % 初始化输出变量
    SigmaF_k = zeros(m,m,Nk);
    SigmaH_k = zeros(m,m,Nk);
    Hmf_k    = H0_k;
    E_k      = E0_k;
    U3d      = U0_3d;

    % ---------- SCF loop ----------
    for it = 1:opts.maxiter
        % 1) Fock Σ_F(k): FFT-based convolution (uniform k-grid)
        SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2);

        % 2) Hartree Σ_H(k): k-independent diagonal term,
        %    使用统一 V(q) 模型的 q->0 极限 (double valley 尺度)
        SigmaH_k = hartree_uniform_smallq(g_dv, g_dv.wpos, rho_k, w_k, Wnorm, opts);

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
        fprintf('SCF(FFT, 2valley) %3d : mu = %.6f, err = %.3e\n', it, mu, err);

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

function rho_orb = kavg_rho(rho_k, w_k, Wnorm)
% K-AVERAGE of rho(k) to per-cell density matrix:
%   ρ_orb = (1/Ω_BZ) ∑_k w(k) ρ(k)
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
    % 使用 pagemtimes 做 page-wise 乘法：
    %   ρ_k = U * F * U^\dagger
    rho_k = pagemtimes( pagemtimes(U3d, F3), permute(conj(U3d), [2 1 3]) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- k-mesh: full BZ (optional) -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_cart, k_frac, w_k, Wnorm, Nk1, Nk2] = ...
    make_kmesh_uniform_rect_fft(g, Nk1, Nk2)
% MAKE_KMESH_UNIFORM_RECT_FFT
%   Uniform rectangular k-mesh over the full BZ parallelogram.
%   Fractional coords in [0,1) x [0,1), k_cart = k_frac * g.b (row-basis).
%
%   g.b: [3 x 3], rows are reciprocal basis vectors (b1; b2; b3).

    % BZ area (2D)：用前三个分量的叉乘
    area_BZ_2d = norm(cross(g.b(1,:), g.b(2,:)));

    % uniform fractional grid in [0,1)
    k1v = (0:Nk1-1) / Nk1;
    k2v = (0:Nk2-1) / Nk2;
    [K1,K2] = ndgrid(k1v, k2v);     % (Nk1,Nk2)

    Nk  = Nk1*Nk2;
    k_frac        = zeros(Nk,3);
    k_frac(:,1)   = K1(:);
    k_frac(:,2)   = K2(:);
    k_frac(:,3)   = 0;              % kz = 0 for 2D

    % cartesian k vectors: (Nk×3)*(3×3) = (Nk×3)
    k_cart = k_frac * g.b;          % row-basis

    % uniform weights: sum(w_k) = area_BZ
    w_k   = ones(Nk,1) * (area_BZ_2d / Nk);
    Wnorm = area_BZ_2d;
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
%   INPUT:
%     g.b   : [3 x 3], rows are b1,b2,b3 (row-basis, Cartesian)
%     Nk1   : #k-points along b1 direction in the patch
%     Nk2   : #k-points along b2 direction in the patch
%     R_frac: half-width of the patch in fractional coords (e.g. 0.05–0.2)
%
%   OUTPUT:
%     k_cart : [Nk x 3], Cartesian k (Nk = Nk1*Nk2)
%     k_frac : [Nk x 3], fractional k (k1,k2,k3=0)
%     w_k    : [Nk x 1], integration weights, sum(w_k) = area_patch
%     Wnorm  : scalar, = area_patch
%     Nk1,Nk2: as input

    % ---- 1. BZ area (2D) ----
    area_BZ_2d = norm(cross(g.b(1,:), g.b(2,:)));

    % ---- 2. valley center in fractional coords ----
    K_frac = [1/3, 2/3];  % K valley; K' would be [2/3,1/3]

    % ---- 3. local coordinates in [-R_frac, +R_frac) ----
    dk1 = ((0:Nk1-1)/Nk1 - 0.5) * 2*R_frac;   % in [-R, +R)
    dk2 = ((0:Nk2-1)/Nk2 - 0.5) * 2*R_frac;
    [D1,D2] = ndgrid(dk1, dk2);               % (Nk1,Nk2)

    k1_patch = K_frac(1) + D1;   % frac coords around K
    k2_patch = K_frac(2) + D2;

    % 注意：这里不再 mod 1，patch 本身就在 BZ 内部。
    Nk = Nk1 * Nk2;
    k_frac        = zeros(Nk,3);
    k_frac(:,1)   = k1_patch(:);
    k_frac(:,2)   = k2_patch(:);
    k_frac(:,3)   = 0;

    % Cartesian k: k_cart = k_frac * g.b
    k_cart = k_frac * g.b;

    % patch area in frac-space = (2R)^2
    area_patch_frac = (2*R_frac)^2;
    area_patch      = area_patch_frac * area_BZ_2d;

    % uniform weights on patch
    w_k   = ones(Nk,1) * (area_patch / Nk);
    Wnorm = area_patch;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Precompute Vr_Fock(q) --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vr_Fock = precompute_Vr_Fock_uniform(g, Nk1, Nk2, opts)
% PRECOMPUTE_VR_FOCK_UNIFORM
%   Build V_ab(q) on a uniform (Nk1,Nk2) q-grid (periodic BZ),
%   then FFT2 -> Vr_Fock, which is the convolution kernel in r-space.
%
%   返回:
%     Vr_Fock : (Nk1, Nk2, m, m)  = FFT2( V_ab(q) )
%               fock_fft_uniform 中会用:
%                   Σ_F(r) = - (w_factor) * Vr_Fock .* ρ(r)
%               然后再 ifft2 回到 k 空间。
%
%   这里可以通过 opts.fock_drop_q0 把 V(q=0) 设为 0，从而"丢掉自交换"。
%   虽然是 FFT，但只要在 q-grid 上把 q=0 对应的那一点清零就相当于 drop_q0。

    m  = size(g.wpos,1);
    Nk = Nk1 * Nk2;

    % 构造 q_frac 网格: [0,1)×[0,1)，再 wrap 到 [-0.5,0.5) 保持对称
    q1v = (0:Nk1-1)/Nk1;
    q2v = (0:Nk2-1)/Nk2;
    [Q1,Q2] = ndgrid(q1v, q2v);         % (Nk1,Nk2)

    q_frac = zeros(Nk,3);
    q_frac(:,1) = Q1(:);
    q_frac(:,2) = Q2(:);
    q_frac(:,3) = 0;

    % wrap 到 [-0.5,0.5)（只影响 q 的"label"，不影响 FFT 结构）
    q_frac = wrap_frac_pm(q_frac);

    % q_cart = q_frac * g.b
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

    % 如果需要，drop q = 0：把 qmag 很小的那个 grid 点的 V_ab(q) = 0
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

    % 数值上保证 kernel 每个 (q,r) 对应矩阵是 Hermitian (可选)
    Vr_Fock = 0.5 * (Vr_Fock + permute(conj(Vr_Fock), [1 2 4 3]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Fock via FFT convolution -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SigmaF_k = fock_fft_uniform(rho_k, w_k, Wnorm, Vr_Fock, Nk1, Nk2)
% FOCK_FFT_UNIFORM
%   FFT-based Fock self-energy on a uniform k-mesh.
%
%   输入:
%     rho_k    : [m x m x Nk]
%     w_k      : [Nk x 1] 均匀权重 (sum = Wnorm)
%     Wnorm    : = sum(w_k) = area_patch (valley patch) 或 area_BZ
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
    %   Σ_F(k) = -(1/Ω) ∑_{k'} w(k') V(q) ρ(k')
    %   将 w_k/Wnorm = 1/Nk (在 patch 或 uniform BZ 情况下一般是常数)
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

function SigmaH_k = hartree_uniform_smallq(g, tau, rho_k, w_k, Wnorm, opts)
% HARTREE_UNIFORM_SMALLQ
%   k-independent, diagonal Hartree shift with layer resolution.
%
%   使用 *完全相同* 的 V(q) 模型 (opts.Vq_fun)，在 q_cart = 0 的位置求
%   V_ab(q->0)，通过 q = max(|q|, qmin) 做 regularization。
%
%   公式上:
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

    % ---- 3) V_ab(q→0)：统一通过 Vq_fun(q_cart=0,opts) 计算 ----
    d = size(g.b,1);
    q_cart0 = zeros(1,d);                     % 1×d
    Vab0    = opts.Vq_fun(q_cart0, opts);     % m×m

    % ---- 4) Diagonal Hartree shift Σ_H = V(q→0) * dn ----
    SigH_diag = Vab0 * dn;                    % m×1

    % 可选：减去整体平均，避免 rigid shift（只留层/轨道间相对势能差）
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
%     vmodel = 'doublegate': same-layer / diff-layer:
%           V_s(q) = 2π e² / [q ε cosh(qD)]
%           V_d(q) = 2π e² / [q ε sinh(qD)]
%     vmodel = 'lattice'   : user-provided V_ab(R) real-space Fourier sum
%
%   这里假设 g.wpos 已经包含 valley × spin × orbital × layer，
%   层信息只通过 z = g.wpos(:,3) 体现。

    m = size(g.wpos,1);

    % 层信息：layer_z / layer_id / same-layer mask
    if ~isfield(opts,'layer_z') || isempty(opts.layer_z)
        opts.layer_z = g.wpos(:,3);
    end
    z = opts.layer_z(:);

    % 用 unique + 容差 做层归类
    if ~isfield(opts,'layer_tol') || isempty(opts.layer_tol)
        opts.layer_tol = 1e-3;   % Å 的量级，按实际 z 的 scale 调
    end
    [z_unique,~,layer_id] = uniquetol(z, opts.layer_tol, 'DataScale', 1);
    opts.layer_id  = layer_id;
    opts.z_unique  = z_unique;
    opts.same_layer_mask = (layer_id(:) == layer_id(:).');  % m×m logical

    % 常用物理参数默认值
    if ~isfield(opts,'e2') || isempty(opts.e2)
        opts.e2 = 14.3996;   % eV·Å (e^2 / 4π ε0)
    end
    if ~isfield(opts,'eps') || isempty(opts.eps)
        opts.eps = 4.0;      % 相对介电常数
    end
    if ~isfield(opts,'qmin') || isempty(opts.qmin)
        opts.qmin = 1e-2;    % small-q regulator
    end

    switch lower(opts.vmodel)
        case '2d'
            % V(q) = (2π e² / ε) * 1/q * exp(-q |z_a - z_b|)
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;
            end
            opts.Vq_fun     = @(q_cart,opts2) ...
                Vq_2d_layered_core(q_cart, opts2.alpha, opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all,opts2) ...
                Vq_2d_layered_vec_core(q_cart_all, opts2.alpha, opts2.layer_z, opts2);

        case 'keldysh'
            % V(q) = (2π e² / ε) / [q(1 + r0 q)] * exp(-q |z_a - z_b|)
            if ~isfield(opts,'alpha') || isempty(opts.alpha)
                opts.alpha = 2*pi*opts.e2/opts.eps;
            end
            if ~isfield(opts,'r0') || isempty(opts.r0)
                error('For vmodel=''keldysh'', please set opts.r0 (Keldysh length).');
            end
            opts.Vq_fun     = @(q_cart,opts2) ...
                Vq_keldysh_layered_core(q_cart, opts2.alpha, opts2.r0, ...
                                        opts2.layer_z, opts2);
            opts.Vq_vec_fun = @(q_cart_all,opts2) ...
                Vq_keldysh_layered_vec_core(q_cart_all, opts2.alpha, opts2.r0, ...
                                            opts2.layer_z, opts2);

        case 'doublegate'
            % Double-gate 模型 (简化版):
            %   same layer: V_s(q) = 2π e² / [q ε cosh(q D)]
            %   diff layer: V_d(q) = 2π e² / [q ε sinh(q D)]
            %   D ~ gate separation parameter.
            if ~isfield(opts,'Dgate') || isempty(opts.Dgate)
                error('For vmodel=''doublegate'', please set opts.Dgate (gate separation D).');
            end
            opts.Vq_fun     = @(q_cart,opts2) ...
                Vq_doublegate_layered_core(q_cart, opts2.e2, opts2.eps, ...
                                           opts2.Dgate, opts2.same_layer_mask, opts2);
            opts.Vq_vec_fun = @(q_cart_all,opts2) ...
                Vq_doublegate_layered_vec_core(q_cart_all, opts2.e2, opts2.eps, ...
                                               opts2.Dgate, opts2.same_layer_mask, opts2);

        case 'lattice'
            % 用户提供 V_ab(R) 的 Fourier sum:
            %   V_ab(q) = ∑_R V_ab(R) e^{-i q·R}
            assert(isfield(opts,'VR_R_cart') && isfield(opts,'VR_ab'), ...
                'For vmodel=''lattice'', provide opts.VR_R_cart and opts.VR_ab.');
            opts.Vq_fun     = @(q_cart,opts2) ...
                Vq_from_VR_core(q_cart, opts2.VR_R_cart, opts2.VR_ab);
            opts.Vq_vec_fun = @(q_cart_all,opts2) ...
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


% ---- Double-gate, layered (core) ---------------------------------------
function V_ab = Vq_doublegate_layered_core(q_cart, e2, epsr, Dgate, same_layer_mask, opts)
    % same-layer: Vs(q) = 2π e² / [q ε cosh(qD)]
    % diff-layer: Vd(q) = 2π e² / [q ε sinh(qD)]
    %
    % 简化假设：同一 z-layer 的任意轨道 pair 使用 Vs，其他使用 Vd。

    q = max(norm(q_cart), get_qmin(opts));
    pref = 2*pi*e2/epsr;

    Vs = pref / (q * cosh(q*Dgate));
    Vd = pref / (q * sinh(q*Dgate));

    same = same_layer_mask;
    V_ab = Vs * same + Vd * (~same);
end

function V_pages = Vq_doublegate_layered_vec_core(q_cart_all, e2, epsr, Dgate, same_layer_mask, opts)
    Nk   = size(q_cart_all,1);
    qmag = vecnorm(q_cart_all,2,2);
    qmag = max(qmag, get_qmin(opts));
    pref = 2*pi*e2/epsr;

    m    = size(same_layer_mask,1);
    same = same_layer_mask;
    V_pages = zeros(m,m,Nk);
    for j=1:Nk
        q = qmag(j);
        Vs = pref / (q * cosh(q*Dgate));
        Vd = pref / (q * sinh(q*Dgate));
        V_pages(:,:,j) = Vs*same + Vd*(~same);
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
        qmin = 1e-2;    % 默认 0.01 (可以在 main 中显式设置)
    end
end


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ wrap_frac_pm -----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xw = wrap_frac_pm(x)
% wrap each component to [-1/2, 1/2)
    xw = x - round(x);
end


function out = valley_block_line(g, k_line_frac)
% VALLEY_BLOCK_LINE
%   Build valley-block Hamiltonian along a given high-symmetry line:
%       H_valley(k) = blkdiag( H(k), H(-k) )
%
% Row-basis convention:
%   - g.a, g.b : [d x d], each ROW is a_i / b_i  (Cartesian)
%   - g.hopr   : [L x d], fractional lattice vectors ΔR
%   - g.ham    : [m x m x L], hopping t(ΔR)
%
% Inputs
%   g           : struct of Wannier TB
%   k_line_frac : [Nk x d], fractional k points along a high-symmetry line
%                 (e.g. from Γ to K, or K to M, etc.)
%                 d = 2 or 3; if 2, we automatically pad k_z = 0
%
% Output (out)
%   .k_frac     : [Nk x d]  fractional k used (after padding / wrap)
%   .k_cart     : [Nk x d]  Cartesian k = k_frac * g.b
%   .Hk         : [m x m x Nk]       H(k)
%   .H_minus_k  : [m x m x Nk]       H(-k)
%   .H_valley   : [2m x 2m x Nk]     blkdiag(H(k), H(-k))
%   .E_valley   : [Nk x 2m]          eigenvalues of H_valley(k)
%   .U_valley   : [2m x 2m x Nk]     eigenvectors of H_valley(k)

    % ---------- basic sizes ----------
    d  = size(g.b,1);            % dimension (2 or 3)
    m  = size(g.ham,1);          % orbitals per valley
    Nk = size(k_line_frac,1);

    % ---------- ensure k_line_frac has d components ----------
    if size(k_line_frac,2) < d
        k_line_frac = [k_line_frac, zeros(Nk, d-size(k_line_frac,2))];
    end

    % 如果你更习惯 [-0.5,0.5) 的 BZ，可以在这里 wrap 一下
    % k_line_frac = wrap_frac_pm(k_line_frac);

    % ---------- compute H(k) & H(-k) ----------
    Hk        = zeros(m,m,Nk);
    H_minus_k = zeros(m,m,Nk);

    for ik = 1:Nk
        kf  = k_line_frac(ik,:);           % 1×d
        kf_ = -kf;                         % -k

        Hk(:,:,ik)        = H0k_single(g, kf);
        H_minus_k(:,:,ik) = conj(H0k_single(g, kf_));
    end

    % ---------- valley-block Hamiltonian ----------
    H_valley = zeros(2*m, 2*m, Nk);
    E_valley = zeros(Nk, 2*m);
    U_valley = zeros(2*m, 2*m, Nk);

    for ik = 1:Nk
        Hv = blkdiag(Hk(:,:,ik), H_minus_k(:,:,ik));
        Hv = (Hv + Hv')/2;                          % Hermitize (数值稳定)
        H_valley(:,:,ik) = Hv;

        [U, D] = eig(Hv, 'vector');
        [e_sorted, idx] = sort(real(D), 'ascend');
        U = U(:, idx);

        E_valley(ik,:)   = e_sorted.';
        U_valley(:,:,ik) = U;
    end

    % ---------- k_cart ----------
    k_cart = k_line_frac * g.b;            % row-basis: (Nk×d)*(d×d)

    % ---------- pack output ----------
    out.k_frac    = k_line_frac;
    out.k_cart    = k_cart;
    out.Hk        = Hk;
    out.H_minus_k = H_minus_k;
    out.H_valley  = H_valley;
    out.E_valley  = E_valley;
    out.U_valley  = U_valley;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- single-k Wannier Hamiltonian --------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hk = H0k_single(g, k_frac)
% H0K_SINGLE  H(k) = sum_R t(R) e^{ i k · R }
% Row-basis convention: k_cart = k_frac * g.b,  R_cart = R_frac * g.a

    d  = size(g.b,1);
    if numel(k_frac) < d
        k_frac = [k_frac(:).', zeros(1,d-numel(k_frac))];
    else
        k_frac = k_frac(:).';
    end

    L  = size(g.hopr,1);
    m  = size(g.ham,1);

    k_cart = k_frac * g.b;        % 1×d

    Hk = zeros(m,m);
    for ell = 1:L
        dR_frac = g.hopr(ell,:);  % 1×d
        dR_cart = dR_frac * g.a;  % 1×d
        phase   = exp(1i * dot(k_cart, dR_cart));
        Hk      = Hk + g.ham(:,:,ell) * phase;
    end

    Hk = (Hk + Hk')/2;
end
