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



%%
opts.bands_active = [nv, nc];  % 相对于对角化后 1..norb 的索引
opts.nspin    = 2;
opts.kT       = 1e-3;
opts.mu       = 0.0;       % 或者设成你想要的化学势
opts.max_iter = 200;
opts.tol      = 1e-6;
opts.mix_eta  = 0.5;
opts.verbose  = true;

out = hfmf_band_valley(g, kmesh, Vq, opts);

%%
function out = hfmf_band_valley(g, kmesh, Vq, opts)
% HFMF_BAND_VALLEY  Band-basis Hartree-Fock for K/K' valleys + spin.
%
% 这是一个"band-projected HFMF"的主程序：
%   1) 从 Wannier TB (g) 构造 H0(k)，并在每个 k 上对角化得到 U(k), eps(k)；
%   2) 只保留费米能附近的一小撮 active bands，在 band basis 上做 HF 自洽；
%   3) 相互作用 (Hartree-Fock) 仍在轨道基底上计算，再投影回 band basis。
%
% 输入:
%   g     : 结构体，Wannier TB 信息（与你原来的 HFMF 一致），至少需要:
%           g.hopr  [nR x dim]   跳跃向量（分数坐标, lattice gauge）
%           g.ham   [norb x norb x nR] 对应每个 hopr 的 hopping 矩阵
%           （可视为 spinless TB；spin 用 HF 自洽来产生）
%
%   kmesh : 结构体，k 网格和 valley patch 信息，例如:
%           kmesh.kpts   [Nk x dim]   k 点 (分数坐标, 比如 (kx,ky))
%           kmesh.weight [Nk x 1]     每个 k 的权重（否则可全 1/Nk）
%           kmesh.valley [Nk x 1]     valley 标签: +1 对应 K, -1 对应 K'
%
%   Vq    : 结构体，相互作用信息。这里为了简单，只用到 q=0 的 V_ab(0):
%           Vq.V0   [norb x norb]    V_ab(q=0)，即实空间 V_ab(R) 的总和
%           （你完全可以扩展成一般 q，用你的 FFT 代码做完整 HF）
%
%   opts  : 结构体，控制参数，例如:
%       opts.bands_active : [Nb_act x 1]  要保留的 active band 索引 (相对于对角化后全部带)
%       opts.nspin        : 自旋数目 (一般 2)
%       opts.kT           : 温度 (能量单位)
%       opts.mu           : 化学势（简单起见这里固定，不做填充自洽）
%       opts.max_iter     : HF 最大迭代数
%       opts.tol          : 收敛阈值 (max|rho_new - rho_old|)
%       opts.mix_eta      : 密度矩阵 mixing 参数 (0<eta<=1)
%       opts.verbose      : 是否输出迭代信息
%
% 输出:
%   out : 结构体，包含:
%       out.eps_band   [Nb_act x Nk]         非相互作用能谱（active bands）
%       out.U_act      [norb x Nb_act x Nk]  对应 U(k) 的列子矩阵
%       out.rho_band   [Nb_act x Nb_act x Nk x nspin]  自洽后的密度矩阵
%       out.Sigma_band [Nb_act x Nb_act x Nk x nspin]  自洽后的 HF 自能
%       out.E_HF       [Nb_act x Nk x nspin] HF 本征能谱
%       out.iters      实际迭代步数
%       out.converged  是否达到收敛
%
% 注意:
%   - 这里的 HF builder orbital_HF_builder 只实现了 Hartree 项，
%     Fock 项留有接口，你可以直接替换成你已有的 FFT + Fock 的轨道基代码；
%   - 这里 valley 只通过 k 属于 K/K' patch 间接体现，真正的 valley 相关
%     序参量取决于你如何在初始 rho 中设置 seed (比如给 K / K' 不同占据)。

%% ============= 0. 预设缺省参数 ===========================
if ~isfield(opts, 'nspin'),      opts.nspin = 2;          end
if ~isfield(opts, 'kT'),        opts.kT     = 1e-3;       end
if ~isfield(opts, 'mu'),        opts.mu     = 0.0;        end
if ~isfield(opts, 'max_iter'),  opts.max_iter = 200;      end
if ~isfield(opts, 'tol'),       opts.tol    = 1e-6;       end
if ~isfield(opts, 'mix_eta'),   opts.mix_eta = 0.5;       end
if ~isfield(opts, 'verbose'),   opts.verbose = true;      end

nspin = opts.nspin;

%% ============= 1. 从 Wannier TB 得到 band basis (spinless) =========
% 目标:
%   对每个 k, 构造 H0(k) 并对角化:
%     h(k) U(k) = U(k) diag(eps(k))
%   然后只保留 active bands (near-EF)。

kpts   = kmesh.kpts;      % [Nk x dim]
Nk     = size(kpts, 1);
hopr   = g.hopr;          % [nR x dim]
hamR   = g.ham;           % [norb x norb x nR]
norb   = size(hamR, 1);
ib_act = opts.bands_active(:);  % active band indices (相对于 full band)

Nb_act = numel(ib_act);

eps_band = zeros(Nb_act, Nk);         % spinless 冷带能量 epsilon_n(k)
U_act    = zeros(norb, Nb_act, Nk);   % 对应的 U(k) 子矩阵

% === 并行 对每个 k 计算 H0(k) 并对角化 ===
% 注意: parfor 需要 Parallel Toolbox.
parfor ik = 1:Nk
    kvec = kpts(ik, :);                      % 当前 k (分数坐标)
    Hk   = build_hk_from_hopr(hopr, hamR, kvec);  % [norb x norb] 轨道基 H(k)
    [Uk, Dk] = eig(Hk, 'vector');           % Uk 列为本征态, Dk 为本征值
    % 选择 active bands
    eps_band(:, ik) = Dk(ib_act);
    U_act(:, :, ik) = Uk(:, ib_act);
end

%% ============= 2. 初始化 band-space 密度矩阵 rho_band = <d^† d> =======
% rho_band 维度: [Nb_act x Nb_act x Nk x nspin]
%   - Nb_act: active bands (例如 conduction & valence)
%   - Nk    : k 网格
%   - nspin : spin ↑/↓ 通道
%
% 默认: 从非相互作用的 Fermi 分布开始 (对角, 无自旋/valley 极化)
rho_band = zeros(Nb_act, Nb_act, Nk, nspin);

kT = opts.kT;
mu = opts.mu;

for is = 1:nspin
    for ik = 1:Nk
        eps_vec = eps_band(:, ik);      % [Nb_act x 1]
        occ     = fermi_function(eps_vec - mu, kT); % [Nb_act x 1]
        rho_band(:, :, ik, is) = diag(occ);        % 初始: 只对角占据
    end
end

% 你可以在这里根据 valley (kmesh.valley) + spin is 来加入 seed，
% 实现自旋极化 / valley 极化 / LAF seeds：
% 例如:
%   valley_sign = kmesh.valley(ik);   % +1 for K, -1 for K'
%   if valley_sign>0 && is==1 (K, spin-up) -> 加一点额外占据等
% 我这里先不给 seed，保持对称，让 HF 自己决定是否自发破缺。

%% ============= 3. HF 自洽循环 =======================================
Sigma_band = zeros(Nb_act, Nb_act, Nk, nspin);  % HF 自能 (band-space)
E_HF       = zeros(Nb_act, Nk, nspin);          % HF 本征能

converged = false;

for it = 1:opts.max_iter
    rho_old = rho_band;
    
    % 3.1 band-space rho -> orbital-space rho
    % rho_orb 维度: [norb x norb x Nk x nspin]
    rho_orb = band2orbital_density(rho_band, U_act);
    
    % 3.2 在轨道基底上构建 HF 自能 (目前实现: Hartree; Fock 留接口)
    % Sigma_orb: [norb x norb x Nk x nspin]
    Sigma_orb = orbital_HF_builder(rho_orb, Vq, kmesh);
    
    % 3.3 将轨道自能投影回 band basis
    % Sigma_band: [Nb_act x Nb_act x Nk x nspin]
    Sigma_band = orbital2band_selfenergy(Sigma_orb, U_act);
    
    % 3.4 构建 band-space HF Hamiltonian 并对角化 -> 更新 rho_band
    rho_new = zeros(size(rho_band));
    parfor ik = 1:Nk   % k 点之间相互独立，可以并行
        for is = 1:nspin
            H_HF = diag(eps_band(:, ik)) + Sigma_band(:, :, ik, is);  % [Nb_act x Nb_act]
            [W, Ehf] = eig(H_HF, 'vector');    % W columns = eigenvectors, Ehf = eigenvalues
            E_HF(:, ik, is) = Ehf;             % 记录 HF 本征能量
            
            occ = fermi_function(Ehf - mu, kT);  % [Nb_act x 1]
            % HF 态下的密度矩阵: rho = W f(E) W^†
            rho_new(:, :, ik, is) = W * (occ .* W');   % 等价于 W*diag(occ)*W'
        end
    end
    
    % 3.5 mixing + 收敛判断
    eta = opts.mix_eta;
    rho_band = (1 - eta) * rho_band + eta * rho_new;
    
    delta = max(abs(rho_band(:) - rho_old(:)));
    
    if opts.verbose
        fprintf('HF iter %4d: max|Δrho| = %.3e\n', it, delta);
    end
    
    if delta < opts.tol
        converged = true;
        break;
    end
end

%% ============= 4. 打包输出 ==========================================
out.eps_band   = eps_band;
out.U_act      = U_act;
out.rho_band   = rho_band;
out.Sigma_band = Sigma_band;
out.E_HF       = E_HF;
out.iters      = it;
out.converged  = converged;

end % ====== 主函数结束 ================================================


%% ====================================================================
%% 子函数 1: 从 hopr + hamR 构造 H(k) (轨道基)
%% ====================================================================
function Hk = build_hk_from_hopr(hopr, hamR, kvec)
% BUILD_HK_FROM_HOPR  构造给定 kvec 下的 H(k)
%
% 输入:
%   hopr  [nR x dim]      跳跃向量 ΔR (分数坐标, lattice gauge)
%   hamR  [norb x norb x nR] 对应每个 ΔR 的 hopping 矩阵
%   kvec  [1 x dim]       当前 k 点 (分数坐标, 与 hopr 的坐标系一致)
%
% 输出:
%   Hk    [norb x norb]   轨道基 H(k) 矩阵

[norb, ~, nR] = size(hamR);
Hk = zeros(norb, norb);
for ir = 1:nR
    dR = hopr(ir, :);                 % ΔR (分数坐标)
    phase = exp(-1i * 2*pi * (dR * kvec.'));  % lattice gauge 相位
    Hk = Hk + hamR(:, :, ir) * phase;
end
% 确保厄米：
Hk = (Hk + Hk')/2;
end


%% ====================================================================
%% 子函数 2: band-space rho -> orbital-space rho
%% ====================================================================
function rho_orb = band2orbital_density(rho_band, U_act)
% BAND2ORBITAL_DENSITY  将 band-basis 密度矩阵投影到 Wannier/orbital 基底
%
% 输入:
%   rho_band  [Nb_act x Nb_act x Nk x nspin]
%   U_act     [norb x Nb_act x Nk]   轨道 -> active band 的变换矩阵
%
% 输出:
%   rho_orb   [norb x norb x Nk x nspin]
%
% 公式:
%   rho_orb_ab(k,s) = Σ_{n,m in active} U_{a n}(k) rho_{n m}(k,s) U^*_{b m}(k)

[norb, Nb_act, Nk] = size(U_act);
[~, ~, ~, nspin]   = size(rho_band);

rho_orb = zeros(norb, norb, Nk, nspin);

for is = 1:nspin
    for ik = 1:Nk
        Uk  = U_act(:, :, ik);                 % norb x Nb_act
        rb  = rho_band(:, :, ik, is);          % Nb_act x Nb_act
        rho_orb(:, :, ik, is) = Uk * rb * Uk'; % norb x norb
    end
end
end


%% ====================================================================
%% 子函数 3: 轨道基 Hartree-Fock builder (此处实现 Hartree, Fock 留接口)
%% ====================================================================
function Sigma_orb = orbital_HF_builder(rho_orb, Vq, kmesh)
% ORBITAL_HF_BUILDER  在 Wannier/orbital 基底上构造 HF 自能
%
% 当前版本:
%   只实现 Hartree 自能: Σ_H_ab(k,s) = δ_ab * Σ_c V_ac(0) * n_c
%   其中:
%       n_c = (1/Nk) Σ_{k,s} rho_orb_cc(k,s) (总的轨道密度)
%
%   Fock 项 (k 依赖, 非对角) 留给你用已有的 FFT + 卷积代码来实现:
%       Σ_F_ab(k,s) = - (1/Nk) Σ_{k',s'} V_ab(k-k') * ρ_ba(k',s')
%   你可以在这里直接替换为你之前写好的 Python/MATLAB 实现。
%
% 输入:
%   rho_orb  [norb x norb x Nk x nspin]
%   Vq       结构体, 至少包含:
%       Vq.V0  [norb x norb]  V_ab(q=0) (Hartree 用)
%   kmesh    结构体, 包含:
%       kmesh.weight [Nk x 1] k 点权重 (若无可统一设为 1/Nk)
%
% 输出:
%   Sigma_orb [norb x norb x Nk x nspin]

[norb, ~, Nk, nspin] = size(rho_orb);

Sigma_orb = zeros(norb, norb, Nk, nspin);

% ---- Hartree: 先计算每个轨道的总密度 n_a ----
if isfield(kmesh, 'weight')
    wk = kmesh.weight(:);
else
    wk = ones(Nk, 1) / Nk;
end

n_orb = zeros(norb, 1);  % 每个轨道的总粒子数 (包括自旋,k加权)
for is = 1:nspin
    for ik = 1:Nk
        rho_k = rho_orb(:, :, ik, is);    % norb x norb
        n_orb = n_orb + wk(ik) * diag(real(rho_k));   % 只取对角 + 实部
    end
end
% 如果需要固定总填充，可以在这里重新规范化 n_orb.

% Hartree potential: V_H(a) = Σ_b V_ab(0) * n_b
V0 = Vq.V0;   % [norb x norb]
V_H = V0 * n_orb;   % [norb x 1]

% 将 Hartree 视作对 orbital onsite 的修正: Σ_H_ab = δ_ab * V_H(a)
for is = 1:nspin
    for ik = 1:Nk
        Sigma_orb(:, :, ik, is) = diag(V_H);
    end
end

% ---- Fock: 留接口 (示意), 你可以完全替换这一段 ----
%{
% 示例 (极简, 不实际): 使用某个近似的 V_Fock_ab(k) ...
for is = 1:nspin
    for ik = 1:Nk
        rho_k = rho_orb(:, :, ik, is);
        % 这里调用你自己的 Fock builder, 例如:
        % Sigma_Fock = build_Fock_qspace(rho_orb, Vq, ik, is, kmesh);
        % 然后:
        % Sigma_orb(:, :, ik, is) = Sigma_orb(:, :, ik, is) + Sigma_Fock;
    end
end
%}

end


%% ====================================================================
%% 子函数 4: orbital-space Σ -> band-space Σ
%% ====================================================================
function Sigma_band = orbital2band_selfenergy(Sigma_orb, U_act)
% ORBITAL2BAND_SELFENERGY  将轨道基自能投影到 band basis
%
% 输入:
%   Sigma_orb [norb x norb x Nk x nspin] 自能(轨道基)
%   U_act     [norb x Nb_act x Nk]       轨道 -> active band 变换
%
% 输出:
%   Sigma_band [Nb_act x Nb_act x Nk x nspin] 自能(带基)
%
% 公式:
%   Σ_band_nm(k,s) = Σ_{a,b} U^*_{a n}(k) Σ_orb_ab(k,s) U_{b m}(k)

[norb, ~, Nk, nspin] = size(Sigma_orb);
[~, Nb_act, ~]       = size(U_act);

Sigma_band = zeros(Nb_act, Nb_act, Nk, nspin);

for is = 1:nspin
    for ik = 1:Nk
        Uk   = U_act(:, :, ik);           % [norb x Nb_act]
        Sigma_k = Sigma_orb(:, :, ik, is);% [norb x norb]
        % U^† Σ U:
        Sigma_band(:, :, ik, is) = Uk' * Sigma_k * Uk;   % [Nb_act x Nb_act]
        % 确保厄米:
        Sb = Sigma_band(:, :, ik, is);
        Sigma_band(:, :, ik, is) = (Sb + Sb')/2;
    end
end

end


%% ====================================================================
%% 辅助: Fermi 函数
%% ====================================================================
function f = fermi_function(E, kT)
% FERMI_FUNCTION  Fermi-Dirac 分布
%   f(E) = 1 / (exp(E/kT) + 1)
% 输入:
%   E  : 能量 (可以是向量)
%   kT : 温度 (能量单位)
if kT <= 0
    f = double(E <= 0);
else
    f = 1 ./ (exp(E./kT) + 1);
end
end
