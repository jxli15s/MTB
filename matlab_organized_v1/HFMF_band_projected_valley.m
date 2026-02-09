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
% % kpts=zeros(100,2);
% % kpts(:,1)=linspace(0,1/3,100);
% % kpts(:,2)=linspace(0,2/3,100);
% % kmesh.kpts=kpts;
% % bands_active=(17:24).';
% % H_seed_orb=[];
% % [eps_band, U_act] = compute_valley_doubled_bands(g, kmesh, bands_active, H_seed_orb);
% % %%
% % figure()
% % for i =1:8
% % hold on;
% % plot(eps_band(i,:))
% % end
%%
% 假设你已经有 g 结构 (g.a, g.b, g.hopr, g.ham, g.wpos...)

Nkx = 51;
Nky = 51;

% 默认中心 K=(1/3,2/3), delta=0.15
kmesh = build_kmesh_K_patch(g, Nkx, Nky);

% ==== 从 g.wpos 构造 single-valley double-gate V ====
params.d_gate      = 300;        % Å, 你设定的 gate 距离
params.e2_over_eps = 14.4/40;    % eV·Å
V_orb0 = build_V_doublegate_from_wpos_singlevalley(g, kmesh, params);

% ==== 提升到 double-valley orbital basis ====
V = lift_V_singlevalley_to_valley(V_orb0);
%
% ==== HF 选取活跃带 ====
opts.bands_active = (17:24).';  % 你自己挑 Fermi 附近的 Nb_act 个 band
opts.kT     = 8.617333262e-5*0.2; %or 0
opts.kT = 0;
% opts.mu     =  0;
opts.max_iter = 200;
opts.tol    = 1e-6;
opts.mix_eta= 0.5;
opts.verbose= true;

% valley seed 例子：
opts.seed.type = 'valley';
opts.seed.amp  = 10e-3;      % ~ 1 meV 的 valley Zeeman seed

out = hfmf_band_valley(g, kmesh, V, opts);
%%

%%

%%
E=reshape(out.E_HF,8,Nkx,Nky);
% E=reshape(eps_seed_band,8,Nkx,Nky);
%%
figure()
for j=1:51
for i =1:8
hold on;
plot(E(i,:,j))
end
end
%%
figure()
for i =1:8
hold on;
surf(squeeze(E(i,:,:)),'EdgeColor','none')
end

%%
function out = hfmf_band_valley(g, kmesh, V, opts)
% HFMF_BAND_VALLEY
% --------------------------------------------------------------
% Band-basis Hartree-Fock starting from *single-valley* lattice
% Wannier TB, with explicit valley doubling:
%
%   H_lat(k) = diag( H_single(k), conj(H_single(-k)) ).
%
% 轨道基自由度 a = 1..norb 已经包含 valley×spin×layer×sublattice 等，
% 其中 norb = 2 * norb0，对应 (valley K 的 norb0 轨道, valley K' 的 norb0 轨道)。
%
% 重要：seed 只用于生成 *初始密度矩阵* ρ_seed，HF 自洽过程中不再包含 seed 外场。
%
% ===================== INPUT =====================
% g  (single-valley Wannier TB):
%   g.hopr   [nR x dim]            hopping 向量 ΔR (fractional)
%   g.ham    [norb0 x norb0 x nR]  单 valley 跳跃矩阵 t(ΔR)
%   g.wpos   [norb0 x 3]           Wannier 轨道位置 (Å)，这里只用 z 分量
%   g.a      [3 x 3]               实空间基矢 (Å)
%   g.b      [3 x 3]               倒空间基矢 (1/Å)
%
% kmesh:
%   kmesh.kpts   [Nk x dim]        fractional k 点，扁平化 Nk = Nkx*Nky
%   kmesh.Nkx    标量              kx 方向 k 点数
%   kmesh.Nky    标量              ky 方向 k 点数
%   (可选) kmesh.weight [Nk x 1]   物理 k 权重 w_k = dkx dky / (2π)^2
%
% V (double-valley orbital basis 相互作用):
%   V.V0   [norb x norb]                 Hartree 核 V^0_ab (q→0)
%   V.V_r  [Nkx x Nky x norb x norb]     Fock 用 V_ab(r)
%   V.A    标量                          归一化 A = Nk * Acell
%
% opts:
%   opts.bands_active [Nb_act x 1]  选出来的 active bands（按能量排序后 index）
%   opts.kT           温度 k_B T (eV)
%   opts.mu           化学势 μ (eV) —— 当前版本 μ 不自洽
%   opts.max_iter     最大 HF 迭代步数
%   opts.tol          收敛阈值 max|Δρ|
%   opts.mix_eta      mixing 参数 0<eta≤1
%   opts.verbose      是否打印迭代信息
%   opts.kT           温度 k_B T (eV)
%   opts.mu           初始化学势 μ 初值 (eV)
%   opts.N_target     若存在且 opts.fix_mu ~= true，则用它作为目标粒子数：
%                     N(μ) = Σ_k w_k Σ_n f(E_HF - μ) = N_target
%                     否则默认 N_target = (Nb_act/2)*Σ_k w_k （半填充）
%
%   (可选) opts.seed:
%       .type = 'none' | 'valley' | 'custom_orbital'
%       .amp  = 小的实数 (eV 量级)，对 'valley' 有效
%       若 type='custom_orbital'，还需:
%       .H_orbital [norb x norb]  你自定义的 seed 哈密顿量（valley⊗orbital）
%
% ===================== OUTPUT =====================
% out.eps_band   [Nb_act x Nk]           无 seed 的 H0 带能 ε_n(k)
% out.U_act      [norb x Nb_act x Nk]    无 seed 的 Bloch 波函数
% out.rho_band   [Nb_act x Nb_act x Nk]  自洽后的 band 密度矩阵
% out.Sigma_band [Nb_act x Nb_act x Nk]  HF 自能 Σ_HF(k)
% out.E_HF       [Nb_act x Nk]           HF 准粒子能量
% out.E_tot      标量                    总能量/胞 (E0+E_int)
% out.E0         标量                    非相互作用部分
% out.E_int      标量                    相互作用部分
% out.iters      实际迭代步数
% out.converged  logical, 是否收敛
% --------------------------------------------------------------

%% 0. 默认参数
if ~isfield(opts, 'kT'),       opts.kT = 1e-3;   end
if ~isfield(opts, 'mu'),       opts.mu = 0.0;    end
if ~isfield(opts, 'max_iter'), opts.max_iter = 200; end
if ~isfield(opts, 'tol'),      opts.tol = 1e-6;  end
if ~isfield(opts, 'mix_eta'),  opts.mix_eta = 0.5; end
if ~isfield(opts, 'verbose'),  opts.verbose = true; end
if ~isfield(opts, 'seed'),     opts.seed.type = 'none'; opts.seed.amp = 0.0; end

%% 1. 先对 *无 seed* 的 H0 做 double-valley 带结构
[eps0_band, U0_act] = compute_valley_doubled_bands(g, kmesh, opts.bands_active, []);
[Nb_act, Nk] = size(eps0_band);
norb = size(U0_act, 1);

kT = opts.kT;
mu = opts.mu;

%% 2. 用 seed 生成初始密度矩阵 rho_band_init
seed = opts.seed;
use_seed = isfield(seed, 'type') && ~strcmpi(seed.type, 'none');

if ~use_seed || (isfield(seed, 'amp') && seed.amp == 0 && ~strcmpi(seed.type, 'custom_orbital'))
    % ---- 无 seed 情况：直接用 H0 的 eps0_band 做 Fermi 填充 ----
    eps_band = eps0_band;
    U_act    = U0_act;

    rho_band = zeros(Nb_act, Nb_act, Nk);
    parfor ik = 1:Nk
        eps_vec = eps_band(:, ik);
        occ     = fermi_function(eps_vec - mu, kT);
        rho_band(:,:,ik) = diag(occ);
    end

else
    % ---- 有 seed：seed 只用来生成初始 rho，不参与后续 HF 外场 ----

    % 2.1 在 valley⊗orbital 轨道基上构造 H_seed_orb [norb x norb]
    H_seed_orb = build_seed_orbital(g, seed);   % Hermitian

    % 2.2 用 H0 + H_seed_orb 再对角化一次，得到 seed 带结构
    [eps_seed_band, U_seed_act] = compute_valley_doubled_bands( ...
        g, kmesh, opts.bands_active, H_seed_orb);

    % 2.3 在 seed 的 band 基底下做 Fermi 填充 → ρ_seed_band
    % 先根据 eps_seed_band 自洽求一个与之对应的 mu_seed
    if isfield(opts, 'fix_mu') && opts.fix_mu
        % 固定化学势模式：seed 也用同一个 μ
        mu_seed = opts.mu;
    else
        % 自洽模式：用 seed 光谱 eps_seed_band 求 μ_seed
        mu_seed = solve_mu_band(eps_seed_band, kT, kmesh, Nb_act, opts);
    end

    rho_seed_band = zeros(Nb_act, Nb_act, Nk);
    parfor ik = 1:Nk
        eps_vec = eps_seed_band(:, ik);
        occ     = fermi_function(eps_vec - mu_seed, kT);
        rho_seed_band(:,:,ik) = diag(occ);
    end

    % 可选：把 HF 循环的初始 μ 也设成 mu_seed，作为一个合理起点
    mu = mu_seed;


    % 2.4 把 ρ_seed_band 转成 轨道基密度矩阵 ρ_orb_seed
    rho_orb_seed = band2orbital_density(rho_seed_band, U_seed_act);  % [norb x norb x Nk]

    % 2.5 再用 *无 seed* 的波函数 U0_act 把 ρ_orb_seed 投到 H0 的 band 基底
    rho_band = orbital2band_density(rho_orb_seed, U0_act);

    % 2.6 后续 HF 一律使用 H0 的本征系作为基底
    eps_band = eps0_band;
    U_act    = U0_act;
end

%% 3. 自洽 HF 循环（不再包含 seed 外场）
Sigma_band = zeros(Nb_act, Nb_act, Nk);
E_HF       = zeros(Nb_act, Nk);
converged  = false;

for it = 1:opts.max_iter
    rho_old = rho_band;

    % 3.1 band → orbital
    rho_orb = band2orbital_density(rho_band, U_act);  % [norb x norb x Nk]

    % 3.2 轨道基：Hartree + Fock（double-gate + valley 结构都内嵌在 V 里）
    Sigma_orb = orbital_HF_fft(rho_orb, V, kmesh);    % [norb x norb x Nk]

    % 3.3 orbital → band
    Sigma_band = orbital2band_selfenergy(Sigma_orb, U_act);

    % 3.4 在 H0 band 基底中对角化 H_HF(k) = diag(eps_band) + Σ_HF(k)
    %     先对所有 k 求 Ehf 和 W_all，之后统一用 Ehf 求新的 μ，再算 rho_new。
    E_HF   = zeros(Nb_act, Nk);
    W_all  = zeros(Nb_act, Nb_act, Nk);

    for ik = 1:Nk
        H0k  = diag(eps_band(:, ik));     % 无 seed 的 H0 (band basis)
        Sigk = Sigma_band(:, :, ik);      % HF 自能
        H_HF = H0k + Sigk;                % 有效 HF Hamiltonian

        [W, D_HF] = eig(H_HF);
        Ehf = real(diag(D_HF));

        % 可选：重新排序，让能量从小到大排列
        [Ehf_sorted, idx] = sort(Ehf, 'ascend');
        W = W(:, idx);

        E_HF(:, ik)      = Ehf_sorted;
        W_all(:, :, ik)  = W;
    end

    % 3.5 根据当前 HF 能谱更新化学势 μ：

    mu = solve_mu_band(E_HF, kT, kmesh, Nb_act, opts);


    % 3.6 用新的 μ 和 Ehf 构造 rho_new
    rho_new = zeros(size(rho_band));
    for ik = 1:Nk
        W   = W_all(:, :, ik);
        occ = fermi_function(E_HF(:, ik) - mu, kT);
        rho_new(:,:,ik) = W * diag(occ) * W';
    end

    % 3.7 mixing + 收敛判断
    eta = opts.mix_eta;
    rho_band = (1 - eta) * rho_band + eta * rho_new;

    delta = max(abs(rho_band(:) - rho_old(:)));
    if opts.verbose
        fprintf('HF iter %4d: max|Δrho| = %.3e, mu = %.6f eV\n', it, delta, mu);
    end

    if delta < opts.tol
        converged = true;
        break;
    end
end


%% 4. 用 k 权重计算总能量（H0 不含 seed）
[E_tot, E0, E_int] = compute_energy_band(eps_band, Sigma_band, rho_band, kmesh);

%% 5. 打包输出
out.eps_band   = eps_band;
out.U_act      = U_act;
out.rho_band   = rho_band;
out.Sigma_band = Sigma_band;
out.E_HF       = E_HF;
out.E_tot      = E_tot;
out.E0         = E0;
out.E_int      = E_int;
out.iters      = it;
out.converged  = converged;

end

function [eps_band, U_act] = compute_valley_doubled_bands(g, kmesh, bands_active, H_seed_orb)
% COMPUTE_VALLEY_DOUBLED_BANDS
% --------------------------------------------------------------
% 对 valley-doubled Hamiltonian 在整张 k 网格上对角化：
%
%   H_full(k) = diag( H_single(k), conj(H_single(-k)) ) + H_seed_orb
%
% H_seed_orb 是 (valley⊗orbital) 轨道基上的 k 无关 seed，
% 若为空 [] 或没给，则认为无 seed。
%
% 输入:
%   g.hopr      [nR x dim]
%   g.ham       [norb0 x norb0 x nR]
%   kmesh.kpts  [Nk x dim]
%   bands_active [Nb_act x 1]
%   H_seed_orb  [2*norb0 x 2*norb0] or []  (可选)
%
% 输出:
%   eps_band [Nb_act x Nk]
%   U_act    [2*norb0 x Nb_act x Nk]

kpts  = kmesh.kpts;
kpts  = kpts + [1/3,2/3];
Nk    = size(kpts, 1);
hopr  = g.hopr;
hamR0 = g.ham;                      % 单 valley
norb0 = size(hamR0, 1);
norb  = 2 * norb0;

ib_act = bands_active(:);
Nb_act = numel(ib_act);

if nargin < 4 || isempty(H_seed_orb)
    H_seed_orb = zeros(norb);
else
    if ~isequal(size(H_seed_orb), [norb, norb])
        error('H_seed_orb size mismatch: expected [%d x %d].', norb, norb);
    end
    H_seed_orb = 0.5 * (H_seed_orb + H_seed_orb');  % Hermitian
end

eps_band = zeros(Nb_act, Nk);
U_act    = zeros(norb, Nb_act, Nk);

parfor ik = 1:Nk
    kvec = kpts(ik, :);                        % fractional k
    Hfull = build_hk_from_hopr(hopr, hamR0, kvec);
    Hfull = Hfull + H_seed_orb;

    [Uk, Dk] = eig(Hfull);
    eps_all  = real(diag(Dk));

    [eps_sorted, idx_sort] = sort(eps_all, 'ascend');
    Uk_sorted = Uk(:, idx_sort);

    eps_band(:, ik) = eps_sorted(ib_act);
    U_act(:, :, ik) = Uk_sorted(:, ib_act);
end

end


function Hfull = build_hk_from_hopr(hopr, hamR, kvec)
% BUILD_HK_FROM_HOPR
%   H_full(k) = diag( H(k), conj(H(-k)) )

Hk   = build_hk_singlevalley(hopr, hamR,  kvec); % H(k)
Hkm  = build_hk_singlevalley(hopr, hamR, -kvec); % H(-k)
Hfull = blkdiag(Hk, Hkm);

end


function Hk = build_hk_singlevalley(hopr, hamR, kvec)
% BUILD_HK_SINGLEVALLEY
%   H(k) = Σ_R t(ΔR) e^{i 2π k·ΔR}  (lattice gauge)

[norb0, ~, nR] = size(hamR);
Hk = zeros(norb0, norb0);
for ir = 1:nR
    dR    = hopr(ir, 1:2);
    phase = exp(1i * 2*pi * (dR * kvec.'));
    Hk    = Hk + hamR(:,:,ir) * phase;
end
Hk = 0.5 * (Hk + Hk');  % enforce Hermitian

end


function rho_orb = band2orbital_density(rho_band, U_act)
% BAND2ORBITAL_DENSITY
%   ρ_orb(k) = U_act(k) ρ_band(k) U_act(k)^\dagger

[norb, Nb_act, Nk] = size(U_act);
rho_orb = zeros(norb, norb, Nk);

parfor ik = 1:Nk
    Uk = U_act(:,:,ik);       % [norb x Nb_act]
    rb = rho_band(:,:,ik);    % [Nb_act x Nb_act]
    rho_orb(:,:,ik) = Uk * rb * Uk';
end

end


function rho_band = orbital2band_density(rho_orb, U_act)
% ORBITAL2BAND_DENSITY
%   ρ_band(k) = U_act(k)^\dagger ρ_orb(k) U_act(k)

[norb, ~, Nk] = size(rho_orb);
[~, Nb_act, ~] = size(U_act);

rho_band = zeros(Nb_act, Nb_act, Nk);

parfor ik = 1:Nk
    Uk  = U_act(:,:,ik);        % [norb x Nb_act]
    rho = rho_orb(:,:,ik);      % [norb x norb]
    rho_band(:,:,ik) = Uk' * rho * Uk;
end

end


function Sigma_band = orbital2band_selfenergy(Sigma_orb, U_act)
% ORBITAL2BAND_SELFENERGY
%   Σ_band(k) = U_act(k)^\dagger Σ_orb(k) U_act(k)

[norb, ~, Nk] = size(Sigma_orb);
[~, Nb_act, ~] = size(U_act);

Sigma_band = zeros(Nb_act, Nb_act, Nk);

parfor ik = 1:Nk
    Uk    = U_act(:,:,ik);
    Sigma = Sigma_orb(:,:,ik);
    Sb    = Uk' * Sigma * Uk;
    Sigma_band(:,:,ik) = 0.5 * (Sb + Sb');  % Hermitian
end

end

function Sigma_orb = orbital_HF_fft(rho_orb, V, kmesh)
% ORBITAL_HF_FFT
% --------------------------------------------------------------
% 轨道基 Hartree + Fock 自能：
%
%   Hartree: Σ_H_ab(k) = δ_ab Σ_c V0_ac * n_c
%            其中 V0_ac = V.V0(a,c) 来自 q->0 的 double-gate V(q;z,z')
%
%   Fock   : Σ_F_ab(k) = -Σ_q V_ab(q) ρ_ab(k - q)
%
% Fock 用 FFT 做卷积：
%   ρ(k) → FFT → ρ(r)
%   V(q)→IFFT→V(r)
%   Σ_F(r) = - V(r)/A * ρ(r)
%   再 FFT 回 k。
%
% 注意：
%   - 这里完全不再手写 VS0/VD0/VI，一切 q→0 信息都已经编码在 V.V0。
%   - V.V_r 由 double-gate + g.wpos 构造（single valley），再 lift 到 double valley。
%
% 输入:
%   rho_orb [norb x norb x Nk]     ρ_ab(k)
%   V.V0    [norb x norb]          Hartree 核
%   V.V_r   [Nkx x Nky x norb x norb]  Fock 核 (r-space)
%   V.A     标量                    归一化 A = Nk * Acell
%   kmesh.Nkx, kmesh.Nky
%   (可选) kmesh.weight [Nk x 1]    k 权重
%
% 输出:
%   Sigma_orb [norb x norb x Nk] = Σ_H + Σ_F

[norb, ~, Nk] = size(rho_orb);
Nkx = kmesh.Nkx;
Nky = kmesh.Nky;

if Nk ~= Nkx * Nky
    error('orbital_HF_fft: Nk != Nkx*Nky, check kmesh flattening.');
end

%% 1. Hartree: Σ_H_ab(k) = δ_ab Σ_c V0_ac n_c
if isfield(kmesh, 'weight')
    wk = kmesh.weight(:);
else
    wk = ones(Nk,1) / Nk;
end

n_orb = zeros(norb,1);
for ik = 1:Nk
    rho_k = rho_orb(:,:,ik);
    n_orb = n_orb + wk(ik) * real(diag(rho_k));
end

if isfield(V, 'V0') && ~isempty(V.V0)
    V0 = V.V0;
else
    V0 = zeros(norb);
end

V_H_vec = V0 * n_orb;       % Σ_c V0_ac n_c
Sigma_H = zeros(norb, norb, Nk);
for ik = 1:Nk
    Sigma_H(:,:,ik) = diag(V_H_vec);
end

%% 2. Fock: Σ_F_ab(k) via FFT
rho_k = zeros(Nkx, Nky, norb, norb);
for ik = 1:Nk
    [ix, iy] = ind2sub([Nkx, Nky], ik);
    rho_k(ix,iy,:,:) = rho_orb(:,:,ik);
end

rho_r = fft2(rho_k);    % [Nkx x Nky x norb x norb]

V_r_raw = V.V_r;
A       = V.A;

% 允许两种形式:
%   (1) [Nkx x Nky]：对所有 (a,b) 相同的 V(r)
%   (2) [Nkx x Nky x norb x norb]：轨道依赖 V_ab(r)
if ndims(V_r_raw) == 2
    V_r = repmat(V_r_raw / A, 1,1,norb,norb);
elseif ndims(V_r_raw) == 4
    [Nkx2, Nky2, n1, n2] = size(V_r_raw);
    if Nkx2 ~= Nkx || Nky2 ~= Nky || n1 ~= norb || n2 ~= norb
        error('orbital_HF_fft: size(V.V_r) incompatible with rho_orb.');
    end
    V_r = V_r_raw / A;
else
    error('orbital_HF_fft: V.V_r must be 2D or 4D.');
end

Sigma_r   = - V_r .* rho_r;      % Σ_F(r) = -V(r)/A * ρ(r)
Sigma_F_k = ifft2(Sigma_r);      % back to k

Sigma_F = zeros(norb, norb, Nk);
for ik = 1:Nk
    [ix, iy] = ind2sub([Nkx, Nky], ik);
    Sigma_F(:,:,ik) = squeeze(Sigma_F_k(ix,iy,:,:));
end

for ik = 1:Nk
    S = Sigma_F(:,:,ik);
    Sigma_F(:,:,ik) = 0.5 * (S + S');
end

Sigma_orb = Sigma_H + Sigma_F;

end

function f = fermi_function(E, kT)
% FERMI_FUNCTION  Fermi-Dirac 占据数
if kT <= 0
    f = double(E <= 0);
else
    f = 1 ./ (exp(E ./ kT) + 1);
end
end


function [E_tot, E0, E_int] = compute_energy_band(eps_band, Sigma_band, rho_band, kmesh)
% COMPUTE_ENERGY_BAND
% --------------------------------------------------------------
% 在 band 基底中计算 HF 总能量：
%
%   H0(k)   = diag(ε_n(k))       （这里 ε_n 是 *无 seed* 的 H0 带能）
%   E0      = Σ_k w_k Tr[ H0(k) ρ(k) ]
%   E_int   = 0.5 Σ_k w_k Tr[ Σ_HF(k) ρ(k) ]
%   E_tot   = E0 + E_int

[Nb_act, Nk] = size(eps_band);

if isfield(kmesh, 'weight')
    wk = kmesh.weight(:);
else
    wk = ones(Nk,1) / Nk;
end

E0_k   = zeros(Nk,1);
Eint_k = zeros(Nk,1);

for ik = 1:Nk
    rho  = rho_band(:,:,ik);
    occ  = real(diag(rho));
    epsk = eps_band(:,ik);

    E0_k(ik) = sum(epsk .* occ);

    Sig = Sigma_band(:,:,ik);
    Eint_k(ik) = real(trace(Sig * rho));
end

E0    = sum(wk .* E0_k);
E_int = 0.5 * sum(wk .* Eint_k);
E_tot = E0 + E_int;

end



function V_orb0 = build_V_doublegate_from_wpos_singlevalley(g, kmesh, pars)
% BUILD_V_DOUBLEGATE_FROM_WPOS_SINGLEVALLEY
% --------------------------------------------------------------
% 从 Wannier 轨道的 z 位置 g.wpos(:,3) 出发，在 *单 valley* 轨道基下
% 构造双 gate 的库仑相互作用核 V(q; z_a, z_b)：
%
%   gate 在 z' = ±d_gate，样品轨道 z'_a 在中间（自动平移到中心）。
%
% 标准双 gate 形式（2312.11617 类公式）为：
%
%   V(q;z,z') = 2π (e^2/ε_r) / q *
%               [sinh(q(z_<+d)) sinh(q(d-z_>))] / sinh(2qd)
%
% 这里为了数值稳定，把上式改写成只含负指数的等价形式：
%
%   记 z',z'' 为平移到样品中心的坐标，z_< = min(z',z''), z_> = max(...)
%   A = 2qd, α = q(z_<+d), β = q(d-z_>)
%
%   R(q;z,z') = sinhα·sinhβ / sinh(2qd)
%             = ½ e^{-q|z'-z''|} · (1 - e^{-2α})(1 - e^{-2β}) / (1 - e^{-4qd})
%
%   V(q;z,z') = 2π e^2/ε_r · (1/q) · R(q;z,z')
%
% 输出：
%   V_orb0.V_r  [Nkx x Nky x norb0 x norb0]  Fock 用 real-space 核（尚未 /A）
%   V_orb0.V0   [norb0 x norb0]             Hartree 用 q→0 核
%   V_orb0.A    标准化面积 A = Nk * Acell（Å^2）
%
% 注意：
%   - 这是 *单 valley* 的核，后面还要用 lift_V_singlevalley_to_valley
%     复制到 double valley。
%   - 这里完全没有 sinh/cosh，只有 exp(负数)，不会 overflow。
% --------------------------------------------------------------

    %% ---------- 0. 基本尺寸与参数 ----------
    norb0 = size(g.wpos, 1);
    Nkx   = kmesh.Nkx;
    Nky   = kmesh.Nky;
    Nk    = Nkx * Nky;

    % 必要参数
    if ~isfield(pars,'d_gate')
        error('pars.d_gate 未设置（gate 与样品中心的距离，单位 Å）');
    end
    d_gate_in = pars.d_gate;             % 你设定的 gate 距离（Å）

    if ~isfield(pars,'e2_over_eps')
        error('pars.e2_over_eps 未设置，应为 e^2/ε_r [eV·Å]，例如 14.4/ε_r');
    end
    e2_over_eps = pars.e2_over_eps;      % e^2/ε_r [eV·Å]

    if isfield(pars,'q_eps')
        q_eps = pars.q_eps;
    else
        q_eps = 1e-8;                    % 防止除以 q=0
    end

    %% ---------- 1. 把 z 平移到样品中心 ----------
    z_raw    = g.wpos(:,3);              % 原始 z [Å]
    z_min    = min(z_raw);
    z_max    = max(z_raw);
    z_center = 0.5 * (z_min + z_max);    % 样品中心

    z_orb    = z_raw - z_center;         % 现在样品位于大致 [-t/2, t/2]
    half_thick = 0.5 * (max(z_orb) - min(z_orb));

    % 确保 gate 在样品外面：|z| < d_gate
    d_gate = max(d_gate_in, half_thick + 1.0);   % 如有必要自动往外放一点

    %% ---------- 2. 构造与 kmesh 匹配的 q-mesh（中心在 q=0） ----------
    % 从 kmesh.kpts 恢复 k 的 fractional 步长
    kx_frac = reshape(kmesh.kpts(:,1), Nkx, Nky);
    ky_frac = reshape(kmesh.kpts(:,2), Nkx, Nky);

    if Nkx > 1
        dkx_frac = kx_frac(1,2) - kx_frac(1,1);
    else
        dkx_frac = 0.0;
    end
    if Nky > 1
        dky_frac = ky_frac(2,1) - ky_frac(1,1);
    else
        dky_frac = 0.0;
    end

    % 把 q-grid 的"零点"放在网格中心（满足卷积的周期假设）
    [ix_grid, iy_grid] = ndgrid( (0:Nkx-1) - floor(Nkx/2), ...
                                 (0:Nky-1) - floor(Nky/2) );
    qx_frac = ix_grid * dkx_frac;
    qy_frac = iy_grid * dky_frac;

    % fractional → 笛卡尔 q (Å^-1): q = f1*b1 + f2*b2
    B2 = g.b(1:2,1:2);     % 行为 b1, b2
    qx = qx_frac * B2(1,1) + qy_frac * B2(2,1);
    qy = qx_frac * B2(1,2) + qy_frac * B2(2,2);

    q      = sqrt(qx.^2 + qy.^2);
    q_safe = q + q_eps;    % 避免除零（q=0 点专门特别处理）

    % 与 d_gate 相关：A = 2qd, 1 - e^{-4qd} = 1 - e^{-2A}
    A        = 2 * q_safe * d_gate;      % 2qd
    exp_m2A  = exp(-2 * A);              % e^{-4qd}，不会 overflow
    denA     = 1 - exp_m2A;              % 1 - e^{-4qd}
    denA(denA == 0) = 1.0;               % 理论上只有 q=0 极限才会 0，这里避免 0/0

    %% ---------- 3. 对每一对轨道 (a,b) 构造 V(q;z_a,z_b) ----------
    V_r = zeros(Nkx, Nky, norb0, norb0);
    V0  = zeros(norb0, norb0);

    for a = 1:norb0
        za = z_orb(a);           % 已经平移后的 z_a
        for b = 1:norb0
            zb = z_orb(b);

            zmin = min(za, zb);
            zmax = max(za, zb);

            % α = q (z_< + d), β = q (d - z_>)
            alpha = q_safe .* (zmin + d_gate);
            beta  = q_safe .* (d_gate - zmax);

            % 只含负指数：exp(-2α), exp(-2β) 均在 (0,1]
            exp_m2alpha = exp(-2 * alpha);
            exp_m2beta  = exp(-2 * beta);

            % F = (1 - e^{-2α})(1 - e^{-2β}) / (1 - e^{-4qd})
            numF = (1 - exp_m2alpha) .* (1 - exp_m2beta);
            F    = numF ./ denA;         % 不会出现 Inf，因为 denA 已避免 0

            % ratio = ½ e^{-q|z-z'|} * F
            ratio = 0.5 * exp(- q_safe .* abs(za - zb)) .* F;

            % V(q;z,z') = 2π (e^2/ε_r) * ratio / q
            Vq = (2*pi*e2_over_eps) * ratio ./ q_safe;

            % Fock 不需要 q=0 模式：直接设为 0，Hartree 单独由 V0 处理
            % Vq(q == 0) = 0.0;

            % real-space 核：V_r = FFT[V(q)]
            V_r(:,:,a,b) = fft2(Vq);

            % Hartree：用最小一圈非零 q 的平均值作为 V(q→0;z_a,z_b)
            q_nonzero = q_safe(q_safe > q_eps);
            if ~isempty(q_nonzero)
                q_min = min(q_nonzero);
                mask_small = (q_safe <= 1.5*q_min) & (q_safe > q_eps);
                vals = real(Vq(mask_small));
                if ~isempty(vals)
                    V0(a,b) = mean(vals);
                end
            end
        end
    end

    %% ---------- 4. 面积归一化 ----------
    % Acell 使用 g.a 前 2×2 的行列式
    % Acell = abs(det(g.a(1:2,1:2)));   % Å^2
    w_uniform = mean(kmesh.weight);
    A_eff = 1/w_uniform;

    V_orb0.V_r = V_r;                 % 注意：还没除以 A，留给 HF 里统一处理
    V_orb0.V0  = V0;                  % Hartree q→0 核
    V_orb0.A   = Nk * A_eff;          % 约定：A = Nk * Acell

    %% ---------- 5. 最后保险：清理非有限值（理论上不会触发） ----------
    V_orb0.V_r(~isfinite(V_orb0.V_r)) = 0.0;
    V_orb0.V0(~isfinite(V_orb0.V0))   = 0.0;
end




function V_full = lift_V_singlevalley_to_valley(V_orb0)
% LIFT_V_SINGLEVALLEY_TO_VALLEY
% --------------------------------------------------------------
% 从单 valley 轨道基 Coulomb 核 V_orb0 生成 double-valley 核：
%
% 单 valley:
%   V_orb0.V0  : [norb0 x norb0]
%   V_orb0.V_r : [Nkx x Nky x norb0 x norb0]
%   V_orb0.A   : Nk * Acell
%
% double valley (v = 1:K, 2:K'):
%   a = (v-1)*norb0 + α, α=1..norb0
%
% 约定：
%   - Hartree 完全 valley-blind:
%       V0_full((v,α),(v',β)) = V0_0(α,β) (任意 v,v')
%   - Fock 只 intravalley:
%       V_r_full((v1,α),(v2,β)) = V_r0(α,β) if v1==v2 else 0

    V0_0  = V_orb0.V0;
    V_r0  = V_orb0.V_r;
    A     = V_orb0.A;

    [Nkx, Nky, norb0, norb0_2] = size(V_r0);
    if norb0 ~= norb0_2
        error('V_orb0.V_r must be [Nkx x Nky x norb0 x norb0].');
    end

    norb = 2 * norb0;

    % Hartree: valley-blind
    V0_full = zeros(norb, norb);
    for v1 = 1:2
        for v2 = 1:2
            for a = 1:norb0
                for b = 1:norb0
                    ia = (v1-1)*norb0 + a;
                    jb = (v2-1)*norb0 + b;
                    V0_full(ia,jb) = V0_0(a,b);
                end
            end
        end
    end

    % Fock: intravalley only
    V_r_full = zeros(Nkx, Nky, norb, norb);
    for v1 = 1:2
        for v2 = 1:2
            for a = 1:norb0
                for b = 1:norb0
                    ia = (v1-1)*norb0 + a;
                    jb = (v2-1)*norb0 + b;
                    if v1 == v2
                        V_r_full(:,:,ia,jb) = V_r0(:,:,a,b);
                    else
                        V_r_full(:,:,ia,jb) = 0.0;
                    end
                end
            end
        end
    end

    V_full = struct();
    V_full.V0  = V0_full;
    V_full.V_r = V_r_full;
    V_full.A   = A;
end

function H_seed_orb = build_seed_orbital(g, seed)
% BUILD_SEED_ORBITAL
% --------------------------------------------------------------
% 在 double-valley 轨道基 (valley⊗orbital) 上构造 seed Hamiltonian:
%
%   H_seed_orb: [norb x norb], norb = 2*norb0
%
% 支持：
%   seed.type = 'none'          → 返回 []
%   seed.type = 'valley'        → H_seed = amp * τ_z
%   seed.type = 'custom_orbital'→ 用户提供 seed.H_orbital
%
% 注意：
%   - norb0 = size(g.ham,1) 是 *单 valley* 轨道数；
%   - double valley 基底假定为 [valley K 的 norb0 轨道; valley K' 的 norb0 轨道]。
%   - H_seed_orb 只用于生成初始密度矩阵，HF 自洽过程中不再显式出现。

hamR0 = g.ham;
norb0 = size(hamR0, 1);
norb  = 2 * norb0;

H_seed_orb = [];

if ~isfield(seed, 'type')
    return;
end
if strcmpi(seed.type, 'none')
    return;
end

if ~isfield(seed, 'amp')
    amp = 0.0;
else
    amp = seed.amp;
end

switch lower(seed.type)
    case 'valley'
        if amp == 0
            H_seed_orb = [];
            return;
        end
        % tau_z = diag([ +ones(norb0,1); -ones(norb0,1) ]);
        sigma_z=[1,0;0,-1];
        layer_z=diag([1,1,0.0,0.0,0,0,0,0,-1,-1]);
        valley_z=kron(sigma_z,layer_z);
        sigma_0=[1,0;0,1];
        tau_z = kron(sigma_0,valley_z);
        H_seed_orb = amp * tau_z;

    case 'custom_orbital'
        if ~isfield(seed, 'H_orbital')
            error('build_seed_orbital: custom_orbital requires seed.H_orbital.');
        end
        H = seed.H_orbital;
        if ~isequal(size(H), [norb, norb])
            error('build_seed_orbital: H_orbital must be [%d x %d].', norb, norb);
        end
        H_seed_orb = 0.5 * (H + H');   % Hermitian

    otherwise
        warning('build_seed_orbital: unknown seed.type=%s, treat as none.', seed.type);
        H_seed_orb = [];
end

end

function kmesh = build_kmesh_K_patch(g, Nkx, Nky, center_frac, delta_frac)
% BUILD_KMESH_K_PATCH
% --------------------------------------------------------------
% 构造以 K = (1/3, 2/3) 为中心的小 patch 的 k-mesh（fractional 坐标），
% 用于 valley-patch HFMF。在这个 patch 上，我们之后会构造
% double-valley Hamiltonian:
%
%   H_full(k) = diag( H_single(k), conj(H_single(-k)) ),
%
% 这样 K 和 K' 两个 valley 都已经自动包含进来了。
%
% INPUT
%   g.a          [3 x 3]   实空间基矢（Å），只用前 2×2 求 Acell
%
%   Nkx, Nky     标量      patch 网格在 kx, ky 方向上的点数
%
%   center_frac  [1 x 2]   patch 中心（fractional），默认 [1/3, 2/3]
%                          例如 K 点为 (1/3, 2/3)
%
%   delta_frac   标量或 [1 x 2]
%                每个方向的半宽度（fractional），即 k ∈ center ± delta。
%                若为标量，则 x,y 同用该值（例如 0.15）。
%
% OUTPUT
%   kmesh 结构体：
%     .kpts        [Nk x 2] fractional k 点 (kx, ky)，flatten 顺序为 (ix,iy)
%     .Nkx, .Nky   网格维度
%     .weight      [Nk x 1] 物理 k 权重 w_k = patch_fraction / (Acell * Nk)
%                           满足
%                               sum_k w_k = patch_fraction / Acell
%                           其中 patch_fraction = Δkx_frac * Δky_frac
%
%     .center_frac [1 x 2]   记录中心
%     .delta_frac  [1 x 2]   记录半宽
%
% 说明：
%   - fractional 坐标是相对于 g.b 的基矢，H(k) 用 e^{i 2π k·R} 构造；
%   - patch_fraction = (kx_max - kx_min) * (ky_max - ky_min)，
%     对应 patch 占整个 BZ（fractional [0,1)×[0,1)）的面积比例；
%   - 对整 BZ 的均匀网格，w_k = 1/(Acell * Nk)；这里 patch 只占一部分，
%     所以 w_k = patch_fraction/(Acell * Nk)。
% --------------------------------------------------------------

    if nargin < 4 || isempty(center_frac)
        center_frac = [0, 0];
    end
    if nargin < 5 || isempty(delta_frac)
        delta_frac = 0.15;
    end

    % 允许 delta_frac 是标量或 [dx, dy]
    if isscalar(delta_frac)
        dx = delta_frac;
        dy = delta_frac;
    else
        dx = delta_frac(1);
        dy = delta_frac(2);
    end

    % ---- 1. 在 fractional 坐标中定义 patch 区间 ----
    Kx0 = center_frac(1);
    Ky0 = center_frac(2);

    kx_min = Kx0 - dx;
    kx_max = Kx0 + dx;
    ky_min = Ky0 - dy;
    ky_max = Ky0 + dy;

    % patch 在 fractional 单元中的面积比例
    frac_x = kx_max - kx_min;   % = 2*dx
    frac_y = ky_max - ky_min;   % = 2*dy
    patch_fraction = frac_x * frac_y;

    % ---- 2. 生成 1D 网格并 wrap 到 [0,1) ----
    kx_line = linspace(kx_min, kx_max, Nkx);
    ky_line = linspace(ky_min, ky_max, Nky);

    % wrap 到 [0,1)
    kx_line = mod(kx_line, 1.0);
    ky_line = mod(ky_line, 1.0);

    % ---- 3. 生成 2D 网格并 flatten 为 [Nk x 2] ----
    [kx_grid, ky_grid] = meshgrid(kx_line, ky_line);
    Nk = Nkx * Nky;

    kpts = [kx_grid(:), ky_grid(:)];   % fractional coords

    % ---- 4. 计算权重 w_k = ΔA_patch / (2π)^2 / Nk ----
    % BZ 面积 A_BZ = |det(b1,b2)| = (2π)^2 / Acell
    % patch 面积 A_patch = patch_fraction * A_BZ
    % 所以
    %   w_k = A_patch / (2π)^2 / Nk = patch_fraction / (Acell * Nk)
    %
    % Acell 用 g.a 的前 2×2 求行列式
    Acell = abs(det(g.a(1:2,1:2)));   % Å^2
    w_k   = patch_fraction / (Acell * Nk);
    weight = w_k * ones(Nk,1);

    % ---- 5. 打包输出 ----
    kmesh.kpts        = kpts;          % [Nk x 2] fractional
    kmesh.Nkx         = Nkx;
    kmesh.Nky         = Nky;
    kmesh.weight      = weight;
    kmesh.center_frac = center_frac(:).';
    kmesh.delta_frac  = [dx, dy];

    %（需要的话可以在这里顺便存一下原始的 min/max，方便 debug）
    kmesh.kx_range = [kx_min, kx_max];
    kmesh.ky_range = [ky_min, ky_max];
end


function mu = solve_mu_band(E_HF, kT, kmesh, Nb_act, opts)
% SOLVE_MU_BAND
% --------------------------------------------------------------
% 根据当前 HF 光谱 E_HF(n,k) 自洽求解化学势 μ，使
%
%   N(μ) = Σ_k w_k Σ_n f(E_HF(n,k) - μ, kT) = N_target
%
% 其中：
%   - w_k = kmesh.weight(ik)
%   - 默认 N_target = (Nb_act/2) * Σ_k w_k  （半填充）
%   - 若 opts.N_target 存在，则用 opts.N_target 作为目标粒子数
%
% 注意：
%   - 这里的 N(μ) 是 "每晶胞的电子数"（如果 w_k 按照你在 build_kmesh_K_patch
%     里那样定义：w_k = patch_fraction / (Acell * Nk)）。
%   - 若你希望指定特定的掺杂密度，可以自己设 opts.N_target。
% --------------------------------------------------------------

    [~, Nk] = size(E_HF);

    % k 权重
    if isfield(kmesh, 'weight')
        wk = kmesh.weight(:);
    else
        wk = ones(Nk,1) / Nk;
    end

    % 目标粒子数：默认半填充；也可以在 opts.N_target 中显式给出
    if isfield(opts, 'N_target')
        N_target = opts.N_target;
    else
        N_target = (Nb_act/2) * sum(wk);  % half-filling of active manifold
    end

    % 定义数目方程 N(μ) - N_target = 0
    Emin = min(E_HF, [], 'all');
    Emax = max(E_HF, [], 'all');

    % 给 fzero 一个比较宽的搜索区间
    dE = max(Emax - Emin, 1.0);
    mu_low  = Emin - 5 * dE;
    mu_high = Emax + 5 * dE;

    f = @(mu_try) number_diff_mu(mu_try, E_HF, kT, wk, N_target);

    mu = fzero(f, [mu_low, mu_high]);
end


function val = number_diff_mu(mu_try, E_HF, kT, wk, N_target)
% NUMBER_DIFF_MU
%   fzero 用的数目方程：N(μ) - N_target

    [Nb_act, Nk] = size(E_HF);

    N_mu = 0.0;
    for ik = 1:Nk
        occ_k = fermi_function(E_HF(:, ik) - mu_try, kT);  % [Nb_act x 1]
        N_mu = N_mu + wk(ik) * sum(occ_k);
    end

    val = N_mu - N_target;
end











