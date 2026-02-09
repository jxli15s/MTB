%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Construct the g.ham                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("Rgra_5s");
a=2.46;
g.a=[sqrt(3)/2,-1/2,0;0,1,0;0,0,1]*a;
% 倒格矢（标准 2D：b = 2π * (A^{-1})^T）
g.b = 2*pi * inv(g.a).';
%%
%% demo_band_5LG.m
% 简单画一下 5-layer graphene 两带有效模型的能带
clear; clc;

% ---- 1. 基本参数 ----
Nlayer = 5;
SWM    = get_SWM_params(Nlayer);

% 有效模型中的 k 是无量纲 (k_phys * a_lat)
kmax = 0.4;       % 无量纲的最大动量（可调）
Nk   = 401;
kxs  = linspace(-kmax, kmax, Nk);
ky   = 0.0;       % 沿 kx 方向扫

Evals = zeros(2, Nk);

% ---- 2. 扫描 kx，计算 2×2 有效模型本征值 ----
for ik = 1:Nk
    kx = kxs(ik);
    H2 = H2B_5LG(kx, ky, +1, SWM, Nlayer);   % valley K (q=+1), 单自旋
    [V, D] = eig(H2);
    e = sort(real(diag(D)));                 % 排序
    Evals(:, ik) = e;
end

% ---- 3. 作图 ----
figure;
plot(kxs, Evals(1, :), 'LineWidth', 1.5); hold on;
plot(kxs, Evals(2, :), 'LineWidth', 1.5);
xlabel('k_x (dimensionless)');
ylabel('Energy (eV)');
title(sprintf('N = %d-layer graphene, 2-band effective model (valley K)', Nlayer));
legend('Band 1','Band 2','Location','Best');
grid on;
%%
% 1. SWM 参数 & 有效模型参数
Nlayer    = 5;
SWM       = get_SWM_params(Nlayer);
pars.N    = Nlayer;
pars.SOC  = 0.0;
pars.SOC_dir    = 0;
pars.SOC_single = 1;

% 2. k-mesh (无量纲 1/a)
dk   = 2.0e-3;
kmax = 0.15;
k1   = linspace(-kmax, kmax, round(2*kmax/dk)+1);
[kx,ky] = meshgrid(k1,k1);
nk   = numel(k1);

% 3. V_r, hf 结构体 (这里你用 VS/VD + fft2 构造好 V_r)
hf.a_lat     = 2.46;         % Å
hf.A       = (2*pi* hf.a_lat/ dk)^2;   % 只是示意, 你可以用之前的 (2π a /dk)^2
hf.beta    = 5.8e4;      % Python 里 beta~58 1/meV, 这里 /1000 -> 1/eV ~0.2K
hf.ke      = 14.4;         % eV·Å
hf.er      = 27;
hf.N_layer = Nlayer;
hf.d_gate  = 369;          % Å
hf.d_lay   = 3.35;         % Å
hf.alp     = 0.03;         % 或者你想要的 alp
% hf.alp     = 0.0;         % 或者你想要的 alp

V_r = build_Vr_5LG(kx, ky, hf);  % [nk x nk x 8 x 8]
hf.V_r     = V_r;          % [nk x nk x 8 x 8]

% 4. seed_V (按照你 LAFz 的方式构造一个 8x8 的矩阵)
sig0 = eye(2);
sigx = [0 1; 1 0];
sigz = [1 0; 0 -1];
m_LAF = 0.000;  % 对应原代码的 5 meV

K_up  = +m_LAF;
Kp_up = +m_LAF;
K_dn  = -m_LAF;
Kp_dn = -m_LAF;

vx  = 0;       % 只要 LAF_z，相干项先关掉
sx  = 0;
vsx = 0;
% 8-band operators: spin × valley × orbital
eye8 = kron(sig0, kron(sig0, sig0));
S_x  = kron(sigx, kron(sig0, sig0));
S_z  = kron(sigz, kron(sig0, sig0));
V_x  = kron(sig0, kron(sigx, sig0));
V_z  = kron(sig0, kron(sigz, sig0));
L_x  = kron(sig0, kron(sig0, sigx));
L_z  = kron(sig0, kron(sig0, sigz));

% ==== 构造 seed_V: 与原 Python 完全同构 ====
seed_V = ...
    + K_up  * (eye8 + V_z) * (eye8 + S_z) * L_z / 4 ...
    + K_dn  * (eye8 + V_z) * (eye8 - S_z) * L_z / 4 ...
    + Kp_up * (eye8 - V_z) * (eye8 + S_z) * L_z / 4 ...
    + Kp_dn * (eye8 - V_z) * (eye8 - S_z) * L_z / 4 ...
    + vx    * V_x ...
    + sx    * S_x ...
    + vsx   * (V_x * S_x);

seed_V = m_LAF * (S_z * L_z);   % valley-identity 的 LAF_z

% seed_V = zeros(8);
% ... 按照 (eye+V_z)(eye+S_z)L_z/4 那套构造 ...

% 5. opts
opts.U_ext    = 0.0;        % 暂时先不用外场
opts.ne       = 0.0e12;     % 零掺杂
opts.seed_V   = seed_V;
opts.mix      = 0.3;
opts.max_iter = 200;
opts.tol      = 1e-6;

% 6. 自洽
res = hfmf_single_point(kx, ky, SWM, pars, hf, opts);
%%

figure()
for i =1:8
    hold on;
surf(squeeze(res.E(:,:,i)))
end
%%
figure()
for j=1:nk
for i =1:8
hold on;
plot(res.E(j,:,i))
end
end






%%
function result = hfmf_single_point(kx, ky, SWM, pars, hf, opts)
%HFMF_SINGLE_POINT  单个 (U_ext, ne) 点的 Hartree-Fock 自洽
%
%   result = hfmf_single_point(kx, ky, SWM, pars, hf, opts)
%
% 输入：
%   kx, ky : [nk x nk]  无量纲动量 (单位 1/a)
%   SWM    : get_SWM_params(N) 的输出 (单位 eV)
%   pars   : 结构体, 至少包含
%            .N           层数
%            .SOC
%            .SOC_dir
%            .SOC_single
%   hf     : HF 参数结构体, 至少包含
%            .V_r    [nk x nk x N_band x N_band] Coulomb 核 (r-space)
%            .A      系统面积 (与 ne 的单位匹配)
%            .beta   1/(k_B T) [1/eV]
%            .ke, .er, .N_layer, .d_gate, .d_lay, .alp   (给 V_HF_builder 用)
%   opts   : 自洽选项
%            .U_ext      外场 U (倍乘 H_U, 这里先设 H_U = 0)
%            .ne         电子密度 (与 hf.A 一致的单位, 例如 e/Å^2)
%            .seed_V     [N_band x N_band] 初始 seed 势矩阵 (LAFz 等)
%            .mix        mixing 参数 (0~1)
%            .max_iter   最大迭代次数
%            .tol        收敛阈值
%
% 输出：
%   result 结构体：
%     .rho      [nk x nk x N_band x N_band] 自洽后的密度矩阵
%     .mu       收敛后的化学势
%     .E        [nk x nk x N_band] 本征值
%     .Evec     [nk x nk x N_band x N_band] 本征态
%     .E_tot    总能量 (包含 HF)
%     .E_tot0   非相互作用能量
%     .error    最终收敛误差
%     .iter     实际迭代步数

    % -------------------- 基本尺寸 --------------------
    [nkx, nky] = size(kx);
    if nkx ~= nky
        error('当前实现假定 nkx = nky 的方形 k-mesh');
    end
    nk = nkx;

    % 先看一个 k 点，确定带数
    Htest   = H0_5LG_eff_cart(kx(1,1), ky(1,1), SWM, pars);
    N_band  = size(Htest, 1);

    % -------------------- 构造 H_0(k) --------------------
    H0_k = complex(zeros(nk, nk, N_band, N_band));
    for ix = 1:nk
        for iy = 1:nk
            H0_k(ix,iy,:,:) = H0_5LG_eff_cart(kx(ix,iy), ky(ix,iy), SWM, pars);
        end
    end

    % -------------------- 外场 H_U(k)（目前先设为 0） --------------------
    H_U_k = zeros(size(H0_k));

    % 如果你以后要加上 python 里那一坨 r_gap 的 H_U, 可以在这里填：
    % H_U_k(ix,iy,:,:) = ...  (照 python 的公式写 r_gap+r_gap2)/r_N 那套)

    U_ext  = opts.U_ext;          % 外场 U
    H_n_k  = H0_k + U_ext * H_U_k;  % 非相互作用 part (不含 seed_V 和 HF)

    % -------------------- seed_V 只用于初始态 --------------------
    seed_V = opts.seed_V;         % [N_band x N_band]
    if isempty(seed_V)
        seed_V = zeros(N_band);
    end

    % 利用 implicit expansion: [nk nk N N] + [1 1 N N]
    H_seed_k = H_n_k + reshape(seed_V, 1,1,N_band,N_band);

    % -------------------- 初始对角化 (H_seed) --------------------
    E_seed    = zeros(nk, nk, N_band);
    Evec_seed = complex(zeros(nk, nk, N_band, N_band));

    for ix = 1:nk
        for iy = 1:nk
            H = squeeze(H_seed_k(ix,iy,:,:));   % [N_band x N_band]
            [V,D] = eig(H);
            e  = real(diag(D));
            [e_sorted, idx] = sort(e, 'ascend');
            V  = V(:, idx);
            E_seed(ix,iy,:)    = e_sorted;
            Evec_seed(ix,iy,:,:) = V;
        end
    end

    % -------------------- 初始 rho_old = rho_new --------------------
    ne   = opts.ne;
    beta = hf.beta;
    A    = hf.A;

    mu0      = get_mu_mat(E_seed, ne, beta, A, N_band);
    occ_seed = fermi_dirac_mat(E_seed, mu0, beta);
    rho_old  = get_rho_mat(Evec_seed, occ_seed);
    rho_new  = rho_old;

    % -------------------- 开始自洽迭代 --------------------
    max_iter = opts.max_iter;
    tol      = opts.tol;
    mix      = opts.mix;

    E_int    = E_seed;
    Evec_int = Evec_seed;
    mu       = mu0;

    for it = 1:max_iter
        % 1. mixing 密度矩阵
        rho_mix = mix * rho_new + (1 - mix) * rho_old;

        % 2. 由 rho_mix 构造 HF 势 Σ_HF(k)
        V_F_k = V_HF_builder_mat(rho_mix, hf);   % [nk x nk x N_band x N_band]

        % 3. 有效 Hamiltonian: H_int(k) = H_n(k) + Σ_HF(k)
        H_int_k = H_n_k + V_F_k;

        % 4. 对每个 k 对角化 H_int(k)
        for ix = 1:nk
            for iy = 1:nk
                H = squeeze(H_int_k(ix,iy,:,:));
                [V,D] = eig(H);
                e  = real(diag(D));
                [e_sorted, idx] = sort(e, 'ascend');
                V  = V(:, idx);
                E_int(ix,iy,:)    = e_sorted;
                Evec_int(ix,iy,:,:) = V;
            end
        end

        % 5. 更新化学势 & 密度矩阵
        mu       = get_mu_mat(E_int, ne, beta, A, N_band);
        occ_int  = fermi_dirac_mat(E_int, mu, beta);
        rho_old  = rho_new;
        rho_new  = get_rho_mat(Evec_int, occ_int);

        % 6. 能量和误差
        [E_tot, E_tot0] = get_total_energies(H_int_k, H_n_k, rho_new);

        err_sum = 0.0;
        for ix = 1:nk
            for iy = 1:nk
                diff = squeeze(rho_new(ix,iy,:,:) - rho_old(ix,iy,:,:));
                err_sum = err_sum + norm(diff, 'fro');
            end
        end
        err = err_sum / ( (N_band * nk)^2 );

        fprintf('Iter %3d: err = %.3e, E_tot = %.6f eV, mu = %.6f eV\n', ...
                it, err, E_tot, mu);

        if err < tol
            break;
        end
    end

    % -------------------- 输出结果 --------------------
    result.rho    = rho_new;
    result.mu     = mu;
    result.E      = E_int;
    result.Evec   = Evec_int;
    result.E_tot  = E_tot;
    result.E_tot0 = E_tot0;
    result.error  = err;
    result.iter   = it;
end


% ==================== 辅助子函数 ====================

function nF = fermi_dirac_mat(E, mu, beta)
%FERMI_DIRAC_MAT  Fermi-Dirac 占据 n_F(E,mu)
% E,mu,beta 单位一致 (eV)
    nF = 1 ./ (1 + exp(beta * (E - mu)));
end


function rho_k = get_rho_mat(Evec, occ)
%GET_RHO_MAT  根据本征态和占据数构造 rho(k)
%
%   Evec: [nk x nk x N_band x N_band], 每个 (k) 列是本征态 |ψ_n(k)>
%   occ : [nk x nk x N_band], Fermi 占据 f_{n,k}
%
%   rho_k(ix,iy,:,:) = sum_n f_{n,k} |ψ_n(k)><ψ_n(k)|

    [nkx, nky, N_band, ~] = size(Evec);
    rho_k = complex(zeros(nkx, nky, N_band, N_band));
    for ix = 1:nkx
        for iy = 1:nky
            V = squeeze(Evec(ix,iy,:,:));     % [N_band x N_band]
            f = squeeze(occ(ix,iy,:));        % [N_band]
            rho_k(ix,iy,:,:) = V * (diag(f) * V');   % V diag(f) V^\dagger
        end
    end
end


function mu = get_mu_mat(E, ne, beta, A, N_band)
%GET_MU_MAT  给定能谱 E(k,n) 和目标密度 ne, 求解化学势 mu
%
% 输入：
%   E   : [nk x nk x N_band]
%   ne  : 电子密度 (e/面积)
%   beta: 1/(k_B T)
%   A   : 系统总面积
%   N_band: 带数
%
% 方程：
%   sum_{k,n} f(E_{k,n}-mu) - nk^2 * (N_band/2) - ne*A = 0

    [nkx, nky, ~] = size(E);
    Nk_tot = nkx * nky;

    % 中心点（大致在 gap 中央）
    E_valence_max = max(E(:,:,1:floor(N_band/2)), [], 'all');
    E_cond_min    = min(E(:,:,floor(N_band/2)+1:end), [], 'all');
    mu_center     = 0.5 * (E_valence_max + E_cond_min);

    mu_low  = mu_center - 100;   % eV
    mu_high = mu_center + 100;   % eV

    % 定义数值方程
    f = @(mu) find_mu(mu, E, ne, beta, A, Nk_tot, N_band);

    mu = fzero(f, [mu_low, mu_high]);
end


function val = find_mu(mu, E, ne, beta, A, Nk_tot, N_band)
    occ = 1 ./ (1 + exp(beta * (E - mu)));
    val = sum(occ, 'all') - Nk_tot * (N_band/2) - ne * A;
end


function [E_tot, E_tot0] = get_total_energies(H_int_k, H_n_k, rho_k)
%GET_TOTAL_ENERGIES  计算总能量与非相互作用能量
%
%   E_tot  = sum_k Tr[ (H_int + H_n)/2 * rho ]
%   E_tot0 = sum_k Tr[ H_n * rho ]

    [nkx, nky, ~, ~] = size(H_int_k);
    E_tot  = 0.0;
    E_tot0 = 0.0;

    for ix = 1:nkx
        for iy = 1:nky
            Hn   = squeeze(H_n_k(ix,iy,:,:));
            Hint = squeeze(H_int_k(ix,iy,:,:));
            rho  = squeeze(rho_k(ix,iy,:,:));
            E_tot  = E_tot  + real(trace(0.5 * (Hint + Hn) * rho));
            E_tot0 = E_tot0 + real(trace(Hn * rho));
        end
    end
end


function V_r = build_Vr_5LG(kx, ky, hf)
%BUILD_VR_5LG  构造 5-layer graphene 的 Coulomb 核 V_r(kx,ky,band,band)
%
%   V_r = build_Vr_5LG(kx, ky, hf)
%
% 输入：
%   kx, ky : [nk x nk]  无量纲动量 (k = k_phys * a_lat)
%   hf     : 结构体，至少包含
%            .ke      库仑常数 (eV·Å)，建议用 14.4
%            .er      介电常数 ε_r
%            .a_lat   晶格常数 a_lat (Å)，比如 2.46
%            .d_gate  gate 距离 (Å)
%            .d_lay   层间距 (Å)，比如 3.35
%            .N_layer 层数 N (2–6)
%
% 输出：
%   V_r : [nk x nk x N_band x N_band]，这里 N_band = 8
%         (spin × valley × 2 orbital)

    [nkx, nky] = size(kx);
    if nkx ~= nky
        error('build_Vr_5LG: 当前实现假定 nkx = nky 的方形 k-mesh');
    end
    nk = nkx;

    ke    = hf.ke;        % eV·Å
    er    = hf.er;
    a_lat = hf.a_lat;     % Å
    d_gate= hf.d_gate;    % Å
    d_lay = hf.d_lay;     % Å
    N     = hf.N_layer;

    % ---------- 1. |q| (Å^-1)，注意 kx,ky 是无量纲 (k_phys * a_lat) ----------
    q = sqrt(kx.^2 + ky.^2) ./ a_lat;  % [Å^-1]
    q = q + 1e-6;                      % 避免 q=0 时 1/q 爆掉

    % ---------- 2. same-layer / different-layer 的 V(q) 核 ----------
    % 对应 Python 中：
    %  V_S_r = 2π ke/er * fftn( 1/qq * ( cosh(2qq d_gate) - cosh(qq (N-1) d_lay) ) / sinh(2qq d_gate) )
    %  V_D_r = 2π ke/er * fftn( 1/qq * ( cosh(2qq d_gate - qq (N-1) d_lay) - 1 ) / sinh(2qq d_gate) )

    kernel_S = (1 ./ q) .* ...
        (cosh(2 * q * d_gate) - cosh(q * (N-1) * d_lay)) ./ ...
         sinh(2 * q * d_gate);

    kernel_D = (1 ./ q) .* ...
        (cosh(2 * q * d_gate - q * (N-1) * d_lay) - 1) ./ ...
         sinh(2 * q * d_gate);

    V_S_r = 2 * pi * ke / er * fft2(kernel_S);   % [eV·Å^2]
    V_D_r = 2 * pi * ke / er * fft2(kernel_D);   % [eV·Å^2]

    % ---------- 3. orbital 2×2 结构: [[V_S, V_D],[V_D, V_S]] ----------
    Norb = 2;
    V_orb = zeros(nk, nk, Norb, Norb);
    V_orb(:,:,1,1) = V_S_r;
    V_orb(:,:,2,2) = V_S_r;
    V_orb(:,:,1,2) = V_D_r;
    V_orb(:,:,2,1) = V_D_r;

    % ---------- 4. 扩展到 4 个 flavor (spin×valley)：N_band = 4 × 2 ----------
    N_flavor = 4;
    N_band   = N_flavor * Norb;   % 8

    V_r = zeros(nk, nk, N_band, N_band);

    for f1 = 0:N_flavor-1
        for f2 = 0:N_flavor-1
            for o1 = 1:Norb
                for o2 = 1:Norb
                    i = f1 * Norb + o1;
                    j = f2 * Norb + o2;
                    V_r(:,:,i,j) = V_orb(:,:,o1,o2);
                end
            end
        end
    end
end




%%
function V_F_k = V_HF_builder_mat(rho_k, hf)
%V_HF_BUILDER_MAT  用 FFT 构造 Hartree-Fock 势 V_F(k)
%
%   V_F_k = V_HF_builder_mat(rho_k, hf)
%
% 输入：
%   rho_k : [nk x nk x N_band x N_band]  频域密度矩阵 rho(k)
%           （对每个 k 是一个 N_band×N_band 矩阵）
%   hf    : 结构体，至少包含
%           .V_r     [nk x nk x N_band x N_band]  事先构造好的 Coulomb 核 (r-space)
%           .A       标量，系统面积
%           .alp     valley interchange strength
%           .ke, .er Coulomb 参数
%           .N_layer 层数 N
%           .d_gate, .d_lay   距离参数（单位与 V_r 构造时一致）
%
% 输出：
%   V_F_k : [nk x nk x N_band x N_band]  Hartree-Fock 势矩阵 Σ(k)
%
% 逻辑：
%   1. Fock:   Σ_F(k) = - sum_q V(q) rho(k-q)  (用 FFT 做卷积)
%   2. VI:     近似 q~0 的 valley-interchange 项 (alp 控制强度)
%   3. Hartree:近似 q~0 的 same-layer/different-layer Hartree 项
%   4. 把 2+3 加到每个 k 上 (k 无关)

    % -------------------- 基本信息 --------------------
    [nkx, nky, N_band, ~] = size(rho_k);
    if nkx ~= nky
        error('当前实现假定 nkx = nky 的正方形 k-mesh');
    end
    nk = nkx;

    V_r   = hf.V_r;     % (nk x nk x N_band x N_band)
    A     = hf.A;
    N     = hf.N_layer;
    ke    = hf.ke;
    er    = hf.er;
    dgate = hf.d_gate;
    dlay  = hf.d_lay;
    alp   = hf.alp;

    % -------------------- 1. Fock 项（用 FFT 做卷积） --------------------
    V_F_k = complex(zeros(nkx, nky, N_band, N_band));

    % 对每个 band (a,b) 做 2D FFT/IFT
    for a = 1:N_band
        for b = 1:N_band
            rho_r_ab = fft2(rho_k(:,:,a,b));           % FFT over (kx,ky)
            V_r_ab   = V_r(:,:,a,b);                   % 对应的 Coulomb 核
            V_F_r_ab = - rho_r_ab .* V_r_ab / A;       % Σ_F(r) = - V(r) * ρ(r) / A
            V_F_k(:,:,a,b) = ifftshift(ifft2(V_F_r_ab));          % 回到 k 空间
        end
    end

    % -------------------- 2. 计算 k-summed 密度矩阵 --------------------
    % rho_k_sum = sum over k of rho_k(k), shape (N_band, N_band)
    rho_k_sum = squeeze(sum(sum(rho_k, 1), 2));  % [N_band x N_band]

    % -------------------- 3. valley-interchange (VI) 项 --------------------
    % 对应 Python 中：
    % tV_VI = - alp * (2π ke/er*(d_gate - ((N-1)*d_lay)^2/(4*d_gate))) * rho_k_sum / A
    % 然后复制到 spin/valley/flavor 上，再做一次 block 互换

    V_VI = zeros(N_band, N_band);
    if abs(alp) > 0
        % 常数因子 (q -> 0)
        C_VI = 2*pi*ke/er * (dgate - ((N-1)*dlay)^2/(4*dgate));
        tV_VI = - alp * C_VI * rho_k_sum / A;   % [N_band x N_band]

        % 复制到 4 个 flavor × 2 orbital
        % 原 Python: tV_VI = tV_VI * kron(ones([4,4]), eye(2))
        % 这里假定 N_band = 8 = 4 flavor × 2 orbital
        if N_band ~= 8
            warning('当前 VI 实现假定 N_band = 8 (4 flavor × 2 orbital)');
        end
        tV_VI = tV_VI .* kron(ones(4,4), eye(2));

        % block 互换 (valley K ↔ K')
        N_BL  = N_band / 4;  % 每个 flavor 的带数 (这里=2)
        rows1 = [1:N_BL,       2*N_BL+1:3*N_BL];   % blocks 0 & 2
        rows2 = [N_BL+1:2*N_BL,3*N_BL+1:4*N_BL];   % blocks 1 & 3

        % 对角块拷贝
        V_VI(rows1, rows1) = tV_VI(rows2, rows2);
        V_VI(rows2, rows2) = tV_VI(rows1, rows1);
    end

    

    % -------------------- 4. Hartree 项 (q ~ 0 same/diff layer) --------------------
    % 完全照 Python:
    %
    % diag_D[::2]  = trace(rho_k_sum[1::2,1::2]).real
    % diag_D[1::2] = trace(rho_k_sum[::2,::2]).real
    %
    % diag_S[1::2] = diag_D[::2]
    % diag_S[::2]  = diag_D[1::2]

    diag_D = zeros(N_band,1);
    % Python diag_D[::2] -> MATLAB 1:2:end
    % rho_k_sum[1::2,1::2] -> MATLAB (2:2:end, 2:2:end)
    diag_D(1:2:end) = real(trace(rho_k_sum(2:2:end, 2:2:end)));
    % Python diag_D[1::2] -> MATLAB 2:2:end
    % rho_k_sum[::2,::2] -> MATLAB (1:2:end, 1:2:end)
    diag_D(2:2:end) = real(trace(rho_k_sum(1:2:end, 1:2:end)));

    diag_S = zeros(N_band,1);
    % Python diag_S[1::2] = diag_D[::2]
    %   -> MATLAB diag_S(2:2:end) = diag_D(1:2:end)
    diag_S(2:2:end) = diag_D(1:2:end);
    % Python diag_S[::2]  = diag_D[1::2]
    %   -> MATLAB diag_S(1:2:end) = diag_D(2:2:end)
    diag_S(1:2:end) = diag_D(2:2:end);

    % V_D_0 和 V_S_0 (q->0)
    V_D0 = 2*pi*ke/er * ( dgate - (N-1)*dlay + ((N-1)*dlay)^2/(4*dgate) );
    V_S0 = 2*pi*ke/er * ( dgate                - ((N-1)*dlay)^2/(4*dgate) );

    V_H = ( V_D0 * diag(diag_D) + V_S0 * diag(diag_S) ) / A;
    % 注释里提到：常数 V_S0*diag(diag_D+diag_S)/A 被减掉，
    % 并且 V_D0 - V_S0 < 0 对应 top/bottom 电子更均匀时能量更低。
    % 如果你想可以改成：
    % V_H = (V_D0 - V_S0) * diag(diag_D) / A;

    % -------------------- 5. 把 VI + Hartree 加到每个 k 上 --------------------
    for ix = 1:nkx
        for iy = 1:nky
            M = squeeze(V_F_k(ix,iy,:,:));  % [N_band x N_band]
            M = M + V_VI + V_H;
            V_F_k(ix,iy,:,:) = M;
        end
    end

    % -------------------- 6. flavor_symmetrize (暂时不做) --------------------
    % 如果后面需要 SU(4)/SU(2) 等对称约束，可以在这里加：
    % V_F_k = flavor_symmetrize_mat(V_F_k, hf.fix_sym, hf.SOC_dir);
end




%%
function SWM = get_SWM_params(N)
%GET_SWM_PARAMS  Slonczewski–Weiss–McClure 参数（单位 eV）
%
%   SWM = GET_SWM_PARAMS(N)
%
%   输入：
%     N : 层数 (2–6)
%
%   输出：
%     SWM.gamma0, gamma1, gamma2, gamma3, gamma4, delta  （单位 eV）
%
%   注：原 Python 代码中 gamma_* 单位为 meV，这里全部除以 1000 变成 eV

    switch N
        case 2
            SWM.gamma0 = 3160 / 1000;   % eV
            SWM.gamma1 =  500 / 1000;
            SWM.gamma2 =    0 / 1000;
            SWM.gamma3 = -280 / 1000;
            SWM.gamma4 = -200 / 1000;
            SWM.delta  = -0.05/1000;
        case 3
            SWM.gamma0 = 3160 / 1000;
            SWM.gamma1 =  460 / 1000;
            SWM.gamma2 =  -17.0 / 1000;
            SWM.gamma3 = -300  / 1000;
            SWM.gamma4 =  -86  / 1000;
            SWM.delta  =  -1.1/1000;
        case 4
            SWM.gamma0 = 3160 / 1000;
            SWM.gamma1 =  445 / 1000;
            SWM.gamma2 =  -18.2 / 1000;
            SWM.gamma3 = -319  / 1000;
            SWM.gamma4 =  -79  / 1000;
            SWM.delta  =  -0.066;
        case 5
            SWM.gamma0 = 3160 / 1000;
            SWM.gamma1 =  435 / 1000;
            SWM.gamma2 =  -18.5 / 1000;
            SWM.gamma3 = -322  / 1000;
            SWM.gamma4 =  -67.5 / 1000;
            SWM.delta  =  -0.147/1000;
        case 6
            SWM.gamma0 = 3160 / 1000;
            SWM.gamma1 =  430 / 1000;
            SWM.gamma2 =  -18.5 / 1000;
            SWM.gamma3 = -325  / 1000;
            SWM.gamma4 =  -73   / 1000;
            SWM.delta  =  -1.6/1000;
        case 15
            SWM.gamma0 = 3160 / 1000;
            SWM.gamma1 =  435 / 1000;
            SWM.gamma2 =  -18.5 / 1000;
            SWM.gamma3 = -322  / 1000;
            SWM.gamma4 =  -67.5   / 1000;
            SWM.delta  =  -0.2/1000;
        otherwise
            error('Number of layers is only supported between 2–6.');
    end
end

function H = H2B_5LG(kx, ky, q, SWM, N)
%H2B_5LG  N-layer graphene 单自旋、单 valley 的 2×2 有效模型
%
%   H = H2B_5LG(kx, ky, q, SWM, N)
%
%   kx, ky : 无量纲动量（k_phys * a_lat）
%   q      : valley 指示 (+1 = K, -1 = K')
%   SWM    : 结构体, get_SWM_params(N) 返回
%   N      : 层数（这里 N=5）

    % ---- 1. SWM -> hv0a, hv3a, hv4a 等 ----
    hv0a = sqrt(3)/2 * SWM.gamma0;   % eV
    hv3a = sqrt(3)/2 * SWM.gamma3;   % eV
    hv4a = sqrt(3)/2 * SWM.gamma4;   % eV
    gm1  = SWM.gamma1;
    gm2  = SWM.gamma2;
    delta= SWM.delta;

    % ---- 2. 复动量 Pi = q kx + i ky ----
    Pi   = q * kx + 1i * ky;      % 复数
    H    = complex(zeros(2,2));   % 确保是复矩阵

    % ---- 3. r_N, r_{N-1} 与 v3 线性修正 ----
    x2   = abs(hv0a * Pi / gm1).^2;   % |hv0a * Pi / gm1|^2

    % 处理 x2 = 1 的极限 (避免 0/0)，简单起见这里直接用公式，
    % 真要非常精确可以加一个 eps 判断。
    r_N   = (x2.^N     - 1) ./ (x2 - 1);
    r_Nm1 = (x2.^(N-1) - 1) ./ (x2 - 1);

    % Pi^3 + conj(Pi)^3
    Pi3_plus = Pi.^3 + conj(Pi).^3;

    % v3 的线性修正 (Pentalayer only!!!!!)
    r_N = r_N - hv0a^2 * hv3a / gm1^3 .* real(Pi3_plus) .* ...
                (1 + 2 * x2 + 3 * x2.^2);

    % ---- 4. H_ch ----
    H(1,2) = H(1,2) + (-gm1) * (hv0a * conj(Pi) / (-gm1)).^N;
    H(2,1) = H(2,1) + (-gm1) * (hv0a * Pi        / (-gm1)).^N;

    % ---- 5. H_s ----
    H(1,1) = H(1,1) + delta - 2 * abs(Pi).^2 .* hv0a .* hv4a ./ gm1 .* r_Nm1;
    H(2,2) = H(2,2) + delta - 2 * abs(Pi).^2 .* hv0a .* hv4a ./ gm1 .* r_Nm1;

    % ---- 6. H_tr ----
    pref_tr = ( (N - 2) * gm2 / 2 - (N - 1) * hv0a * hv3a .* abs(Pi).^2 ./ gm1 );

    H(1,2) = H(1,2) + pref_tr .* (hv0a * conj(Pi) / (-gm1)).^(N - 3);
    H(2,1) = H(2,1) + pref_tr .* (hv0a * Pi        / (-gm1)).^(N - 3);

    % ---- 7. v3 的二阶修正（只在 pentalayer 中保留）----
    pref_v3_2 = ( 3 * hv0a * hv3a^2 .* abs(Pi).^2 ./ gm1.^2 - hv3a * gm2 / gm1 );

    H(1,2) = H(1,2) + pref_v3_2 .* Pi;
    H(2,1) = H(2,1) + pref_v3_2 .* conj(Pi);

    % 对角修正
    pref_diag = hv0a * hv3a * hv4a ./ gm1.^2 .* real(Pi3_plus) .* ...
                (1 + 3 * x2 + 5 * x2.^2);

    H(1,1) = H(1,1) + pref_diag;
    H(2,2) = H(2,2) + pref_diag;

    % ---- 8. 整体除以 r_N ----
    H = H ./ r_N;
end

function Hk = H0_5LG_eff_cart(kx, ky, SWM, pars)
%H0_5LG_EFF_CART  8×8 有效哈密顿量 (N-layer graphene)
%
%   Hk = H0_5LG_eff_cart(kx, ky, SWM, pars)
%
%   kx, ky : 无量纲动量（k_phys * a_lat）
%   SWM    : get_SWM_params(N)
%   pars   : 结构体
%            .N           层数
%            .SOC         Ising SOC 强度 (eV)
%            .SOC_dir     0 -> Sz, 1 -> Sx
%            .SOC_single  1 -> 单侧 TMD, 0 -> 双侧 TMD
%
%   输出：
%     Hk [8×8]，基底为 spin(↑,↓) ⊗ valley(K,K') ⊗ orbital(1A, NB)

    N          = pars.N;
    SOC        = pars.SOC;
    SOC_dir    = pars.SOC_dir;
    SOC_single = pars.SOC_single;

    % ---- Pauli & 张量结构 ----
    sig0 = eye(2);
    sigx = [0 1; 1 0];
    sigz = [1 0; 0 -1];

    eye8 = kron(sig0, kron(sig0, sig0));   % 8×8 单位
    Sx   = kron(sigx, kron(sig0, sig0));   % spin Pauli
    Sz   = kron(sigz, kron(sig0, sig0));

    Vz   = kron(sig0, kron(sigz, sig0));   % valley Pauli
    Lz   = kron(sig0, kron(sig0, sigz));   % orbital Pauli

    % ---- 两个 valley 的 2×2 有效模型 ----
    H2_K  = H2B_5LG(kx, ky, +1, SWM, N);   % K valley
    H2_Kp = H2B_5LG(kx, ky, -1, SWM, N);   % K' valley

    % spin ⊗ valley(2band)
    H_spin_valley = kron(eye(2), blkdiag(H2_K, H2_Kp));  % 8×8

    % ---- Ising SOC ----
    Ising_SOC = SOC * Vz * (SOC_dir * Sx + (1 - SOC_dir) * Sz) * ...
               ((1 - SOC_single) * Lz + SOC_single * (eye8 - Lz) / 2);

    % ---- 总哈密顿量 ----
    Hk = H_spin_valley + Ising_SOC;
end

function A = get_sample_area(a_lat, dk)
%GET_SAMPLE_AREA  给定无量纲动量步长 dk, 返回实空间面积 A (长度单位^2)
%
%   A = (2π a_lat / dk)^2

    A = (2*pi * a_lat / dk)^2;
end
