b1=g.b(1,1:2);
b2=g.b(2,1:2);
% b2=[2.9494        , 0];
% b1=[-1.4747 ,   2.5543];
Ncoarse=20;
Ndense=100;
R_dense=0.1;
[k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense)
% [k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense);
%%
k_cart=reshape(k_cart)
%%
figure()
% plot(k_cart(:,1),k_cart(:,2),'o')
plot(k_frac(:,1),k_frac(:,2),'o')
%%
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

    % -------- 6. 标记 region_id: coarse(1), K_patch(2), Kp_patch(3) --------
    % 我们再次用 frac 空间的窗口来检查
    Nk = size(kfrac, 1);
    region_id = zeros(Nk, 1);

    % coarse region: 先设为 1
    region_id(:) = 1;

    % K patch
    deltaK_all  = abs(kfrac - K_frac);
    inK_all     = (deltaK_all(:,1) <= R_dense) & (deltaK_all(:,2) <= R_dense);

    % K' patch
    deltaKp_all = abs(kfrac - Kp_frac);
    inKp_all    = (deltaKp_all(:,1) <= R_dense) & (deltaKp_all(:,2) <= R_dense);

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
