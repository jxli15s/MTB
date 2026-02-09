%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the geomtery information                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
N=3;
a0=2.46;
c = 62;
d0=3.36;
% d0=3.40;
g = build_NGra(N, a0, c, d0);
V_pi = -2.810;
V_sigam = 0.48;
l0=3.364; q_pi=3.1451; a_pi=1.418; q_sigma=7.428;
a_sigma=3.349; r_c=6.14; l_c=0.265;
parameters_tsk=[V_pi,V_sigam,l0,q_pi,a_pi,q_sigma,a_sigma,r_c,l_c];
shells=[2,5,6,21];
t2=-0.007;
visp_val=5e-5;
[g, hopr, result_matrices] = get_ham_SK(g, parameters_tsk, shells, t2, visp_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.9,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;
efermi=-0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
g.iniham=g.ham+0;
ylim([-0.1,0.1])
% bands_file = 'data/Graphene/dft-bands/dft-bands-3s.dat';
bands_file = 'data/Graphene/dft-bands/dft-bands-3s-small.dat';
bands_data = readmatrix(bands_file);  % [nq × nmode]
bands_data(:,2)=bands_data(:,2)+2.8066;
kpath_vasp=bands_data(:,1);
bands_vasp=bands_data(:,2:end);
kpath_vasp=reshape(kpath_vasp,[],96);
bands_vasp=reshape(bands_vasp,[],96);
bands_vasp=bands_vasp(:,10:15);
kpath_vasp(end/2,:)=[];
bands_vasp(end/2,:)=[];
hold on;
% plot(bands_data(:,1),bands_data(:,2),'r.',LineWidth=2)
for i=1:size(bands_vasp,2)
    plot(kpath_vasp(:,i),bands_vasp(:,i),'r-',LineWidth=2)
end
ylim([-1,1])
%%
p=[-2.810,0.48];
E_DFT=bands_vasp;
%% ====================== 1. 拟合设置 ======================
% 初值、边界（如无先验可放宽）
V_pi = -2.810;
V_sigam = 0.48;
l0=3.364; q_pi=3.1451; a_pi=1.418; q_sigma=7.428;
a_sigma=3.349; r_c=6.14; l_c=0.265;

params0 = [-2.810, 0.48,-0.007,1e-5];              % 初始猜测
lb = [-4, 0.3, -0.02,0];                    % 下界
ub = [-2, 0.6,   0.02,1e-2];                    % 上界

% —— 构造"权重矩阵" weights: Nk×Nb —— 
% 你可以按能量窗口、k 点、或某些 band 定制
weights = ones(size(E_DFT));

% (a) 费米能附近加权：|E- Ef| < 0.1 eV 权重*10
win = 0.015;
mask_E = abs(E_DFT - 0) < win;
weights(mask_E) = 50;

% (b) 高对称点附近加权（这里以 k=0 为例，±5% 区间放大）
% % k_center = 0;
% % mask_k = abs(klist - k_center) < 0.05 * (kmax - kmin);
% % weights(mask_k, :) = weights(mask_k, :) * 2;

% (c) band-selective：例如第一条带更重要
% % band_weight = [2.0, 1.0];   % 对每条带的全局权重
% % weights = weights .* band_weight;  % 自动按列广播到 Nk×Nb

% —— 残差函数（加权）——
% 注意：lsqnonlin 期望"向量残差"，我们把 Nk×Nb 展平成一维向量
residual_fun = @(p) reshape( weights .* (get_sk_bands(p) - E_DFT), [], 1 );

% 选项
opts = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5e4, ...
    'MaxIterations', 2e3, ...
    'StepTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12);

%% ====================== 2. 运行拟合 ======================
[params_fit, ~, ~, exitflag, output] = lsqnonlin(residual_fun, params0, lb, ub, opts); %#ok<ASGLU>
disp('===== 拟合完成 =====');
disp(['params_fit = [t, Delta] = [', num2str(params_fit(1),'%.6f'), ', ', num2str(params_fit(2),'%.6f'), ']']);
disp(['exitflag = ', num2str(exitflag)]);
disp(output.message);
%
bands_sk=get_sk_bands(params_fit);
% plot(bands_data(:,1),bands_data(:,2),'r.',LineWidth=2)
figure()
hold on;
for i=1:size(bands_vasp,2)
    plot(kpath_vasp(:,i),bands_vasp(:,i),'b-',LineWidth=2)
    plot(kpath_vasp(:,i),bands_sk(:,i),'r--',LineWidth=2)
end
ylim([-1,1])





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi Surface-1             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate plane bands
%%
% clc;
% clear;
knum=401;
kxline=[0,0.5];
kyline=[-0.5,0];
% kxline=[0,0.5];
% kyline=[0,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Write the plane eigenvalue               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Occ=2;
% filename='NbSe2_1sEnk-501x501.dat';
% writeEnk(Enk,Kx,Ky,Occ,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx_list, ky_list: 分别是 nkx × 1 和 nky × 1 的向量
% load('E_data.mat'); % 或你已经在工作区中

band_index = nbands/2;  % 要绘制的能带索引
% Ef = 0;          % 费米能级（假设为0）

% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序
kx=Kx;
ky=Ky;
% 获取目标能带的能量值（转置以匹配 meshgrid）
Ez = squeeze(Enk(:,:,band_index))';          

% 绘图
figure;
hold on;

% 1. 三维能带
% s = surf(kx, ky, Ez, 'EdgeColor', 'none');
% colormap turbo
shading interp
alpha(0.9);
Ez1 = squeeze(Enk(:,:,band_index))';
% surf(kx, ky, Ez1, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
Ez2 = squeeze(Enk(:,:,band_index+1))';
surf(kx, ky, Ez2-Ez1, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
% colormap blue;        % 色图
colormap(slanCM('RdBu'))
colorbar;
axis equal tight;
% % lighting phong    % 或 gouraud，phong 更平滑
% shading interp    % 表面平滑
% for band = 1:2
%     Ez = squeeze(Enk(:,:,band))';
%     surf(kx, ky, Ez, 'EdgeColor', 'none', 'FaceAlpha', 0.6,'FaceColor', [1,0,0]);
% end



% Ef=0.0;
% 2. 费米面等高线（Ef）
% contour3(kx, ky, Ez, [Ef Ef], 'k', 'LineWidth', 2);  % Fermi contour

% Ef_list = -0.2:0.2:0.6;  % 多个等高值
% Ef_list =[-0.2 -0.1]
Ef_list=[2.0,2.0];
% 2. 费米面等高线（Ef）
contour3(kx, ky, Ez1+2, Ef_list, 'LineColor', 'blue', 'LineWidth', 3);
contour3(kx, ky, Ez2+2, Ef_list, 'LineColor', 'blue' ,'LineWidth', 1.5);

% 3. 视图与标签
view(3);
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
title(['Bands and Fermi Surface ']);
colorbar;
% axis tight;
axis equal;
box on
view(0,90)
%%
function Amp=sk_appro(R)
V_pi = -2.7;
V_sigam = 0.48;
delta = 0.045; 
d=norm(R)/10;
acc=0.142;
d0=0.335;
Amp=V_pi*exp((acc-d)/delta)*(1-(R(3)/d)^2)+V_sigam*exp((d0-d)/delta)*(R(3)/d)^2;
end

function tsk=sk_appro_bab(parameters_tsk,R)
V_pi = parameters_tsk(1);
V_sigam = parameters_tsk(2);
l0=parameters_tsk(3);
q_pi=parameters_tsk(4);
a_pi=parameters_tsk(5);
q_sigma=parameters_tsk(6);
a_sigma=parameters_tsk(7);
r_c=parameters_tsk(8);
l_c=parameters_tsk(9);
z=[0,0,1]*R';
r=norm(R);
tsk=V_pi*(1-z^2/r^2)*exp(q_pi*(1-r/a_pi))/(1+exp(r-r_c)/l_c)+...
    V_sigam*z^2/r^2*exp(q_sigma*(1-r/a_sigma))/(1+exp(r-r_c)/l_c);
end


%%

function bands_sk=get_sk_bands(p)
% V_pi = -2.810;
% V_sigam = 0.48;
V_pi=p(1);
V_sigam=p(2);
l0=3.364; q_pi=3.1451; a_pi=1.418; q_sigma=7.428;
a_sigma=3.349; r_c=6.14; l_c=0.265;
parameters_tsk=[V_pi,V_sigam,l0,q_pi,a_pi,q_sigma,a_sigma,r_c,l_c];
shells=[2,5,6,21];
t2=p(3);
visp_val=p(4);

N=3;
a0=2.46;
c = 62;
d0=3.35;

g = build_NGra(N, a0, c, d0);
[g, ~, ~] = get_ham_SK(g, parameters_tsk, shells, t2, visp_val);
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.9,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=501;
efermi=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);

bands_sk=Energy';
% MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
% % figure()
% % hold on;
% % for i=1:size(bands_vasp,2)
% %     plot(kpath_vasp(:,i),bands_sk(:,i),'k-',LineWidth=2)
% %     plot(kpath_vasp(:,i),bands_vasp(:,i),'r--',LineWidth=2)
% % end
ylim([-1,1])
end

function g = build_NGra(N, a0, c, d0)
%BUILD_NGRA  构造 N 层（ABC 位移）石墨烯几何（MTB.geometry）
%
%   g = build_NGra(N, a0, c, d0)
%
% 参数（均有默认值，可留空）：
%   N   : 层数，默认 3
%   a0  : 面内基矢长度参数，默认 2.46 (Å)
%   c   : 垂直方向晶格常数（含真空厚度），默认 62 (Å)
%   d0  : 相邻层间距（Å），默认 3.35 (Å)
%
% 说明：
% - 与你原始脚本完全一致：两原子/层，下一层在 a1/a2 方向各平移 1/3（即 ABC 堆叠式位移），
%   z 方向按 d0 间距放置，然后整体沿 z 居中到 a3(3)/2。
% - 返回值 g 为 MTB.geometry("NGra") 对象，并写入 g.a/g.b/g.atoms/g.wpos。
%
% 示例：
%   g = build_NGra(5);                       % 5 层，默认 a0=2.48, c=62, d0=3.35
%   g = build_NGra(8, 2.46, 60, 3.35);       % 自定义参数

    if nargin < 1 || isempty(N),  N = 3;     end
    if nargin < 2 || isempty(a0), a0 = 2.46; end
    if nargin < 3 || isempty(c),  c  = 62.0; end
    if nargin < 4 || isempty(d0), d0 = 3.35; end

    % 几何对象
    g = MTB.geometry("NGra");

    % 直接照你脚本定义基矢
    a1 = [sqrt(3)/2, -1/2, 0.0] * a0;
    a2 = [sqrt(3)/2,  1/2, 0.0] * a0;
    a3 = [0.0, 0.0, c];

    g.a = [a1; a2; a3];

    % 体积与倒格矢
    omega = dot(a1, cross(a2, a3));
    b1 = 2*pi*cross(a2, a3) / omega;
    b2 = 2*pi*cross(a3, a1) / omega;
    b3 = 2*pi*cross(a1, a2) / omega;
    g.b = [b1; b2; b3];

    % 层间距转为分数坐标（相对于 c）
    d0_frac = d0 / c;

    % 原子分数坐标（两原子/层），保持你的位移方案 (i-1)/3
    g.atoms = zeros(2*N, 3);
    for i = 1:N
        zf = (i-1) * d0_frac;   % 分数坐标的 z
        s  = (i-1) / 3;         % 每层在 a1/a2 方向的 1/3 平移
        % g.atoms(2*i-1, :) = [0.0   + s, 0.0   + s, zf];
        % g.atoms(2*i,   :) = [1.0/3 + s, 1.0/3 + s, zf];

        % 如果需要避免正好落在 0 或 1/3，可改用你注释里的"mod"写法：
        g.atoms(2*i-1,:) = [mod(0.001   + s, 1), mod(0.001   + s, 1), zf];
        g.atoms(2*i,  :) = [mod(1.003/3 + s, 1), mod(1.003/3 + s, 1), zf];
    end

    % 分数 -> 笛卡尔
    g.wpos = g.atoms * g.a;

    % 沿 z 居中到 a3(3)/2（与原脚本一致）
    g.wpos(:,3) = g.wpos(:,3) - mean(g.wpos(:,3)) + a3(3)/2;

    % 再转回分数坐标（保持你的 inv 写法，也可用右除 / 更稳健）
    g.atoms = g.wpos / g.a;  % 等价于 g.wpos*inv(g.a)
end

function [g, hopr, result_matrices] = get_ham_SK(g, parameters_tsk, shells, t2_value, z_slope)
%ASSEMBLE_HAM_SK  基于邻居搜索与Slater–Koster近似装配哈密顿量
%
%   [g, hopr, result_matrices] = get_ham_SK(g, parameters_tsk, shells, t2_value, z_slope)
%
% 输入：
%   g              : 含有 g.wpos (笛卡尔坐标) 与 g.a 的几何对象（如 MTB.geometry）
%   parameters_tsk : 传给 sk_appro_bab 的参数（由你外部定义）
%   shells         : 邻居"壳层"索引数组（用于 find_neighbor_data 的第3个参数）
%                    默认 [2 5 6 21]，其中前 3 个用 SK 计算，最后一个用常数 t2_value
%   t2_value       : 常数跃迁值（对应 shells(end)），默认 -0.007
%   z_slope        : 垂直势场系数（|z|*z_slope），默认 5e-5
%
% 输出：
%   g              : 增添了 g.hopr (R 向量列表) 与 g.ham (nbands×nbands×nR) 的结构
%   hopr           : 所有去重后的平移向量（分数坐标），每行 [R1 R2 R3]
%   result_matrices: 每个 shell 的邻居搜索结果（与原脚本一致）
%
% 依赖：
%   find_neighbor_data(wpos, a, shell_id, 3)
%   sk_appro_bab(parameters_tsk, R_cart)  % 返回标量跃迁 t
%
% 备注：
%   - 行为与原脚本等价：shells(1:3) 用 SK；shells(end) 加常数 t2_value。
%   - 自动确保 hopr 包含 [0 0 0]，以便叠加 on-site 势 V(z) = -|z|*z_slope。
%   - 使用 'unique(...,"rows","stable")' 保持首出现顺序稳定。

    if nargin < 3 || isempty(shells),   shells   = [2 5 6 21]; end
    if nargin < 4 || isempty(t2_value), t2_value = -0.007;     end
    if nargin < 5 || isempty(z_slope),  z_slope  = 5e-5;       end

    % 1) 为每个 shell 寻找邻居
    nshell = numel(shells);
    result_matrices = cell(1, nshell);
    hopr = [];
    for k = 1:nshell
        result_matrices{k} = find_neighbor_data(g.wpos(1:end, :), g.a, shells(k), 3);
        hopr = [hopr; result_matrices{k}(:, 3:5)]; %#ok<AGROW>
    end

    % 2) 去重（保持出现顺序稳定）
    hopr = unique(hopr, 'rows', 'stable');

    % 3) 预分配 Hamiltonian
    nbands = size(g.wpos, 1);
    ham = zeros(nbands, nbands, size(hopr, 1));

    % 4) 前 nshell-1 个 shell 用 SK 计算（与你原脚本的前三个一致）
    n_sk = max(0, nshell - 1);
    for n_idx = 1:n_sk
        result_matrix = result_matrices{n_idx};
        for i = 1:size(result_matrix, 1)
            % hopr 索引
            [~, raw_index] = ismember(result_matrix(i, 3:5), hopr, 'rows');

            orbital_1 = result_matrix(i, 1);
            orbital_2 = result_matrix(i, 2);

            % 位移矢量（笛卡尔）
            R = result_matrix(i, 3:5) * g.a + ...
                g.wpos(orbital_2, :) - g.wpos(orbital_1, :);

            % SK 近似的跃迁
            tsk = sk_appro_bab(parameters_tsk, R);

            if ~isnan(tsk)
                ham(orbital_1, orbital_2, raw_index) = ...
                    ham(orbital_1, orbital_2, raw_index) + tsk;
            end
        end
    end

    % 5) 最后一个 shell 用常数 t2_value（与你原脚本的第 4 个 shell 一致）
    if nshell >= 1
        result_matrix = result_matrices{end};
        for i = 1:size(result_matrix, 1)
            [~, raw_index] = ismember(result_matrix(i, 3:5), hopr, 'rows');
            orbital_1 = result_matrix(i, 1);
            orbital_2 = result_matrix(i, 2);
            ham(orbital_1, orbital_2, raw_index) = ...
                ham(orbital_1, orbital_2, raw_index) + t2_value;
        end
    end

    % 6) 写回 hopr/ham
    g.hopr = hopr;
    g.ham  = ham;

    % 7) 叠加 on-site 垂直势 V(z) = -|z| * z_slope 到 R=[0 0 0] 通道
    [tf0, raw0] = ismember([0 0 0], hopr, 'rows');
    if ~tf0
        % 若未出现 on-site 通道，则追加
        g.hopr(end+1, :) = [0 0 0];
        g.ham(:, :, end+1) = 0;
        raw0 = size(g.hopr, 1);
    end
    Visp = -diag(abs(g.wpos(:, 3)) * z_slope);
    g.ham(:, :, raw0) = g.ham(:, :, raw0) + Visp;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Function to find the n-th NN neighbor pairing          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result_matrix = find_neighbor_data(coords_cartesian, lattice_vectors, n, dimensionality)
    % 动态处理二维或三维晶格，计算第 n 近邻的原子对信息
    %
    % 输入:
    % coords_cartesian: 原胞内的笛卡尔坐标 (n_atoms x 3)
    % lattice_vectors: 晶格基矢量 (3 x 3)
    % n: 第 n 近邻
    % dimensionality: 2 表示二维晶格，3 表示三维晶格
    %
    % 输出:
    % result_matrix: (num_pairs x 6) 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 1. 将笛卡尔坐标转换为分数坐标
    inv_lattice = inv(lattice_vectors);
    coords_fractional = coords_cartesian * inv_lattice;

    % 2. 构造超胞
    if dimensionality == 2
        % 仅扩展 x 和 y
        [super_coords, super_indices] = construct_supercell_2d(coords_fractional, lattice_vectors, 1);
    elseif dimensionality == 3
        % 扩展 x, y, z
        [super_coords, super_indices] = construct_supercell(coords_fractional, lattice_vectors, 1);
    else
        error('Dimensionality must be 2 or 3');
    end

    % 3. 在 Non-PBC 条件下计算超胞的距离矩阵
    if dimensionality == 2
        % 仅考虑 x-y 平面距离
        dist_matrix_non_pbc = compute_distance_matrix_2d(super_coords, lattice_vectors, false);
    elseif dimensionality == 3
        % 考虑完整三维距离
        dist_matrix_non_pbc = compute_distance_matrix(super_coords, lattice_vectors, false);
    end

    % 4. 提取 unique 距离（非零，考虑浮点误差）
    raw_distances = triu(dist_matrix_non_pbc);
    % raw_distances = raw_distances(raw_distances > 0); % 去掉 0 距离
    rounded_distances = round(raw_distances, 3); % 保留两位小数，分组距离
    unique_distances = unique(rounded_distances, 'sorted');

    % 5. 找到第 n 近邻的距离
    if n > length(unique_distances)
        error('第 %d 近邻超出最大可能的距离范围', n);
    end
    nth_distance = unique_distances(n);

    % 6. 找到满足第 n 近邻距离的原子对
    [pair_i, pair_j] = find(abs(dist_matrix_non_pbc - nth_distance) < 1e-2); % 容忍浮点误差

    % 7. 直接生成结果矩阵
    result_matrix=find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors);
end

function unique_pairs = find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors)
    % 找到所有第 n 近邻的原子对并去除重复
    %
    % 输入:
    % pair_i, pair_j: 满足第 n 近邻条件的原子对索引
    % super_coords: 超胞中原子的分数坐标
    % super_indices: 超胞中原子的索引和周期性偏移
    % nth_distance: 第 n 近邻的距离
    % lattice_vectors: 晶格基矢量
    %
    % 输出:
    % unique_pairs: 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 初始化结果存储
    num_pairs = length(pair_i);
    all_pairs = zeros(num_pairs, 6); % [i, j, frax, fray, fraz, dis]

    % 遍历所有原子对
    for k = 1:num_pairs
        i = pair_i(k); % 原子 i
        j = pair_j(k); % 原子 j

        % 原胞内的原子编号
        atom_i = super_indices(i, 1);
        atom_j = super_indices(j, 1);

        % 计算周期性偏移（分数坐标差）
        delta_r = super_coords(i, :) - super_coords(j, :);

        % 将分数坐标差转换为笛卡尔坐标，用于计算实际距离
        delta_cartesian = delta_r * lattice_vectors;
        distance = sqrt(sum(delta_cartesian.^2));

        % 检查是否接近第 n 近邻距离
        if abs(distance - nth_distance) < 1e-1
            % 记录原子对信息
            all_pairs(k, :) = [atom_i, atom_j, super_indices(j,2:end)-super_indices(i,2:end), nth_distance];
        end
    end

    % 去除重复的原子对（如 i->j 和 j->i）
    % 按原子对的编号排序，并提取唯一值
    [~, unique_rows] = unique(all_pairs, 'rows');
    unique_pairs = all_pairs(unique_rows, :); % 提取唯一行
end

function dist_matrix = compute_distance_matrix(coords, lattice_vectors, pbc)
    % 计算三维距离矩阵，可选是否使用周期性边界条件 (PBC)
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            delta_cartesian = delta_r * lattice_vectors; % 转为笛卡尔坐标
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function dist_matrix = compute_distance_matrix_2d(coords, lattice_vectors, pbc)
    % 计算二维距离矩阵（仅考虑 x 和 y 方向），可选 PBC
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            % 仅保留 x 和 y 方向的笛卡尔坐标
            delta_cartesian = delta_r * lattice_vectors; 
            delta_cartesian = delta_cartesian(:, 1:2); % x 和 y 方向
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function [super_coords, super_indices] = construct_supercell(coords, lattice_vectors, scale)
    % 构造三维超胞
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            for k = -scale:scale
                offset = [i, j, k];
                super_coords = [super_coords; coords + offset];
                for a = 1:n_atoms
                    super_indices = [super_indices; a, i, j, k];
                end
            end
        end
    end
end

function [super_coords, super_indices] = construct_supercell_2d(coords, lattice_vectors, scale)
    % 构造二维超胞（仅扩展 x 和 y 方向）
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            offset = [i, j, 0]; % z 方向固定为 0
            super_coords = [super_coords; coords + offset];
            for a = 1:n_atoms
                super_indices = [super_indices; a, i, j, 0];
            end
        end
    end
end

function fit_band_weighted_demo()
%% ====================== 0. 数据准备 ======================
% 假设：klist 为 Nk×dim 的 k 点坐标（这里只做 1D 示例，dim=1）；
%       E_DFT 为 Nk×Nb 的 DFT 能带矩阵（Nb 条带）。
% 你可以把本节替换为实际的加载：
%   load my_dft.mat  % 里面有 klist, E_DFT, Ef
%
% 这里生成一份"伪 DFT 数据"来演示（真实使用时删除本节，用你的数据）:
Nk   = 201;
kmin = -pi; kmax = pi;
klist = linspace(kmin, kmax, Nk).';   % Nk×1
Nb   = 2;                              % 假设有两条带
Ef   = 0;                              % 费米能（如有请用你的）

% 设定"真实参数"（用来生成伪DFT）：params_true = [t, Delta]
params_true = [1.30, 0.12];  % 最近邻 hopping 和 子晶格势举例
E_clean = E_KS_example(params_true, klist);          % Nk×2
rng(1);
noise_level = 0.003;                                % 小噪声
E_DFT = E_clean + noise_level * randn(size(E_clean));

%% ====================== 1. 拟合设置 ======================
% 初值、边界（如无先验可放宽）
params0 = [1.0, 0.05];              % 初始猜测
lb = [0.0, 0.0];                    % 下界
ub = [3.0, 0.5];                    % 上界

% —— 构造"权重矩阵" weights: Nk×Nb —— 
% 你可以按能量窗口、k 点、或某些 band 定制
weights = ones(size(E_DFT));

% (a) 费米能附近加权：|E- Ef| < 0.1 eV 权重*10
win = 0.10;
mask_E = abs(E_DFT - Ef) < win;
weights(mask_E) = 10;

% (b) 高对称点附近加权（这里以 k=0 为例，±5% 区间放大）
k_center = 0;
mask_k = abs(klist - k_center) < 0.05 * (kmax - kmin);
weights(mask_k, :) = weights(mask_k, :) * 2;

% (c) band-selective：例如第一条带更重要
band_weight = [2.0, 1.0];   % 对每条带的全局权重
weights = weights .* band_weight;  % 自动按列广播到 Nk×Nb

% —— 残差函数（加权）——
% 注意：lsqnonlin 期望"向量残差"，我们把 Nk×Nb 展平成一维向量
residual_fun = @(p) reshape( weights .* (E_KS_example(p, klist) - E_DFT), [], 1 );

% 选项
opts = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5e4, ...
    'MaxIterations', 2e3, ...
    'StepTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12);

%% ====================== 2. 运行拟合 ======================
[params_fit, ~, ~, exitflag, output] = lsqnonlin(residual_fun, params0, lb, ub, opts); %#ok<ASGLU>
disp('===== 拟合完成 =====');
disp(['params_fit = [t, Delta] = [', num2str(params_fit(1),'%.6f'), ', ', num2str(params_fit(2),'%.6f'), ']']);
disp(['exitflag = ', num2str(exitflag)]);
disp(output.message);

%% ====================== 3. 结果可视化 ======================
E_fit = E_KS_example(params_fit, klist);
figure; hold on; box on; grid on;
% 画 DFT（点）与拟合（线）
plot(klist, E_DFT(:,1), '.', 'MarkerSize', 8); 
plot(klist, E_DFT(:,2), '.', 'MarkerSize', 8);
plot(klist, E_fit(:,1), '-', 'LineWidth', 1.8);
plot(klist, E_fit(:,2), '-', 'LineWidth', 1.8);
yline(Ef, '--'); 
xlabel('k'); ylabel('Energy (eV)');
legend({'DFT band 1','DFT band 2','Fit band 1','Fit band 2','E_F'}, 'Location','best');
title('加权最小二乘：KS 参数拟合 DFT 能带');

% 可选：显示权重热度（按 k 合并两条带的平均权重）
w_mean = mean(weights,2);
figure; plot(klist, w_mean, 'o-'); grid on; box on;
xlabel('k'); ylabel('mean weight'); title('权重分布（按 k 平均）');

end

%% ====================== 你的 KS/TB 能带函数 ======================
% 把这个函数替换为你的 KS/TB 模型即可。
% 输入:
%   params : 参数向量（例如最近邻 t 和子晶格势 Delta）
%   klist  : Nk×1（或 Nk×dim）k 点
% 输出:
%   E      : Nk×Nb 的能带矩阵
function E = E_KS_example(params, klist)
    % 一个极简的两带模型示例（1D Dirac-like + cos 带）
    % H(k) = [ Delta,   2t cos(k) ;
    %          2t cos(k),  -Delta ]
    % 特征值: E = ± sqrt(Delta^2 + (2t cos k)^2)
    t     = params(1);
    Delta = params(2);
    c = 2*t*cos(klist);
    Ek = sqrt(Delta.^2 + c.^2);
    E = [-Ek, +Ek];   % Nk×2

end