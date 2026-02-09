clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the geomtery information                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = MTB.geometry("NGra");

a0=2.46;
a1=[sqrt(3)/2, -1/2, 0.0]*a0;
a2=[sqrt(3)/2, 1/2, 0.0]*a0;
a3=[0.0, 0.0, 62];
omega=dot(a1,cross(a2,a3));
% Calculate the reciprocal lattice vectors
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
b=[b1;b2;b3];
g.a = [a1;a2;a3];
g.b = b;
d0=3.35;
d0=d0/62;

N = 7; % Define the number of iterations for the loop
g.atoms = []; % Initialize the atoms property of g

for i=1:N
    % g.atoms=[g.atoms;...
    %     0.0+(i-1)/3,0.0+(i-1)/3,(i-1)*d0;...
    %     1.0/3+(i-1)/3,1.0/3+(i-1)/3,(i-1)*d0];
        g.atoms=[g.atoms;...
        mod(0.001+(i-1)/3,1),mod(0.001+(i-1)/3,1),(i-1)*d0;...
        mod(1.003/3+(i-1)/3,1),mod(1.003/3+(i-1)/3,1),(i-1)*d0];
        %         g.atoms=[g.atoms;...
        % mod(0.00+(i-1)/3,1),mod(0.00+(i-1)/3,1),(i-1)*d0;...
        % mod(1.00/3+(i-1)/3,1),mod(1.00/3+(i-1)/3,1),(i-1)*d0];
end
g.wpos=g.atoms*g.a;

g.wpos(:,3)=g.wpos(:,3)-mean(g.wpos(:,3))+a3(:,3)/2;
g.atoms=g.wpos*inv(g.a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the real space Ham                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=22;
% result_matrices = cell(1,n);
hopr=[]
% for i=1:n
%     result_matrices{i}=find_neighbor_data(g.wpos(1:end,:),g.a, i+1,3);
%     hopr=[hopr;result_matrices{i}(:,3:5)];
% end
n=4;
result_matrices = cell(1,4);
result_matrices{1}=find_neighbor_data(g.wpos(1:end,:),g.a, 2,3); %nn vf
hopr=[hopr;result_matrices{1}(:,3:5)];
result_matrices{2} = find_neighbor_data(g.wpos(1:end,:), g.a, 5, 3); % t1
hopr=[hopr;result_matrices{2}(:,3:5)];
result_matrices{3} = find_neighbor_data(g.wpos(1:end,:), g.a, 6, 3); % v3,v4
hopr=[hopr;result_matrices{3}(:,3:5)];
result_matrices{4} = find_neighbor_data(g.wpos(1:end,:), g.a, 21, 3); % t2
hopr=[hopr;result_matrices{4}(:,3:5)];

% result_matrices{4} = find_neighbor_data(g.wpos(1:end,:), g.a, 6, 3); % t4
%
%
hopr=unique(hopr,'rows');
nbands = size(g.wpos,1);
ham = zeros(nbands,nbands,size(hopr,1));

for n_idx=1:3
    result_matrix=result_matrices{n_idx};
for i =1:size(result_matrix,1)
    [~,raw_index]=ismember(result_matrix(i,3:5),hopr,"rows");
    orbital_1 = result_matrix(i,1);
    orbital_2 = result_matrix(i,2);
    % ham(orbital_1*2-1:orbital_1*2,orbital_2*2-1:orbital_2*2,raw_index) = [1,0;0,1]*t1;
    result_matrix(i,3:5)
    R=result_matrix(i,3:5)*g.a+g.wpos(orbital_2,:)-g.wpos(orbital_1,:);
    tsk=sk_appro_bab(R);
    % tsk=sk_appro(R);
    if ~isnan(tsk)
        % ham(orbital_1*2-1:orbital_1*2, orbital_2*2-1:orbital_2*2, raw_index) = ham(orbital_1*2-1:orbital_1*2, orbital_2*2-1:orbital_2*2, raw_index)+...
            % [1, 0; 0, 1] * tsk;
        ham(orbital_1, orbital_2, raw_index) = ham(orbital_1, orbital_2, raw_index)+...
            tsk;
    end
end
end
%
% result_matrices{4} = find_neighbor_data(g.wpos(1:end,:), g.a, 21, 3); % t2

result_matrix=result_matrices{4};
for i =1:size(result_matrix,1)
    [~,raw_index]=ismember(result_matrix(i,3:5),hopr,"rows");
    orbital_1 = result_matrix(i,1);
    orbital_2 = result_matrix(i,2);
    t2=-0.007;
    % ham(orbital_1*2-1:orbital_1*2, orbital_2*2-1:orbital_2*2, raw_index) = ham(orbital_1*2-1:orbital_1*2, orbital_2*2-1:orbital_2*2, raw_index)+...
    %     [1, 0; 0, 1] * t2;
    ham(orbital_1, orbital_2, raw_index) = ham(orbital_1, orbital_2, raw_index)+...
        t2;
end
g.hopr = hopr;
g.ham = ham;
g.wpos=kron(g.wpos,ones(1,1));
% Visp=-diag(abs(g.wpos(:,3)-d0/2)*0.005);
% % Visp=-diag(abs(g.wpos(:,3)-d0/2)*0.000001);
% [~,raw_index]=ismember([0,0,0],hopr,"rows");
% g.ham(:,:,raw_index)=g.ham(:,:,raw_index)+Visp;




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[1/3,2/3,0.0]*0.9,...
          [1/3,2/3,0.0],...
          [1/3,2/3,0.0]+([0.5,0.5,0.0]-[1/3,2/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;
efermi=0.73;
efermi=-0.01545;
efermi=-0.0
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
g.iniham=g.ham+0;
ylim([-0.1,0.1])
%%

% bands_file = 'data/Graphene/dft-bands/dft-bands-3s.dat';
bands_file = 'data/Graphene/dft-bands/dft-bands-3s-small.dat';
bands_data = readmatrix(bands_file);  % [nq × nmode]
hold on;
plot(bands_data(:,1),bands_data(:,2)+2.8077,'r.',LineWidth=2)
ylim([-1,1])
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi Surface-1             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate plane bands
%%
% clc;
% clear;
knum=801;
kxline=[-1,1];
kyline=[-1,1];
% x=1/3;
% y=2/3;
% kxline=[-0.01+x,0.01+x];
% kyline=[-0.01+y,0.01+y];

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
% surf(kx, ky, Ez2-Ez1, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
surf(kx, ky, Ez2-Ez1);
% colormap blue;        % 色图
colormap(slanCM('RdBu'))
colorbar;
shading interp
% clim([0,0.001])
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
% Ef_list=[2.0,2.0];
% 2. 费米面等高线（Ef）
% contour3(kx, ky, Ez1+2, Ef_list, 'LineColor', 'blue', 'LineWidth', 3);
% contour3(kx, ky, Ez2+2, Ef_list, 'LineColor', 'blue' ,'LineWidth', 1.5);

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

function tsk=sk_appro_bab(R)
V_pi = -2.810;
% V_pi = -2.75;
V_sigam = 0.48;
l0=3.364;
q_pi=3.1451;
a_pi=1.418;
q_sigma=7.428;
a_sigma=3.349;
r_c=6.14;
l_c=0.265;
z=[0,0,1]*R';
r=norm(R);
% Amp=V_pi*exp((acc-d)/delta)*(1-(R(3)/d)^2)+V_sigam*exp((d0-d)/delta)*(R(3)/d)^2;
tsk=V_pi*(1-z^2/r^2)*exp(q_pi*(1-r/a_pi))/(1+exp(r-r_c)/l_c)+...
    V_sigam*z^2/r^2*exp(q_sigma*(1-r/a_sigma))/(1+exp(r-r_c)/l_c);
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