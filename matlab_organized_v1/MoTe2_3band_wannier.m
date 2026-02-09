clear;
clear all;
%p=parpool(8)
g = MTB.geometry("3band_MoTe2");
g = MTB.read_poscar(g,"data/MoTe2/3band_wannier/hr_file/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/MoTe2/3band_wannier/hr_file/wannier90_hr_p1.dat','data/MoTe2/3band_wannier/hr_file/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [2/3,1/3,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("3band_MoTe2");
g = MTB.read_poscar(g,"data/MoTe2/3band_wannier/hr_file/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/MoTe2/3band_wannier/hr_file/wannier90_hr_p1.dat','data/MoTe2/3band_wannier/hr_file/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;

result_matrices = find_neighbor_data(g.wpos(1:1:end, :), g.a, 2, 3);
%
gpair=MTB.geometry("Pair");
gpair.a=g.a;
gpair.b=g.b;
gpair.wpos=g.wpos;
gpair.atoms=g.atoms;
gpair.hopr=unique(result_matrices(:,3:5),'rows');
gpair.ham=zeros(size(g.ham,1),size(g.ham,1),size(gpair.hopr,1));
vector_21=(gpair.atoms(2,:)-gpair.atoms(1,:))*gpair.a;
delta=0;
phase=zeros(1,length(result_matrices));
for row_idx=1:length(result_matrices)
    orbital_idx_start=result_matrices(row_idx,1);
    orbital_idx_end=result_matrices(row_idx,2);
    vector_delta=gpair.atoms(orbital_idx_end,:)-gpair.atoms(orbital_idx_start,:);
    vector_R=result_matrices(row_idx,3:5);
    vector=vector_R+vector_delta;
    vector_card=vector*gpair.a;
    phase_tmp=atan2(vector_card(2),vector_card(1));
    hopr_idx=find(ismember(gpair.hopr,result_matrices(row_idx,3:5),'rows'));
    phase(row_idx)=phase_tmp;
    gpair.ham(orbital_idx_end,orbital_idx_start,hopr_idx)=delta*exp(1j*(phase_tmp));
end

mu=0;
[nbands,~,nrpts]=size(g.ham);
labels={'-K','\Gamma','K'}; % labels for k
hkpoints={[-0.333333,-0.333333,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG_pwave(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b,mu,gpair);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.05,0.05])
%%
%
nslab=20;
gs=MTB.ham.get_supercell(g,nslab,1);
row_idx=find(~ismember(gs.hopr(:,1),0,'rows'));

gs.ham(:,:,row_idx)=[];
gs.hopr(row_idx,:)=[];
gspair=MTB.ham.get_supercell(gpair,nslab,1);
gspair.ham(:,:,[1,2,3,7,8,9])=[];
gspair.hopr([1,2,3,7,8,9],:)=[];

mu=0.0;
[nbands,~,nrpts]=size(gs.ham);
labels={'-M','\Gamma','M'}; % labels for k
hkpoints={[-0.0,-0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG_pwave(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b,mu,gspair);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.05,0.05])
%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("3band_MoTe2");
g = MTB.read_poscar(g,"data/MoTe2/3band_wannier/hr_file/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/MoTe2/3band_wannier/hr_file/wannier90_hr_p1.dat','data/MoTe2/3band_wannier/hr_file/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [2/3,1/3,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
file = 'data/MoTe2/G0/delta_R.h5';
info = h5info(file);

% 数据集数量
nDataset = numel(info.Datasets);

% 用 cell 数组存储数据
data = cell(nDataset,1);
names = cell(nDataset,1);

% 循环读取
for i = 1:nDataset
    names{i} = info.Datasets(i).Name;                % 数据集名字
    data{i} = h5read(file, ['/' names{i}]);          % 读取数据
end
%%
clc;
clear;
g = MTB.geometry("3band_MoTe2");
g = MTB.read_poscar(g,"data/MoTe2/3band_wannier/hr_file/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/MoTe2/3band_wannier/hr_file/wannier90_hr_p1.dat','data/MoTe2/3band_wannier/hr_file/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;

result_matrices = find_neighbor_data(g.wpos(1:1:end, :), g.a, 2, 3);
%
gpair=MTB.geometry("Pair");
gpair.a=g.a;
gpair.b=g.b;
gpair.wpos=g.wpos;
gpair.atoms=g.atoms;
gpair.hopr=g.hopr;
gpair.ham=zeros(size(g.ham));
%
file = 'data/MoTe2/G0/delta_R.h5';
info = h5info(file);
% 数据集数量
nDataset = numel(info.Datasets);
% 用 cell 数组存储数据
data = cell(nDataset,1);
names = cell(nDataset,1);
% 循环读取
for i = 1:nDataset
    names{i} = info.Datasets(i).Name;                % 数据集名字
    data{i} = h5read(file, ['/' names{i}]);          % 读取数据
end

for r_idx = 1:nDataset
    v = gpair.hopr(r_idx,:);
    str_v = sprintf('(%d, %d, %d)', v);
    data_tem = h5read(file, ['/' str_v]); 
    ham=data_tem.r+data_tem.i*1j;
    gpair.ham(:,:,r_idx)=ham;
end

mu=0.0;
[nbands,~,nrpts]=size(g.ham);
labels={'-K','\Gamma','K'}; % labels for k
hkpoints={[-0.333333,-0.333333,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG_pwave(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b,mu,gpair);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.05,0.05])
%%
nslab=20;
gs=MTB.ham.get_supercell(g,nslab,1);
row_idx=find(~ismember(gs.hopr(:,1),0,'rows'));

gs.ham(:,:,row_idx)=[];
gs.hopr(row_idx,:)=[];
gspair=MTB.ham.get_supercell(gpair,nslab,1);
row_idx=find(~ismember(gspair.hopr(:,1),0,'rows'));
gspair.ham(:,:,row_idx)=[];
gspair.hopr(row_idx,:)=[];


mu=0.00;
[nbands,~,nrpts]=size(gs.ham);
labels={'-M','\Gamma','M'}; % labels for k
hkpoints={[-0.0,-0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG_pwave(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b,mu,gspair);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.05,0.05])
%%

%%


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
    rounded_distances = round(raw_distances, 2); % 保留两位小数，分组距离
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
