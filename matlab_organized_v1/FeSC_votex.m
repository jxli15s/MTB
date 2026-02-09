%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       Read FeSC Hamiltonian                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%p=parpool(8)
g = MTB.geometry("FeSC");
g = MTB.read_poscar(g,"data/FeSe_TI/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%%
% result_matrix = find_neighbor_data(g.atoms, g.a, 1, 2);
result_matrices = find_neighbor_data(g.wpos(1:2:end, :), g.a, 2, 3);
%% ===Calculate bulk bands======
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Z'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kindex,"FeSe",Electric_field_in_evpA*10000);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs = MTB.geometry("FeSC");
gs = MTB.read_poscar(gs,"data/FeSe_TI/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
gs.wpos=gs.atoms*gs.a;

MillerIndices=[0,0,1];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=50;

[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5],...
          [0.0000000000,0.0000000000],...
          [0.0,0.5]};% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("FeSC");
gs = MTB.read_poscar(gs,"data/FeSe_TI/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
gs.wpos=gs.atoms*gs.a;
MillerIndices=[0,0,1];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

[nbands,~,nrpts]=size(gs.ham);
hkpoints={[0.5,0.5000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
labels={'M','\Gamma','X'};
Np=2;
omegamin=-1;
omegamax=1;
omeganum=500;
omegas=linspace(omegamin,omegamax,omeganum);
nk=201;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);
%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('magma')); %magma plasma inferno cividis inferno hot heat

shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_r)
colormap(slanCM('ice'))
shading interp
caxis([1, 50])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_bulk)
colormap(slanCM('heat'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       Check open_xy model bands                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===create the tb model======
clc;
clear all;
g = MTB.geometry("FeSC");
g = MTB.read_poscar(g,"data/FeSe_TI/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%%
%===Construct xy opened model======
n1=10;
n2=10;
tic;
gs = MTB.ham.get_xy_open_wannier(g,n1,n2);
toc;
%% =====Calculate bands along z by sparse=====
% [nbands,~,nrpts]=size(gs.ham);
nbands=size(gs.ham{1},1);
nrpts=size(gs.ham,2);
labels={'Z','\Gamma','Z'}; % labels for k
hkpoints={[0.0,0.0,0.5],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=20;
Electric_field_in_evpA=0.0;
numEigs=600;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_sparse(gs.ham,gs.hopr,nbands,nrpts,efermi,hkpoints,nk,numEigs,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%% =====Calculate bands along z by full ham=====
ham=zeros(10*10*6,10*10*6,3);
hopr=zeros(3,3);
ham(:,:,1)=gs.ham{1,1};
ham(:,:,2)=gs.ham{1,2};
ham(:,:,3)=gs.ham{1,3};
gs.ham=ham;
%
[nbands,~,nrpts]=size(gs.ham);
labels={'Z','\Gamma','Z'}; % labels for k
hkpoints={[0.0,0.0,0.5],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=20;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%% =====Calculate bands along z by get_supercell_wannier way=====
%===create the tb model======
clc;
clear all;
g = MTB.geometry("FeSC");
g = MTB.read_poscar(g,"data/FeSe_TI/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%===Construct xy opened model======
n1=10;
n2=10;
tic;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
toc;
gs.hopr=gs.hopr(13:15,:)
gs.ham=gs.ham(:,:,13:15)
% =====Calculate bands along z=====
[nbands,~,nrpts]=size(gs.ham);
labels={'Z','\Gamma','Z'}; % labels for k
hkpoints={[0.0,0.0,0.5],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=20;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       BdG Ham with s wave pairing                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===create the tb model======
clc;
clear all;
g = MTB.geometry("FeSC");
g = MTB.read_poscar(g,"data/FeSe_TI/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/FeSe_TI/wannier90_hr_p1.dat','data/FeSe_TI/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
%%
%===Construct xy opened model======
n1=150;
n2=150;
tic;
gs = MTB.ham.get_xy_open_wannier(g,n1,n2);
toc;
%%
delta=0.08;
sigma_y=sparse([0, -1j; 1j, 0]);
orbital=speye(3);
orbital(2,2)=-1;
h_delta=sparse(kron(orbital,1j*sigma_y*delta));
h_delta=sparse(kron(eye(n1*n2),h_delta));
%%
nbands=size(gs.ham{1},1);
nrpts=size(gs.ham,2);
labels={'Z','\Gamma','Z'}; % labels for k
hkpoints={[0.0,0.0,0.5],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=20;
numEigs=50;
mu=0;

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG_sparse(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b,mu,h_delta,numEigs);
% MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",0);
figure()
plot(kpath,Energy,'ro')
%%
nbands=size(gs.ham{1},1);
nrpts=size(gs.ham,2);
center = mean(gs.wpos(:,1:2),1);
[theta,r]=cart2pol(gs.wpos(1:nbands/n1/n2:end,1)-center(1),gs.wpos(1:nbands/n1/n2:end,2)-center(2));
xi0=1;
theta_vec = spdiags(tanh(r/xi0).* exp(1i * theta),0,n1*n2,n1*n2);
theta_vec = spdiags(exp(1i * theta),0,n1*n2,n1*n2);
%%
delta=0.08;
sigma_y=sparse([0, -1j; 1j, 0]);
orbital=speye(3);
orbital(2,2)=-1;
h_delta=sparse(kron(orbital,1j*sigma_y*delta));
h_delta=sparse(kron(theta_vec,h_delta));
kpoint=[0,0,0.2];
mu=0;
numEigs=20;
%Energy =MTB.ham.get_votex_sparse_at_q(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b,mu,h_delta,numEigs);
%%
mus=linspace(-2,2,101);

parfor idx=1:length(mus)
    kpoint=[0,0,0.5];
    mu=mus(idx);
    tic;
    Energys(:,idx)= MTB.ham.get_votex_sparse_at_q(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b,mu,h_delta,numEigs);
    toc;
end
%%
figure()
plot(mus,Energys,'b-')
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
