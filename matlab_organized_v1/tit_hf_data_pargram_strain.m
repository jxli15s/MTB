clear;
clear all;
%%
% Step 0: Add strain to the Ham
g0 = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
% g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g0.ham = transformBasis(g0.ham, T);

% Step 3: Set Wannier Position
g0.wpos = setWannierPosition(g0);

% Step 5: Neighbor Search
[result_matrices0, ~, ~, ~] = findNeighbors(g0);
%%
strain_list=0.98:0.002:1.05;
% Start total timer
total_time = tic;
epsilons=[2 3 4 5 6 7  8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 50 55 60 70 80 90 100];
band_gaps=zeros(length(strain_list),length(epsilons));
for s = 1:length(strain_list)
    fprintf("Processing on strain %f\n",strain_list(s))
    for eps_idx=1:length(epsilons)
        fprintf("Processing on eps %d\n",epsilons(eps_idx))
        % Step 1: Initialize Geometry and Hamiltonian
        g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
        % Step 2: Basis Transformation
        T = getBasisTransformMatrix();
        g.ham = transformBasis(g.ham, T);
        strain=[1/strain_list(s)*1.00,0,0;
            0.0,strain_list(s),0.0;
            0.0,0.0,1.0];
        g.a=strain*g.a;
        % Step 3: Set Wannier Position
        g.wpos = setWannierPosition(g);
        % 计算邻近矩阵
        [result_matrices1, ~, ~, ~] = findNeighbors(g);
        for i = 2:size(result_matrices1,2)-1
            delta_ij = result_matrices1{i}(:,6) - result_matrices0{i}(:,6); % 确保 result_matrices1{i} 正确
            t_ij = exp(-1.5 * delta_ij ./ result_matrices0{i}(:,6)); % 这里 beta = 1.5, 可调
            % 添加修正因子到矩阵
            result_matrices1{i} = [result_matrices1{i}, t_ij];

            % 遍历所有邻接项
            for j = 1:length(result_matrices1{i})
                bandindex_i = result_matrices1{i}(j,1);
                bandindex_j = result_matrices1{i}(j,2);
                scale = result_matrices1{i}(j,7);
                index=find(ismember(g.hopr,result_matrices1{i}(j,3:5),'rows'));
                if ~isempty(index)
                    g.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) = ...
                        g.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) * scale;
                end
            end
        end
        g.iniham = g.ham + 0;
        % Step 4: Construct Supercell
        n1 = 15; n2 = 1; % Supercell dimensions
        gs = constructSupercell(g, n1, n2);
        gs.iniham = gs.ham + 0; % Store the initial Ham
        [nbands,~,nrpts]=size(gs.ham);
        labels={'Y','\Gamma','X','R'}; % labels for k
        hkpoints={[0.0,0.5,0.0],...
            [0.0,0.0,0.0],...
            [0.5,0.0,0.0],...
            [0.5,0.5,0.0]};% hkpoints-high symmetry k points
        % % % % Initialize Hartree-Fock States
        result=load("data/tit_hf/15x1/10nmd/U-epsilon-f/strain/U_V1_large/rs_13/"+int2str(eps_idx)+"-15x1"+"/result_s1_eps"+int2str(epsilons(eps_idx))+".00_s"+int2str(s)+".mat");
        U0=result.U0;
        U=result.V;
        V=result.V;
        xinitial=result.xinitial;
        pairsU=result.pairsU;
        pairsV=result.pairsV;
        efermi=result.efermi;
        modifyHam(gs, xinitial, V, V, pairsU, pairsV)
        Electric_field_in_evpA=0.00*0.529177; nk=101;
        [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
        En=Energy';
        % En=Enk;
        % valence_band_index=62;
        % conduction_band_index=63;
        valence_band_index=60;
        conduction_band_index=61;
        [band_gap, VBM, CBM, VBM_k, CBM_k] = compute_band_gap_1d(En, valence_band_index, conduction_band_index);
        band_gaps(s,eps_idx)=band_gap;
        % if eps_idx == length(epsilons)
        %     Ua=U0(1);
        %     v1=V(3);
        %     U_V=[U_V,Ua/v1];
        % end
    end
end
% Record and display total runtime
total_runtime = toc(total_time);
fprintf('Total runtime for all scans: %.2f seconds.\n', total_runtime);

% Save results to a .mat file for later analysis
%save('RandomUV_HartreeFock_results.mat', 'results', '-v7.3');

disp('Random U and V Hartree-Fock calculations completed.');
%%
strain_list=0.98:0.002:1.05;
% save("tit_hf_data_pargram_10nm_15x1_f.mat","U_V","band_gaps","epsilons","rs")
% save("tit_hf_data_pargram_10nm_15x1_f_vb.mat","U_V","band_gaps","epsilons","rs")
% save("tit_hf_data_pargram_10nm_15x1_new2_vb.mat","U_V","band_gaps","epsilons","rs")
% load("tit_hf_data_pargram_10nm_15x1_f.mat")
% load("tit_hf_data_pargram_10nm_15x1_f_vb.mat")
% load("tit_hf_data_pargram_10nm_15x1_new2_vb.mat")
% load("tit_hf_data_pargram_10nm_15x1.mat")
% save("tit_hf_data_pargram_10nm_15x1_f_strain.mat","rs","band_gaps","epsilons")
% load("tit_hf_data_pargram_10nm_15x1_f_strain.mat")
% save("tit_hf_data_pargram_10nm_15x1_f_strain_vb.mat","rs","band_gaps","epsilons")
% load("tit_hf_data_pargram_10nm_15x1_f_strain_vb.mat")
% save("tit_hf_data_pargram_10nm_15x1_f_strain_vb_largeU.mat","rs","band_gaps","epsilons")
% load("tit_hf_data_pargram_10nm_15x1_f_strain_vb_largeU.mat")
% save("tit_hf_data_pargram_10nm_15x1_f_strain_largeU_r13.mat","rs","band_gaps","epsilons")
% load("tit_hf_data_pargram_10nm_15x1_f_strain_largeU_r13.mat")
%%
% plot the pargram on the ori grid
figure;
imagesc(epsilons, strain_list, band_gaps); % 使用 imagesc 绘制色块图
% colormap(jet); % 使用 jet 色彩映射
colorbar; % 添加颜色条
xlabel('\epsilon (Dielectric Constant)', 'FontSize', 12); % 横轴标签
ylabel('r (Å)', 'FontSize', 12); % 纵轴标签
title('2D Band Gap Distribution', 'FontSize', 14); % 标题
% caxis([])
% 调整图像方向和显示范围
set(gca, 'YDir', 'normal'); % 确保 y 轴从小到大显示
axis tight; % 自动调整轴范围

%% plot the pargram with finer grid
% 构建 finer grid（插值网格）
rs=0.98:0.002:1.05;
% rs=1:length(strain_list);
fine_eps = linspace(min(epsilons), max(epsilons), 10000); % 更细的 epsilons 网格
fine_rs = linspace(min(rs), max(rs), 10000); % 更细的 rs 网格
[eps_grid, rs_grid] = meshgrid(epsilons, rs);
[fine_eps_grid, fine_rs_grid] = meshgrid(fine_eps, fine_rs);

% 插值 band_gaps 数据到 finer grid
fine_band_gaps = interp2(eps_grid, rs_grid, band_gaps, fine_eps_grid, fine_rs_grid, 'cubic');

% 绘制 2D 等高线分布图

figure;
imagesc(fine_eps, fine_rs, fine_band_gaps); % 使用 imagesc 绘制平滑的色块图
set(gca, 'YDir', 'normal'); % 确保 y 轴从小到大显示

colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp
% 自定义颜色条范围
caxis([-0.03, 0.03]); % 设置颜色范围（单位: eV）
% caxis([-0.15, 0.15]); % 设置颜色范围（单位: eV）
% caxis([-0.1, 0.1]); % 设置颜色范围（单位: eV）
h = colorbar;
% h.Ticks = [0, 1, 2, 3, 4, 5]; % 自定义颜色条刻度
h.Label.String = 'Band Gap (eV)'; % 设置颜色条标签

% 添加等高线
% hold on;
% contour(fine_eps_grid, fine_rs_grid, fine_band_gaps, 10, 'LineColor', 'k', 'LineWidth', 0.5); % 黑色等高线

% 设置图形标签与标题
xlabel('\epsilon (Dielectric Constant)', 'FontSize', 12);
% ylabel('r (Å)', 'FontSize', 12);
ylabel('strain', 'FontSize', 12);
title('2D Band Gap Distribution (Interpolated)', 'FontSize', 14);

% 调整显示
set(gca, 'YDir', 'normal'); % 确保 y 轴从小到大显示
xlim([9, 22]);
ylim([0.98,1.05])
% ylim([1.5,  4.7])
% axis tight;
% axis tight;

print('tit_pargram_strain.png','-dpng','-r600')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                              初始几何和Ham                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = initializeGeometry(materialName, dataPath)
    % initializeGeometry Initializes the material geometry and Hamiltonian
    %
    % Inputs:
    %   materialName - String, the name of the material (e.g., "TaIrTe4")
    %   dataPath     - String, the path to the directory containing input files
    %
    % Outputs:
    %   g            - Struct, contains the geometry and Hamiltonian information
    %
    % This function reads the material geometry from a POSCAR file and
    % constructs the Hamiltonian using Wannier90 HR files.

    % Create geometry object
    g = MTB.geometry(materialName);
    % Read geometry from POSCAR file
    g = MTB.read_poscar(g, fullfile(dataPath, "POSCAR"));
    % Read Wannier Hamiltonian data
    [g.ham, g.hopr] = MTB.wannier.read_hr(...
        fullfile(dataPath, "wannier90_hr_p1.dat"), ...
        fullfile(dataPath, "wannier90_hr_p2.dat"));
end

function T = getBasisTransformMatrix()
    % getBasisTransformMatrix Returns the basis transformation matrix
    %
    % Outputs:
    %   T - Matrix, the basis transformation matrix for Hamiltonian
    %
    % The transformation matrix `T` is used to change the basis of the
    % Hamiltonian matrices.

    %T = [1,0,0,0,0,0,0,0; ...
    %     0,0,0,0,0,1,0,0; ...
    %     0,0,1,0,0,0,0,0; ...
    %     0,0,0,0,0,0,0,1; ...
    %     0,1,0,0,0,0,0,0; ...
    %     0,0,0,0,1,0,0,0; ...
    %     0,0,0,1,0,0,0,0; ...
    %     0,0,0,0,0,0,1,0];

    T = [1,0,0,0,0,0,0,0; ...
         0,0,0,0,1,0,0,0; ...
         0,1,0,0,0,0,0,0; ...
         0,0,0,0,0,1,0,0; ...
         0,0,1,0,0,0,0,0; ...
         0,0,0,0,0,0,1,0; ...
         0,0,0,1,0,0,0,0; ...
         0,0,0,0,0,0,0,1];
end

function ham = transformBasis(ham, T)
    % transformBasis Applies a basis transformation to the Hamiltonian
    %
    % Inputs:
    %   ham - 3D array, the original Hamiltonian in the initial basis
    %   T   - Matrix, the basis transformation matrix
    %
    % Outputs:
    %   ham - 3D array, the Hamiltonian in the transformed basis

    % Apply the transformation for each Hamiltonian slice
    for i = 1:size(ham, 3)
        ham(:, :, i) = T * ham(:, :, i) * inv(T);
    end
end

function wpos = setWannierPosition(g)
    % setWannierPosition Sets the Wannier positions for the geometry
    %
    % Inputs:
    %   g    - Struct, contains the geometry and atomic positions
    %
    % Outputs:
    %   wpos - Array, the Wannier positions for each orbital
    %
    % This function computes the Wannier positions based on the atomic
    % positions and lattice vectors.

    % Compute Wannier positions based on atomic positions
    wpos = g.atoms * g.a;
    wpos = [wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :); ...
            wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :)];
end

function gs = constructSupercell(g, n1, n2)
    % constructSupercell Constructs the supercell Hamiltonian
    %
    % Inputs:
    %   g  - Struct, contains the geometry and Hamiltonian of the material
    %   n1 - Integer, number of cells along the first lattice vector
    %   n2 - Integer, number of cells along the second lattice vector
    %
    % Outputs:
    %   gs - Struct, the supercell Hamiltonian and related properties

    % Generate supercell Hamiltonian
    gs = MTB.ham.get_supercell_wannier(g, n1, n2);
    % Store the initial Hamiltonian for reference
    gs.iniham = gs.ham+0;
end


function [result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs)
    % findNeighbors Finds neighboring atoms for interactions
    %
    % Inputs:
    %   gs - Struct, contains the supercell Hamiltonian and Wannier positions
    %
    % Outputs:
    %   result_matrices - Cell array, contains neighbor data for various distances
    %   pairsU          - Cell array, contains pairs for on-site and off-site interactions (U)
    %   pairsV          - Cell array, contains pairs for on-site and off-site interactions (V)

    result_matrices = cell(1, 7);
    % Loop through different neighbor distances
    for i = 1:7
        result_matrices{i} = find_neighbor_data(gs.wpos(1:2:end, :), gs.a, i, 2);
    end

    % Filter and group pairs
    pairs_onsite = result_matrices{1}(result_matrices{1}(:, 1) == result_matrices{1}(:, 2), :);
    pairs_onsite_nn = result_matrices{1}(result_matrices{1}(:, 1) ~= result_matrices{1}(:, 2), :);

    pairs_onsite_1 = pairs_onsite(mod(pairs_onsite(:,1),4)==1,:);
    pairs_onsite_2 = pairs_onsite(mod(pairs_onsite(:,1),4)==2,:);
    pairs_onsite_3 = pairs_onsite(mod(pairs_onsite(:,1),4)==3,:);
    pairs_onsite_4 = pairs_onsite(mod(pairs_onsite(:,1),4)==0,:);
    pairs_onsite_nn_13 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==1,:);
    pairs_onsite_nn_24 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==0,:);


    pairs_offsite_nn = result_matrices{2};
    pairs_offsite_nn_13=result_matrices{2}(mod(result_matrices{2}(:, 1),2)~=0 & mod(result_matrices{2}(:, 2),2)~=0, :);
    pairs_offsite_nn_24=result_matrices{2}(mod(result_matrices{2}(:, 1),2)==0 & mod(result_matrices{2}(:, 2),2)==0, :);


    pairs_offsite_nnn = result_matrices{3};
    pairs_offsite_nnn_12=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==2,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_14=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_34=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==3,:)];
    pairs_offsite_nnn_32=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==3,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==2,:)];

    pairs_offsite_nnnn = result_matrices{4};
    pairs_offsite_nnnn_12=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==2,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_14=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_34=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==3,:)];
    pairs_offsite_nnnn_32=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==3,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==2,:)];

    pairs_offsite_nnnnn = result_matrices{5};
    pairs_offsite_nnnnn_13=result_matrices{5}(mod(result_matrices{5}(:, 1),2)~=0 & mod(result_matrices{5}(:, 2),2)~=0, :);
    pairs_offsite_nnnnn_24=result_matrices{5}(mod(result_matrices{5}(:, 1),2)==0 & mod(result_matrices{5}(:, 2),2)==0, :);


    % Construct pairs for U and V
    pairsU0 = {pairs_onsite_1, pairs_onsite_2, pairs_onsite_3, pairs_onsite_4};
    pairsU = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34,...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24};
    pairsV = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24};
end

function [xinitial, ni, si, efermi] = runHartreeFock(gs, xinitial, U, V, knum, stepmax, critial)
    % runHartreeFock Runs the Hartree-Fock self-consistent calculation
    %
    % Inputs:
    %   gs       - Struct, the supercell Hamiltonian and related properties
    %   xinitial - Cell array, initial guesses for Hartree-Fock states
    %   U        - Array, on-site Coulomb interaction values
    %   V        - Array, off-site Coulomb interaction values
    %   knum     - Integer, number of k-points along each direction
    %   stepmax  - Integer, maximum number of Hartree-Fock steps
    %   critial  - Float, convergence criterion for self-consistency
    %
    % Outputs:
    %   xinitial - Cell array, converged Hartree-Fock states
    %   ni       - Array, converged occupation numbers
    %   si       - Array, spin polarization values
    %   efermi   - Float, Fermi energy of the system

    % Generate 2D k-mesh
    [Kx, Ky, Kz] = gs.get_Bulk2Dkmesh([0, 1], [0, 1], knum);
    kpoints = [Kx(:), Ky(:), Kz(:)];

    % Set electric field to zero (default)
    Electric_field_in_evpA = 0;

    % Run the Hartree-Fock solver
    [xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, ...
        Electric_field_in_evpA, xinitial, stepmax, critial, U, V, 1e-10, 4);
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   Function to run hartree                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xinitial_U, xinitial_V] = initializeStates(U, V, pairsU, pairsV, nbands,xinitial_0)

    % initializeStates initializes the xinitial_U and xinitial_V arrays
    % Inputs:
    %   U       - Cell array or vector for U
    %   V       - Cell array or vector for V
    %   pairsU  - Cell array of pair matrices corresponding to U
    %   pairsV  - Cell array of pair matrices corresponding to V
    %   nbands  - Number of bands
    % Outputs:
    %   xinitial_U - Cell array containing initialized 3D matrices for U
    %   xinitial_V - Cell array containing initialized 3D matrices for V

    % Initialize empty cell arrays
    xinitial_U = {};
    xinitial_V = {};
    
    % Populate xinitial_U
    if ~isempty(U)
        for i = 1:length(U)
            pair = pairsU{i}; % Get the (i-1)th pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_idx = 1:num_pair
                % 提取原子对索引和分数坐标
                atom_i = pair(pair_idx, 1); % 原子 i
                atom_j = pair(pair_idx, 2);
                x=U(i);
                nj=xinitial_0(2*atom_j-1,2*atom_j-1)+xinitial_0(2*atom_j,2*atom_j);
                correlation(pair_idx,2*atom_i-1, 2*atom_i-1) = x*nj;%
                correlation(pair_idx,2*atom_i, 2*atom_i) = x*nj;%
            end
            xinitial_U = [xinitial_U, correlation.*0.1];
        end
    end
    
    % Populate xinitial_V
    if ~isempty(V)
        for i = 1:length(V)
            pair = pairsV{i}; % Get the ith pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_index = 1:num_pair
                % 提取原子对索引  
                atom_i = pair(pair_index, 1); % 原子 i
                atom_j = pair(pair_index, 2); % 原子 j
                row_up_i = 2 * atom_i - 1;
                row_dn_i = 2 * atom_i;
                col_up_j = 2 * atom_j - 1;
                col_dn_j = 2 * atom_j;
                x=V(i)*0;
                correlation(pair_index,row_up_i, col_up_j) = -x; % 上自旋部分
                correlation(pair_index,row_dn_i, col_dn_j) = -x; % 下自旋部分
                correlation(pair_index,row_dn_i, col_up_j) = 0;
                correlation(pair_index,row_up_i, col_dn_j) = 0;
                correlation(pair_index,col_up_j, row_up_i) = -x; % 上自旋部分
                correlation(pair_index,col_dn_j, row_dn_i) = -x; % 下自旋部分
                correlation(pair_index,col_dn_j, row_up_i) = 0;
                correlation(pair_index,col_up_j, row_dn_i) = 0;

            end
            xinitial_V = [xinitial_V, correlation.*0.1];
        end
    end
end


function potentials = calculate_dual_gate_potentials(r, d, epsilon, n_max)
    % calculate_dual_gate_potentials: Compute U and V values for various distances.
    % 
    % Inputs:
    %   d       - Distance between the gates (in meters)
    %   epsilon - Dielectric constant
    %   n_max   - Maximum number density
    % 
    % Output:
    %   potentials - Struct containing distances and their corresponding potentials

    % Define atomic distances (in meters)
    distances = [r, 3.77, 4.4, 6.91, 7.54, 8.66, 10.17] * 10^(-10);

    % Initialize results
    potentials = struct();
    potentials.distances = distances;
    potentials.values = zeros(size(distances));

    % Calculate dual gate potential for each distance
    for i = 1:length(distances)
        potentials.values(i) = dual_gate_potential(distances(i), d, epsilon, n_max);
    end

    % Display results
    fprintf('Distance (m) \t Potential (V)\n');
    for i = 1:length(distances)
        fprintf('%.2e \t %.2f\n', distances(i), potentials.values(i));
    end
end

function V = dual_gate_potential(r, d, epsilon, n_max)
    % Calculate dual-gate Coulomb potential using image charge method
    % r: radial distance
    % d: distance to the gates
    % q: charge
    % epsilon: dielectric constant
    % n_max: number of image charges to consider
    e=1.602176634*10^(-19);
    epsilon_0=8.854187817*10^(-12);
    k_e=1/(4*pi*epsilon_0);
    V = 0; % Initialize potential
    for n = -n_max:n_max
        V = V + ((-1)^n) / sqrt(r^2 + (2*n*d)^2);
    end
    % V = e^2/(4*pi*epsilon_0)/e/epsilon * V; % Final potential meV*m
    V=e^2*k_e*1/epsilon/e * V; %to eV
end


function [xinitial,ni,si,efermi]=runhartreev8(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)

        [xinitial, metaData] = flattenNestedCell(xinitial);        % 展平操作
        objective = @(xinitial) one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV);
        xinitial = quasi_newton(objective, xinitial, 7, critial, stepmax);
        xinitial = restoreNestedCell(xinitial, metaData);  % 复原操作
        modifyHam(gs, xinitial, U, V, pairsU, pairsV)  %修改Ham
        nbands=size(gs.ham,1);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function xnew = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)
    % 1. 基本初始化
    gs.ham=gs.iniham+0; 
    nbands = size(gs.ham, 1);
    % 恢复操作
    xinitial = restoreNestedCell(xinitial, metaData);
    xinitial_0=xinitial{1};
    if ~isempty(U)
         xinitial_U=xinitial{2};
    end 
    if ~isempty(V)
        xinitial_V=xinitial{3};
    end 

    gs.onsite_modify(xinitial_0);

    if ~isempty(U)
        for i = 1:length(U)
            pairs = pairsU{i};
            xinitial_1 = xinitial_U{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_1(j,:,:)));
            end
        end
    end

    if ~isempty(V)
        for i = 1:length(V)
            pairs = pairsV{i};
            xinitial_2 = xinitial_V{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_2(j,:,:)));
            end
        end
    end


    
    % 2. 计算波函数和能量 (含电场影响)
    [Unk, Enk] = MTB.ham.get_bulk_plane_bands_add_electric(gs, Electric_field_in_evpA, Kx, Ky, Kz);
    Unk = reshape(Unk, [nbands, nbands, knum^2]);
    Enk = reshape(Enk, [knum^2, nbands]);
    C_k_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);
    C_k = calculate_C_k_u1u2(Unk, Enk, u1, u2);

    % 3. 计算 U 和 V 的关联函数
    onsite_correlation_U =  zeros(size(xinitial_0));
    offsite_correlation_U = {};
    offsite_correlation_V = {};
    total_correlation = {};
   
    % 3.1 计算 U[0] 的关联函数 (onsite)
    if ~isempty(U0)
        for i = 1:length(U0)
            onsite_U = calculate_onsite_matrix(C_k_avg,pairsU0{i}, U0(i));
            % onsite_U = onsite_U*0.8+xinitial_0*0.2;
            onsite_correlation_U = onsite_correlation_U + onsite_U;
        end
        % To keep the TR symmetry
        % % diag_terms = diag(onsite_correlation_U)+diag(C_k_avg);
        % % onsite_correlation_U(1:nbands+1:end)=diag_terms/2;
        % % onsite_correlation_U=(onsite_correlation_U+onsite_correlation_U')/2;
    end
    
    % 3.2 计算 U[1:] 的关联函数 (offsite)
    if ~isempty(U)
        for i = 1:length(U)           
            offsite_U = calculate_offsite_U(C_k, kpoints, pairsU{i}, gs.a, U(i));
            % offsite_U = offsite_U*0.8 + xinitial{2}{i}*0.2;
            offsite_correlation_U = [offsite_correlation_U,offsite_U];
        end
    end
    
    % 3.3 计算 V 的关联函数 (offsite_V)
    if ~isempty(V)
        for i = 1:length(V)
            offsite_V = calculate_offsite_V(C_k, kpoints, pairsV{i}, gs.a, V(i));
            % offsite_V = offsite_V*0.8 + xinitial{3}{i}*0.2;         
            offsite_correlation_V = [offsite_correlation_V, offsite_V];
        end
    end
    total_correlation = {onsite_correlation_U,offsite_correlation_U,offsite_correlation_V};

    [xnew, ~] = flattenNestedCell(total_correlation);% 展平成列向量
end

function modifyHam(gs, xinitial, U, V, pairsU, pairsV)
    % modifyStates modifies onsite and offsite states using gs object
    % Inputs:
    %   gs        - Object containing methods `onsite_modify` and `offsite_modify`
    %   xinitial  - Cell array containing initial state matrices
    %   U         - Cell array or vector for U
    %   V         - Cell array or vector for V
    %   pairsU    - Cell array of pair matrices corresponding to U
    %   pairsV    - Cell array of pair matrices corresponding to V

    % Handle onsite modification for xinitial{1}
    xinitial_0 = xinitial{1};
    gs.onsite_modify(xinitial_0);

    % Handle offsite modification for U
    if ~isempty(U)
        xinitial_U = xinitial{2}; % Extract second cell for U states
        for i = 1:length(U)
            pairs = pairsU{i};  % Get pairs for the (i-1)th entry
            xinitial_1 = xinitial_U{i}; % Get initial state for U
            for j = 1:size(pairs,1)
                % Modify offsite states for U
                gs.offsite_modify([0,0,0], squeeze(xinitial_1(j, :, :)));
            end
        end
    end

    % Handle offsite modification for V
    if ~isempty(V)
        xinitial_V = xinitial{3}; % Extract third cell for V states
        for i = 1:length(V)
            pairs = pairsV{i};    % Get pairs for the ith entry
            xinitial_2 = xinitial_V{i}; % Get initial state for V
            for j = 1:size(pairs, 1)
                % Modify offsite states for V
                gs.offsite_modify(pairs(j, 3:5), squeeze(xinitial_2(j, :, :)));
            end
        end
    end
end


function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
    paulix=[0,1;1,0];pauliy=[0,-1i;1i,0];pauliz=[1,0;0,-1];
    parfor i=1:size(kindex,1)
        nki=abs(Unk(:,bandindex(i),kindex(i))).^2+nki;
        psik=reshape(Unk(:,bandindex(i),kindex(i)),[2,sitenum]);
        sxki=psik'*paulix*psik./2+sxki;
        syki=psik'*pauliy*psik./2+syki;
        szki=psik'*pauliz*psik./2+szki;
    end
    knum=size(Unk,3);
    ni=nki/knum;
    ni=reshape(ni,[2,size(ni,1)/2]);%first row spin up, second row spin down
    si=real([diag(sxki).';diag(syki).';diag(szki).']./knum);
end

function [Etot,kindex,bandindex,efermi]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    efermi = max(a,[],"all");
    % rule out the bands below the fermi level
    % % a=a(ceil(size(Enk(:),1)*0.5)+1:end);
    % % b=b(ceil(size(Enk(:),1)*0.5)+1:end);

    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end

function C_total_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)，每行是一个 k 点的能带
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    % u1: 填充比例下限 (0 <= u1 <= 1)
    % u2: 填充比例上限 (0 <= u2 <= 1)
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)

    % 获取维度信息
    [knum2, nbands] = size(Enk);  % Enk 的维度
    C_total = zeros(nbands, nbands); % 初始化总关联函数矩阵

    % 1. 计算 u1 和 u2 对应的费米能级
    E_low = calculate_ef(Enk, u1); % u1 对应的费米能级
    E_high = calculate_ef(Enk, u2); % u2 对应的费米能级

    %fprintf('能级范围: %.6f eV 到 %.6f eV\n', E_low, E_high);

    % 2. 遍历每个 k 点，累加所有 C^(m)
    for m = 1:knum2
        % 提取第 m 个 k 点的本征值和本征矢量
        E_k = Enk(m, :);          % 第 m 个 k 点的能带本征值，维度 (1, nbands)
        U_k = Unk(:, :, m);       % 第 m 个 k 点的本征矢量矩阵，维度 (nbands, nbands)

        % 筛选能级：处于 [E_low, E_high] 区间的态置为 1，其余置为 0
        W_k = diag(double(E_k >= E_low & E_k <= E_high)); % 维度 (nbands, nbands)

        % 计算关联函数矩阵 C^(m)
        C_m = conj(U_k) * W_k * transpose(U_k);

        % 累加当前 k 点的关联函数矩阵
        C_total = C_total + C_m;
    end

    % 3. 求平均：除以总的 k 点数目 knum^2
    C_total_avg = C_total / knum2;

    %fprintf('关联函数计算完成，结果已对 k 点求平均。\n');
end

function C_k = calculate_C_k_u1u2(U_k, E_k, u1, u2)
    % 计算关联函数矩阵 C_k，使用 u1 和 u2 提取占据能级范围
    %
    % 输入:
    % U_k: 本征矢量矩阵，维度 (nbands, nbands, nk)
    % E_k: 本征值矩阵，维度 (nk, nbands)
    % u1, u2: 占据比例范围
    %
    % 输出:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)

    % 获取维度
    [nbands, ~, nk] = size(U_k);
    C_k = zeros(nbands, nbands, nk);

    % 确定占据能级范围 [E_low, E_high]
    total_states = numel(E_k);
    E_sorted = sort(E_k(:)); % 将所有能级排序
    num_u1 = ceil(total_states * u1);
    num_u2 = ceil(total_states * u2);
    E_low = E_sorted(num_u1);
    E_high = E_sorted(num_u2);

    % 遍历每个 k 点
    for k = 1:nk
        % 提取当前 k 点的本征值和本征矢量
        E_k_point = E_k(k, :);
        U_k_point = U_k(:, :, k);

        % 构建占据函数 f(E_k) 范围
        f_k = double(E_k_point >= E_low & E_k_point <= E_high);

        % 计算关联函数矩阵 C_k
        C_k(:, :, k) = conj(U_k_point) * diag(f_k) * transpose(U_k_point);
    end
end

function onsite_U = calculate_onsite_matrix(correlation, pairs, U)

    % 计算 Onsite 矩阵，仅保留对角线上的 2x2 block 并进行变换
    %
    % 输入:
    % correlation: 关联函数矩阵，维度 (2*nbands, 2*nbands)
    % U: Hubbard U 参数
    %
    % 输出:
    % onsite_matrix: Onsite 矩阵，维度 (2*nbands, 2*nbands)

    % 获取轨道数目
    total_bands = size(correlation, 1);
    if mod(total_bands, 2) ~= 0
        error('关联函数矩阵的行数必须是偶数，每个轨道对应上下自旋');
    end
    nbands = total_bands / 2;

    % 定义变换矩阵 T
    T = [0, -1; 1, 0];

    % 初始化输出矩阵
    onsite_matrix = zeros(total_bands, total_bands);
    num_pairs = size(pairs,1);

    % 遍历每个对角线的 2x2 block
       % 遍历每个原子对
    for pair_idx = 1:num_pairs
        i = pairs(pair_idx,1);
        % 提取对角线上的 2x2 block
        block = correlation(2*i-1:2*i, 2*i-1:2*i);

        % 应用变换矩阵 T
        transformed_block = T * block * T';

        % 保存到输出矩阵的对应位置
        onsite_matrix(2*i-1:2*i, 2*i-1:2*i) = U * transformed_block;
        % To keep the TR symmetry
    end
    onsite_U=onsite_matrix;
end

function offsite_U = calculate_offsite_U(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_U = zeros(num_pairs, nbands, nbands);
    % offsite_U_total = zeros(nbands,nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义行/列索引
        index_i = 2 * atom_i - 1; 
        index_j = 2 * atom_j - 1; 


        % 遍历所有 k 点，计算相位修正的关联值
        ni_up_sum = 0;
        ni_dn_sum = 0;
        nj_up_sum = 0;
        nj_dn_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            %phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));
	    phase_factor = 1;
            % 提取 C_k 中的元素
            % ni_up = C_k(index_i, index_i, k) * phase_factor; % 上自旋
            % ni_dn = C_k(index_i+1, index_i+1, k) * phase_factor; % 下自旋
            nj_up = C_k(index_j, index_j, k) * phase_factor; % 上自旋
            nj_dn = C_k(index_j+1, index_j+1, k) * phase_factor; % 下自旋
            % 累加相位修正的贡献
            % ni_up_sum = ni_up_sum + ni_up;
            % ni_dn_sum = ni_dn_sum + ni_dn;
            nj_up_sum = nj_up_sum + nj_up;
            nj_dn_sum = nj_dn_sum + nj_dn;
        end

        % 平均化，并乘以 U
        % ni_up_sum_avg = real(ni_up_sum/nk) * U;
        % ni_dn_sum_avg = real(ni_dn_sum/nk) * U;
        nj_up_sum_avg = nj_up_sum/nk * U;
        nj_dn_sum_avg = nj_dn_sum/nk * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(index_i, index_i) = nj_dn_sum_avg+nj_up_sum_avg; %
        correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg+nj_dn_sum_avg; %
        % correlation_temp(index_j, index_j) = ni_dn_sum_avg+ni_up_sum_avg; %
        % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg+ni_dn_sum_avg; %
        % To keep the PT symmetry
        % % % correlation_temp(index_i, index_i) = nj_dn_sum_avg/2+nj_up_sum_avg/2; %
        % % % correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg/2+nj_dn_sum_avg/2; %
        % % % correlation_temp(index_j, index_j) = ni_dn_sum_avg/2+ni_up_sum_avg/2; %
        % % % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg/2+ni_dn_sum_avg/2; %       
        % 将当前原子对的矩阵存储到输出列表中
        offsite_U(pair_idx, :, :) = correlation_temp;
        % offsite_U_total = offsite_U_total+correlation_temp;
    end
end

function offsite_V = calculate_offsite_V(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_V = zeros(num_pairs, nbands, nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义上下自旋的行/列索引
        col_dn_i = 2 * atom_i;       % 下自旋 (dn) 行
        col_up_i = 2 * atom_i - 1;   % 上自旋 (up) 行
        row_dn_j = 2 * atom_j;       % 下自旋 (dn) 列
        row_up_j = 2 * atom_j - 1;   % 上自旋 (up) 列

        % 遍历所有 k 点，计算相位修正的关联值
        dndn_sum = 0;
        upup_sum = 0;
        updn_sum = 0;
        dnup_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

            % 提取 C_k 中的元素
            dndn = C_k(row_dn_j, col_dn_i, k) * phase_factor; % 下自旋-下自旋
            upup = C_k(row_up_j, col_up_i, k) * phase_factor; % 上自旋-上自旋
            updn = C_k(row_up_j, col_dn_i, k) * phase_factor;
            dnup = C_k(row_dn_j, col_up_i, k) * phase_factor;

            % 累加相位修正的贡献
            dndn_sum = dndn_sum + dndn;
            upup_sum = upup_sum + upup;
            updn_sum = updn_sum + updn;
            dnup_sum = dnup_sum + dnup;
        end

        % 平均化，并乘以 U
        dndn_avg = dndn_sum / nk * U;
        upup_avg = upup_sum / nk * U;
        updn_avg = updn_sum / nk * U;
        dnup_avg = dnup_sum / nk * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(col_up_i, row_up_j) = -upup_avg; % 上自旋部分
        correlation_temp(col_dn_i, row_dn_j) = -dndn_avg; % 下自旋部分
        correlation_temp(col_dn_i, row_up_j) = -updn_avg; % 上自旋部分
        correlation_temp(col_up_i, row_dn_j) = -dnup_avg; % 下自旋部分

        % 将当前原子对的矩阵存储到输出列表中
        offsite_V(pair_idx, :, :) = correlation_temp;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Function for scf convergency                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xn = quasi_newton(func, x0, q, tol, maxsteps, history)
    % 默认参数设置
    if nargin < 3, q = 10; end
    if nargin < 4, tol = 1e-7; end
    if nargin < 5, maxsteps = 100; end
    if nargin < 6, history = {}; end

    % 初始化变量
    f0 = func(x0).'; % 初始梯度，确保为列向量
    n = length(f0); % 确保 x0 和 func 输出长度一致

    alpha = 1e-3; % 初始步长
    U = zeros(n, q); % 初始化 U 矩阵
    V = zeros(n, q); % 初始化 V 矩阵

    xn = x0(:); % 确保 x0 是列向量
    fn = f0; % 初始化梯度为列向量

    disp('Starting quasi-Newton...');

    % === 预热阶段：初始化 U 和 V ===
    for i = 1:q
        fn = func(xn).'; % 确保梯度为列向量
        xn_1 = xn - alpha * fn; % 初步更新变量

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 更新 U 和 V 矩阵
        U(:, i) = fn - xn; % 梯度差
        V(:, i) = func(fn).' - fn; % 梯度映射差
        xn = xn_1; % 更新变量
    end

    % % % % % % disp('Warmed up.');

    % === 主迭代循环 ===
    for i = 1:maxsteps
        fn = func(xn).'; % 确保梯度为列向量

        % 检查维度一致性
        xn = xn(:); % 确保 xn 是列向量
        fn = fn(:); % 确保 fn 是列向量
        if size(U, 1) ~= length(xn)
            error('Dimension mismatch: U should have %d rows, but got %d.', length(xn), size(U, 1));
        end

        % 构造矩阵 C 和向量 b
        C = transpose(U) * U - transpose(U) * V; % 计算 C 矩阵
        b = transpose(U) * (xn - fn); % 计算 b 向量

        % 判断 C 是否病态，必要时正则化
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        % 求解更新方向 delta
        delta = V * (C \ b);

        % 更新变量
        xn_1 = fn - delta;

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 存储历史记录
        history{end+1} = xn_1;

        % 更新 U 和 V 矩阵
        fn_1 = func(xn_1).'; % 确保梯度为列向量
        U = [U(:, 2:end), fn_1 - xn_1]; % 滑动窗口更新 U
        V = [V(:, 2:end), func(fn_1).' - fn_1]; % 滑动窗口更新 V

        % 检查收敛条件
        dx = max(abs(xn_1 - xn)); % 计算变量更新幅度
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Converged at iteration %d\n', i);
            xn = xn_1.'; % 最终结果转置为行向量
            return;
        end

        xn = xn_1; % 更新当前变量
    end

    % 未收敛提示
    disp('Did not converge within the maximum steps.');
    xn = xn.'; % 如果未收敛，也转置为行向量
end

function [flattenedVector, metaData] = flattenNestedCell(C)
    flattenedVector = []; % 用于存储展平后的行向量
    metaData = {};        % 用于存储路径和尺寸信息

    % 递归处理单元数组
    function processElement(element, path)
        if iscell(element)
            % 如果是单元数组，递归处理每个子元素
            for i = 1:numel(element)
                processElement(element{i}, [path, i]);
            end
        else
            % 如果是普通数组，展平并记录信息
            flattenedVector = [flattenedVector, element(:)']; % 展平并拼接
            metaData{end+1} = struct('Path', {path}, 'Size', size(element)); % 记录路径和尺寸
        end
    end

    % 开始递归处理
    processElement(C, []);
end

function restoredCell = restoreNestedCell(flattenedVector, metaData)
    restoredCell = {};    % 初始化恢复后的单元数组
    currentIndex = 1;     % 当前展平向量的索引

    % 递归赋值函数
    function target = assignElement(target, path, sizeInfo)
        if isempty(path)
            % 如果路径为空，说明已经到达最终节点
            numElements = prod(sizeInfo); % 当前数组的总元素数
            reshapedArray = reshape(flattenedVector(currentIndex:currentIndex + numElements - 1), sizeInfo); % 恢复形状
            currentIndex = currentIndex + numElements; % 更新索引位置
            target = reshapedArray; % 直接返回恢复后的数组
        else
            % 处理嵌套单元数组
            idx = path(1);
            if numel(target) < idx || ~iscell(target{idx})
                target{idx} = {}; % 确保路径上的每一级是单元数组
            end
            target{idx} = assignElement(target{idx}, path(2:end), sizeInfo); % 递归处理下一层路径
        end
    end

    % 遍历元数据并恢复每个元素
    for i = 1:numel(metaData)
        restoredCell = assignElement(restoredCell, metaData{i}.Path, metaData{i}.Size);
    end
end

function [band_gap, VBM, CBM, VBM_k, CBM_k] = compute_band_gap_1d(En, valence_band_index, conduction_band_index)
% compute_band_gap_1d: 计算指定价带和导带之间的能隙及对应的 k 点位置
%
% 输入参数:
%   En - 能量矩阵，维度为 (knum, nbands)，包含能量数据
%   valence_band_index - 价带索引 (整数)
%   conduction_band_index - 导带索引 (整数)
%
% 输出参数:
%   band_gap - 带隙 (CBM - VBM)
%   VBM - 价带最高点能量值
%   CBM - 导带最低点能量值
%   VBM_k - 价带最高点对应的 k 点位置 (整数)
%   CBM_k - 导带最低点对应的 k 点位置 (整数)

    % 获取输入矩阵的维度
    [knum, nbands] = size(En);
    assert(valence_band_index <= nbands && conduction_band_index <= nbands, ...
        '价带或导带索引超出能量矩阵的范围');

    % Step 1: 计算 VBM 和 CBM 及其索引
    [VBM, VBM_k] = max(En(:, valence_band_index)); % 价带最高点
    [CBM, CBM_k] = min(En(:, conduction_band_index)); % 导带最低点

    % Step 2: 计算带隙
    band_gap = CBM - VBM;
end

