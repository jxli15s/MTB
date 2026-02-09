%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
%[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

% T=zeros(size(g.ham,1),size(g.ham,2));
% for i=1:size(g.ham,1)
%     if mod(i,2)==1
%     T(i,ceil(i/2))=1;
%     else
%     T(i,4+i/2)=1;
%     end
% end

T=[1,0,0,0,0,0,0,0;...
   0,0,0,0,0,1,0,0;...
   0,0,1,0,0,0,0,0;...
   0,0,0,0,0,0,0,1;...
   0,1,0,0,0,0,0,0;...
   0,0,0,0,1,0,0,0;...
   0,0,0,1,0,0,0,0;...
   0,0,0,0,0,0,1,0];

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)];
% g.wpos=g.wpos
g.wpos(:,3)=0
n1=15;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;

stepmax=100;
minstepmax=11;
step=1;
knum=51;
u=(4*n1+2)/nbands;
% u=1/3;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
kpoints=[Kx(:),Ky(:),Kz(:)];
critial=10^-10;
U=0.3;
V=1;

xinitial=zeros(nbands,nbands);
nsite=nbands/8;
nec1=(4*n1+4)/n1/2/2/2;
nec2=(4*n1)/n1/2/2/2;

%save("tmp.dat","tmp")
 % load("tmp.mat")
 % xinitial=tmp;
    
% load('data.mat')
xinitial=diag(kron(ones(1,n1),[nec1-0.4,nec1-0.4,nec2+0.4,nec2+0.4,nec1+0.4,nec1+0.4,nec2-0.4,nec2-0.4]))*U;
% xinitial=diag(kron(ones(1,15),[nec1,nec1,nec2,nec2,nec1,nec1,nec2,nec2]+rand(1,8)./3))*U;
% xinitial=zeros(nbands,nbands);


result_matrix_1 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 1,2);
result_matrix_1 = result_matrix_1(result_matrix_1(:,1)~=result_matrix_1(:,2),:);
result_matrix_2 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 2,2); %% Nearest neighbor atom1-atom1 atom2-atom2
result_matrix_3 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 3,2); %% Nearest neighbor atom1-atom2
result_matrix_4 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 4,2); %%
result_matrix_5 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 5,2); %% Nearest neighbor atom1-atom1 atom2-atom2
result_matrix_6 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 6,2); %% Nearest neighbor atom1-atom2
result_matrix_7 = find_neighbor_data(gs.wpos(1:2:end,:),gs.a, 7,2); %% Nearest neighbor atom1-atom2
result_matrix=[result_matrix_1;result_matrix_3];
result_matrix=result_matrix_1;
% result_matrix=result_matrix(result_matrix(:,1)~=result_matrix(:,2),:)
% result_matrix=result_matrix(result_matrix(:,1)==result_matrix(:,2),:);



% % % % % % result_matrix = result_matrix(mod((result_matrix(:,2)-result_matrix(:,1)),4)==0,:);
% % % % % C_k = calculate_C_k_u1u2(Unk, Enk, u1, u2);
% % % % % % offsite_U_list = calculate_offsite_U_modified(C_k, kpoints, result_matrix, gs.a, U);
% % % % % offsite_U = calculate_offsite_U(C_k, kpoints, result_matrix, gs.a, U);
% % % % % %%
% % % % % offsite_V = calculate_offsite_V(C_k, kpoints, result_matrix, gs.a, U);
% % % %%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%                          Run HF  version5                         %%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % u1=0.5000001;u2=u;
% % % U=1.0;
% % % %[xinitial,ni,si,efermi]=runhartreev4(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u);
% % % % [xinitial,ni,si,efermi]=runhartreev6(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2,result_matrix);
% % % [xinitial,ni,si,efermi]=runhartreev5(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Run HF new version7                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1=1e-10;u2=u;
% % U=[1.0,0.5,0.4,0.2];
% % V=[-0.5,-0.4,-0.2];
U=[1.0,0.5,0.4,0.2];
V=[-0.5,-0.4,-0.2];
% result_matrix=[result_matrix;result_matrix_1(result_matrix_1(:,1)~=result_matrix_1(:,2),:)];
pairs_onsite_nn=result_matrix_1(result_matrix_1(:,1)~=result_matrix_1(:,2),:); %On site Hubbard U for different orbitals
pairs_onsite_nn_1= pairs_onsite_nn(~mod(pairs_onsite_nn(:,1),2)==0,:);% for 1-3
pairs_onsite_nn_2= pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==0,:); % for 2-4
pairs_offsite_nn=result_matrix_2;
pairs_offsite_nnn=result_matrix_3;
pairs_offsite_nnn_1=pairs_offsite_nnn(mod(pairs_offsite_nnn(:,1),4)==1 & mod(pairs_offsite_nnn(:,2),4)==2,:); %for 1-2
pairs_offsite_nnn_2=pairs_offsite_nnn(mod(pairs_offsite_nnn(:,1),4)==3 & mod(pairs_offsite_nnn(:,2),4)==0,:); %for 3-4

pairs={pairs_onsite_nn_1,pairs_onsite_nn_2,pairs_offsite_nnn_1,pairs_offsite_nnn_2};
pairs={result_matrix_1,result_matrix_1};
pairsU={pairs_onsite_nn,pairs_offsite_nn,pairs_offsite_nnn};
pairsV={pairs_onsite_nn,pairs_offsite_nn,pairs_offsite_nnn};
% Initial states

xinitial_U={};
xinitial_V={};
if length(U) > 1
    for i = 2:length(U)
            pair=pairsU{i-1};
            xinitial_U=[xinitial_U,zeros(size(pair,1),nbands,nbands)];
    end
end

if ~isempty(V)
    for i = 1:length(V)
        pair = pairsV{i};
        xinitial_V=[xinitial_V,zeros(size(pair,1),nbands,nbands)];
    end
end
xinitial_0=xinitial;
xinitial={xinitial_0,xinitial_U,xinitial_V};

[xinitial,ni,si,efermi]=runhartreev8(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial,U,V,u1,u2,pairsU,pairsV);

%%
% save('data.mat', 'xinitial');
% load('data.mat', 'xinitial');
% xinitial=zeros(nbands,nbands);
% xinitial=diag(diag(xinitial))
gs.onsite_modify(xinitial_0);
%%
% % [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
% % Unk=reshape(Unk,[nbands,nbands,knum^2]);
% % Enk=reshape(Enk,[knum^2,nbands]);
% % efermi=calculate_ef(Enk,0.5);
%%
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
 Unk=reshape(Unk,[nbands,nbands,knum^2]);
Enk=reshape(Enk,[knum^2,nbands]);
[~,kindex,bandindex,efermi]=Total_energy(Enk,u);
[ni,si]=calonsite(Unk,kindex,bandindex);
%%
Electric_field_in_evpA=0.00*0.529177;
nk=101;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
% ylim([-0.4,0.4])
hold on;
plot(kpath,Energy(4*n1,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+1,:)-efermi,"Color",'red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
plot(kpath,Energy(4*n1+3,:)-efermi,"Color",'yellow','LineWidth',2);
% Energy_ori=Energy;
% % save('energyori.mat',"Energy_ori")
% 
% load('energyori.mat')
% % for i=1:size(Energy_ori,1)
%     plot(kpath,Energy_ori(i,:)-0.1562,'Color','red','LineWidth',2);
%     % hold on
% % end
%% 0.0931 for n=0.2
%% 0.0747 for n=0.1
%% 0.1212 for n=0.3
%% 0.1562 for n=0.4
%% 0.3908 for n =0.6
%% 0.6378 for n =1

% Energy_ori=Energy;

load('energyori.mat')
for i=1:size(Energy_ori,1)
    plot(kpath,Energy_ori(i,:)-0.678,'Color','red','LineWidth',2);
    % hold on
end
%%
knum=101;
band1=1;
band2=4*n1;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);

knum=101;
band1=4*n1+1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);



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
function [xinitial,T_energy,ni,si,efermi]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    xorders=zeros(size(xinitial,1),size(xinitial,2),5);

    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end

        tap=sum(abs(xnew-xinitial),"all");
 
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end

        xinitial=xnew*0.8+0.2*xinitial;

        % % % step_mod=mod(step-1,4)+1; 
        % % % xorders(:,:,step_mod)=xnew;
        % % % if step<5
        % % %     xinitial=xnew*0.8+0.2*xinitial;
        % % % else
        % % %     xinitial=(xinitial+xorders(:,:,1)+xorders(:,:,2)+xorders(:,:,3)+xorders(:,:,4)+xorders(:,:,5))/6;
        % % % end
        step=step+1;
    end
end

function [xinitial,T_energy,ni,si,efermi]=runhartreev2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,V,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial)); 
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U*1.2;
        end

        Vfock=caldensity(Unk,kindex,bandindex);
        % size(Vfock)
        for k=1:nsite/4
            Vfock(8*(k-1)+1,8*(k-1)+1)=ni(2,4*(k-1)+1+2);
            Vfock(8*(k-1)+2,8*(k-1)+2)=ni(1,4*(k-1)+1+2);
            Vfock(8*(k-1)+3,8*(k-1)+3)=ni(2,4*(k-1)+1+3);
            Vfock(8*(k-1)+4,8*(k-1)+4)=ni(1,4*(k-1)+1+3);
            Vfock(8*(k-1)+5,8*(k-1)+5)=ni(2,4*(k-1)+1);
            Vfock(8*(k-1)+6,8*(k-1)+6)=ni(1,4*(k-1)+1);            
            Vfock(8*(k-1)+7,8*(k-1)+7)=ni(2,4*(k-1)+1+1);
            Vfock(8*(k-1)+8,8*(k-1)+8)=ni(1,4*(k-1)+1+1);
        end
        % size(Vfock)
        
        xnew=xnew+Vfock.*V;

        tap=sum(abs(xnew-xinitial),"all");

    
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end
        xinitial=xnew*0.8+0.2*xinitial;
        step=step+1;
    end
end

function [xinitial,ni,si,efermi]=runhartreev3(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u)
        nbands=size(gs.ham,1);
        xinitial = [real(xinitial(:));imag(xinitial(:))];
        % 优化选项 for fminunc
        options = optimoptions('fminunc', ...
            'Algorithm', 'quasi-newton', ... % 使用拟牛顿方法
            'Display', 'iter', ...           % 显示迭代信息
            'OptimalityTolerance', critial, ...
            'MaxIterations', stepmax);
        objective = @(xinitial) one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        [xinitial_opt,fval] = fminunc(objective,xinitial,options);


        % % swarm_size=50;
        % % lb = -1;
        % % ub = 1;
        % % num_variables=2*nbands^2;
        % % lower_bound=lb*ones(num_variables,1);
        % % upper_bound=ub*ones(num_variables,1);
        % % options = optimoptions('particleswarm', ...
        % % 'SwarmSize', swarm_size, ...        % 粒子数量
        % % 'MaxIterations', stepmax, ...% 最大迭代次数
        % % 'Display', 'iter', ...              % 显示每次迭代的信息
        % % 'PlotFcn', 'pswplotbestf');         % 绘制优化过程
        % % objective = @(xinitial) one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        % % [xinitial_opt,fval] = particleswarm(objective,num_variables,lower_bound,upper_bound,options);



        xinitial_real=reshape(xinitial_opt(1:nbands^2),nbands,nbands);
        xinitial_imag=reshape(xinitial_opt(1+nbands^2:end),nbands,nbands);
        xinitial=xinitial_real+1j*xinitial_imag;


        disp('Relaxed xinitial:');
        disp(xinitial);

        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function [xinitial,ni,si,efermi]=runhartreev4(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u)

        nbands=size(gs.ham,1);
        xinitial = xinitial(:);
        objective = @(xinitial) one_step_hf_v2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        xinitial = quasi_newton(objective, xinitial);

        % disp('Relaxed xinitial:');
        % disp(xinitial);
        xinitial=reshape(xinitial,nbands,nbands);

        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function [xinitial,ni,si,efermi]=runhartreev5(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2)

        nbands=size(gs.ham,1);
        xinitial = xinitial(:);
        objective = @(xinitial) one_step_hf_v3(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u1,u2);
        xinitial = quasi_newton(objective, xinitial);
        % disp('Relaxed xinitial:');
        % disp(xinitial);
        xinitial=reshape(xinitial,nbands,nbands);
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function [xinitial,ni,si,efermi]=runhartreev6(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial,U,u1,u2,pairs)

        nbands=size(gs.ham,1);
        xinitial = xinitial(:);
        objective = @(xinitial) one_step_hf_v4(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,U,u1,u2,pairs);
        xinitial = quasi_newton(objective, xinitial);
        % disp('Relaxed xinitial:');
        % disp(xinitial);
        xinitial=reshape(xinitial,nbands,nbands);
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end


function [xinitial,ni,si,efermi]=runhartreev7(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial,U,V,u1,u2,pairsU,pairsV)
        nbands=size(gs.ham,1);
        % 展平操作
        [xinitial, metaData] = flattenNestedCell(xinitial);
        objective = @(xinitial) one_step_hf_v5(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA, xinitial,metaData,U,V,u1,u2, pairsU,pairsV);
        xinitial = quasi_newton(objective, xinitial, 7, critial, stepmax);
        
        % disp('Relaxed xinitial:');
        % disp(xinitial);
        xinitial=reshape(xinitial,nbands,nbands);
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function [xinitial,ni,si,efermi]=runhartreev8(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial,U,V,u1,u2,pairsU,pairsV)
        nbands=size(gs.ham,1);
        % 展平操作
        [xinitial, metaData] = flattenNestedCell(xinitial);
        objective = @(xinitial) one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData,U,V,u1,u2, pairsU,pairsV);
        xinitial = quasi_newton(objective, xinitial, 7, critial, stepmax);


        % disp('Relaxed xinitial:');
        % disp(xinitial);

        xinitial = restoreNestedCell(xinitial, metaData);
        xinitial_0=xinitial{1};
        xinitial_U=xinitial{2};
        if ~isempty(V)
            xinitial_V=xinitial{3};
        end

        gs.onsite_modify(xinitial_0);

        if length(U) > 1
            for i = 2:length(U)
                pairs = pairsU{i-1};
                xinitial_1 = xinitial_U{i-1};
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
        % gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function fval = one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial_0,U,u)
        nbands=size(gs.ham,1);
        xinitial_real=reshape(xinitial_0(1:nbands^2),nbands,nbands);
        xinitial_imag=reshape(xinitial_0(1+nbands^2:end),nbands,nbands);
        xinitial=xinitial_real+1j*xinitial_imag;
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        % T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        nsite=nbands/2;
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end       
        xnew=xnew*0.8+0.2*xinitial;
        fval = norm(xnew - xinitial, 'fro')^2;
        % fval=sum(abs(xnew-xinitial),"all");
end

function xnew = one_step_hf_v2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial_0,U,u)
        nbands=size(gs.ham,1);
        xinitial_0=reshape(xinitial_0,nbands,nbands);
        gs.onsite_modify(xinitial_0);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        % T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        nsite=nbands/2;
        xnew=zeros(size(xinitial_0));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end       
        xnew=xnew*0.8+0.2*xinitial_0;
        xnew=xnew(:);
end

function xnew = one_step_hf_v3(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial_0,U,u1,u2)
        nbands=size(gs.ham,1);
        xinitial_0=reshape(xinitial_0,nbands,nbands);
        gs.onsite_modify(xinitial_0);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        correlation = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);
        % onsite_U = calculate_onsite_U(correlation, U);
        onsite_U = calculate_onsite_matrix(correlation, U);
        xnew=onsite_U*0.8+0.2*xinitial_0;
        xnew=xnew(:);
end

function xnew = one_step_hf_v4(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial_0,U,u1,u2,pairs)
        nbands=size(gs.ham,1);
        xinitial_0=reshape(xinitial_0,nbands,nbands);
        gs.onsite_modify(xinitial_0);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        
        % correlation = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);
        % onsite_U = calculate_onsite_U(correlation, U);
        C_k = calculate_C_k_u1u2(Unk, Enk, u1, u2);
        offsite_U = calculate_offsite_U(C_k, kpoints, pairs, gs.a, U);
        
        % xnew=(onsite_U+offsite_U)*0.8+0.2*xinitial_0;
        xnew=offsite_U*0.8+0.2*xinitial_0;
        xnew=xnew(:);
end

function xnew = one_step_hf_v5(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData,U,V,u1,u2, pairsU,pairsV)
    % 1. 基本初始化
    gs.ham=gs.iniham+0; 
    nbands = size(gs.ham, 1);
    % 恢复操作
    xinitial = restoreNestedCell(xinitial, metaData);
    xinitial_0=xinitial{1};
    xinitial_U=xinitial{2};
    if ~isempty(V)
        xinitial_V=xinitial{3};
    end 

    gs.onsite_modify(xinitial_0);

    if length(U) > 1
        for i = 2:length(U)
            pairs = pairsU{i-1};
            xinitial_1 = xinitial_U{i-1};
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
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_2(j,:,:)))
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
    if ~isempty(U)
        onsite_U = calculate_onsite_matrix(C_k_avg, U(1));
        onsite_U = onsite_U*0.8 + xinitial_0*0.2;
        onsite_correlation_U = onsite_U;
    end
    
    % 3.2 计算 U[1:] 的关联函数 (offsite)
    if length(U) > 1
        for i = 2:length(U)           
            offsite_U = calculate_offsite_U(C_k, kpoints, pairsU{i-1}, gs.a, U(i));
            offsite_U = offsite_U*0.8 + xinitial_U{i-1}*0.2;
            offsite_correlation_U = [offsite_correlation_U,offsite_U];
        end
    end
    
    % 3.3 计算 V 的关联函数 (offsite_V)
    if ~isempty(V)
        for i = 1:length(V)
            offsite_V = calculate_offsite_V(C_k, kpoints, pairsV{i}, gs.a, V{i});
            offsite_V = offsite_V*0.8 + xinitial_V{i-1}*0.2;
            offsite_correlation_V = [offsite_correlation_V, offsite_V];
        end
    end

    total_correlation = {onsite_correlation_U,offsite_correlation_U,offsite_correlation_V};

    [total_correlation, ~] = flattenNestedCell(total_correlation);% 展平成列向量
    [xinitial,~] = flattenNestedCell(xinitial);% 展平成列向量
    % 4. 混合更新 Hartree-Fock 参数
    xnew = 0.8 * total_correlation + 0.2 * xinitial;
end


function xnew = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData,U,V,u1,u2, pairsU,pairsV)
    % 1. 基本初始化
    gs.ham=gs.iniham+0; 
    nbands = size(gs.ham, 1);
    % 恢复操作
    xinitial = restoreNestedCell(xinitial, metaData);
    xinitial_0=xinitial{1};
    xinitial_U=xinitial{2};
    if ~isempty(V)
        xinitial_V=xinitial{3};
    end 

    gs.onsite_modify(xinitial_0);

    if length(U) > 1
        for i = 2:length(U)
            pairs = pairsU{i-1};
            xinitial_1 = xinitial_U{i-1};
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
    if ~isempty(U)
        onsite_U = calculate_onsite_matrix(C_k_avg, U(1));
        onsite_correlation_U = onsite_U;
    end
    
    % 3.2 计算 U[1:] 的关联函数 (offsite)
    if length(U) > 1
        for i = 2:length(U)           
            offsite_U = calculate_offsite_U(C_k, kpoints, pairsU{i-1}, gs.a, U(i));
            offsite_correlation_U = [offsite_correlation_U,offsite_U];
        end
    end
    
    % 3.3 计算 V 的关联函数 (offsite_V)
    if ~isempty(V)
        for i = 1:length(V)
            offsite_V = calculate_offsite_V(C_k, kpoints, pairsV{i}, gs.a, V(i));
            offsite_correlation_V = [offsite_correlation_V, offsite_V];
        end
    end
    total_correlation = {onsite_correlation_U,offsite_correlation_U,offsite_correlation_V};

    [xnew, ~] = flattenNestedCell(total_correlation);% 展平成列向量
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

function Vfock=caldensity(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    Vfock=zeros(size(Unk,1),size(Unk,1));
    knum=size(Unk,3);
    for i =1:size(Unk,1)/8
        for j=1:2
            t11=zeros(size(Unk,1),size(Unk,1));
            t11(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5)=1;
            % t11_2=t11';
            t12=zeros(size(Unk,1),size(Unk,1));
            t12(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5)=1;
            % t12_2=t11';
            t21=zeros(size(Unk,1),size(Unk,1));
            t21(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5+1)=1;
            % t21_2=t11';
            t22=zeros(size(Unk,1),size(Unk,1));
            t22(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5+1)=1;
            % t22_2=t11';
            h11=0;h12=0;h21=0;h22=0;
            parfor k=1:size(kindex,1)
                unk=Unk(:,bandindex(k),kindex(k))
                h11=h11+unk'*t11*unk
                h12=h12+unk'*t12*unk
                h21=h21+unk'*t21*unk
                h22=h22+unk'*t22*unk
            end
            hh=[h11,h12;h21,h22]./knum;
            Vfock(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5)=hh(1,1);
            Vfock(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5+1)=hh(1,2);
            Vfock(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5)=hh(2,1);
            Vfock(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5+1)=hh(2,2);
        end
    end
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

function correlation_modified = calculate_correlation_with_phase(C_k, kpoints, R_vector)
    % 修正 correlation 矩阵，乘以相位因子 e^{ikR}
    %
    % 输入:
    % C_k: 每个 k 点的关联函数矩阵，维度 (2*nbands, 2*nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % R_vector: 原子对之间的距离向量，维度 (1, 3)
    %
    % 输出:
    % correlation_modified: 修正后的 correlation 矩阵，维度 (2*nbands, 2*nbands)

    % 获取 k 点数量和矩阵大小
    [total_bands, ~, nk] = size(C_k);
    correlation_modified = zeros(total_bands, total_bands);

    % 遍历所有 k 点，累加带相位因子的贡献
    for k = 1:nk
        % 计算相位因子 e^{ikR}
        phase_factor = exp(1i * dot(kpoints(k, :), R_vector));

        % 加权累加到修正的 correlation 矩阵
        correlation_modified = correlation_modified + C_k(:, :, k) * phase_factor;
    end

    % 取实部，确保最终矩阵是实数
    correlation_modified = (correlation_modified+correlation_modified')/2;
end

function onsite_U_old = calculate_onsite_U(correlation, U)
    % 使用块对角矩阵进行基矢变换，计算 Onsite 自旋关联矩阵
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

    % 定义小的变换矩阵 T (2x2)
    T = [0, -1; 1, 0];

    % 构建大块对角矩阵 T_large = kron(I_nbands, T)
    I_nbands = eye(nbands);           % 单位矩阵，维度为 (nbands, nbands)
    T_large = kron(I_nbands, T);      % Kronecker 积
    mask_matrix=kron(I_nbands,ones(2));
    % 基矢变换: T_large * correlation * T_large'
    onsite_U_old = U * (T_large * correlation * T_large').*mask_matrix;
    % To keep the TR symmetry
    diag_terms = diag(onsite_U_old)+diag(correlation);
    onsite_U_old(1:total_bands+1:end)=diag_terms/2;
    % onsite_U_old=diag(diag_terms/2);
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

function onsite_U = calculate_onsite_matrix(correlation, U)

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

    % 遍历每个对角线的 2x2 block
    for i = 1:nbands
        % 提取对角线上的 2x2 block
        block = correlation(2*i-1:2*i, 2*i-1:2*i);

        % 应用变换矩阵 T
        transformed_block = T * block * T';

        % 保存到输出矩阵的对应位置
        onsite_matrix(2*i-1:2*i, 2*i-1:2*i) = U * transformed_block;
        % To keep the TR symmetry
    end
    onsite_U=onsite_matrix;
    % To keep the TR symmetry
    diag_terms = diag(onsite_U)+diag(correlation);
    onsite_U(1:total_bands+1:end)=diag_terms/2;
    onsite_U=(onsite_U+onsite_U')/2;
    % onsite_U=diag(diag_terms/2);
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
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

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
        nj_up_sum_avg = real(nj_up_sum/nk) * U;
        nj_dn_sum_avg = real(nj_dn_sum/nk) * U;
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
        correlation_temp(col_up_i, row_up_j) = upup_avg; % 上自旋部分
        correlation_temp(col_dn_i, row_dn_j) = dndn_avg; % 下自旋部分
        correlation_temp(col_dn_i, row_up_j) = updn_avg; % 上自旋部分
        correlation_temp(col_up_i, row_dn_j) = dnup_avg; % 下自旋部分

        % 将当前原子对的矩阵存储到输出列表中
        offsite_V(pair_idx, :, :) = correlation_temp;
    end
end









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          Function to add E-field and Potential                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj=add_elec(obj,Electric_field_in_evpA)
    % obj.wpos(:,3)=round(obj.wpos(:,3));
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        % obj.wpos(i,3)*Electric_field_in_evpA
    end 
end

function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(1,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,1),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=0;
       V=Vamp.*(cos(2*pi/2/a*x+phi));      
 end
 gs.iniham=gs.ham+0;  
end
%%
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

    disp('Warmed up.');

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



function xn = quasi_newton_old(func, x0, q, tol, maxsteps, history)
    if nargin < 3
        q = 10;
    end
    if nargin < 4
        tol = 1e-7;
    end
    if nargin < 5
        maxsteps = 100;
    end
    if nargin < 6
        history = {};
    end

    f0 = func(x0);
    n = length(f0);

    alpha = 1e-3;
    U = zeros(n, q);
    V = zeros(n, q);

    xn = x0;
    fn = f0;

    disp('Starting quasi-Newton...');

    % Warm-up: Build U and V
    for i = 1:q
        fn = func(xn);
        xn_1 = xn - alpha * fn;

        % Check if vector can be reshaped into a square matrix
        N = sqrt(length(xn_1));

        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        U(:, i) = fn - xn;
        V(:, i) = func(fn) - fn;
        xn = xn_1;
    end

    disp('Warmed up.');

    % Main loop
    for i = 1:maxsteps
        fn = func(xn);

        % Compute matrix C
        C = transpose(U) * U - transpose(U) * V; % 使用普通转置

        b = transpose(U) * (xn - fn); % 使用普通转置

        % Regularization if singular
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        delta = V * (C \ b);
        xn_1 = fn - delta;

        % Check if vector can be reshaped into a square matrix
        N = sqrt(length(xn_1));

        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        history{end+1} = xn_1;

        fn_1 = func(xn_1);
        U = [U(:, 2:end), fn_1 - xn_1];
        V = [V(:, 2:end), func(fn_1) - fn_1];

        dx = max(abs(xn_1 - xn));
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Finished at iteration %d\n', i);
            return;
        end

        xn = xn_1;
    end

    disp('DIDN''T CONVERGE!!!!');
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
