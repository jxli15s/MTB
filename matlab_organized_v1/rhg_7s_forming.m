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
hopr=[];
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

Electric_field_in_evpA=0.001; %0.001*20.1 for20.1meV
g=add_elec(g,Electric_field_in_evpA);
efermi=get_ef(g);
%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[1/3,2/3,0.0]*0.8,...
          [1/3,2/3,0.0],...
          [1/3,2/3,0.0]+([0.5,0.5,0.0]-[1/3,2/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;
% efermi=0.0088
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
knum=1001;
% kxline=[-1,1];
% kyline=[-1,1];
x=1/3;
y=2/3;
kxline=[-0.015+x,0.015+x];
kyline=[-0.015+y,0.015+y];

% kxline=[0,0.5];
% kyline=[0,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);

% efermi = calculate_ef(Enk(:), 0.5);
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
% surf(kx, ky, Ez1, 'EdgeColor', 'none');
Ez2 = squeeze(Enk(:,:,band_index+1))';
% surf(kx, ky, Ez2, 'EdgeColor', 'none');
surf(kx, ky, Ez2-Ez1);
% colormap blue;        % 色图
colormap(slanCM('RdBu'))
colorbar;
shading interp
clim([0,0.001])
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
% zlim([-0.05,0.05])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Get the band structure by E_field                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = build_NGra_slab_tb('N',7,'t2',-0.007,'use_offset',true);
% Electric_field_in_evpA=0.001; %0.001*20.1 for20.1meV
Electric_field_in_evpA=0.1; %0.001*20.1 for20.1meV
g=add_elec(g,Electric_field_in_evpA);
efermi=get_ef(g);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.8,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
labels={'Kl','K','Kr'}; % labels for k
hkpoints={[-1/3,1/3,0.0]-[0.5,0.5,0.0]*0.2,...
          [-1/3,1/3,0.0],...
          [-1/3,1/3,0.0]+[0.5,0.5,0.0]*0.2,...
          };% hkpoints-high symmetry k points
nk=251;
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% Electric_field_in_evpA=0.0;
[Enk,Unk,kpath,kindex]=MTB.ham.get_bulk_bands_psi_add_electric(g.ham,g.hopr,g.wpos,0,nbands,nrpts,hkpoints,nk,g.a,g.b);
Energy=Enk;
Energy=Energy;

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
%%
hold on
plot(kpath,Energy(7,:)-efermi,'Color','blue','LineWidth',2);
plot(kpath,Energy(8,:)-efermi,'Color','red','LineWidth',2);
ylabel("Energy (meV)")
%%
nlayer=7;
nband=2*nlayer;
M = zeros(nlayer,2*nlayer);
for L = 1:nlayer
    M(L, 2*L-1:2*L) = 1;   % 每层两个轨道
end
A = abs(Unk).^2;                 % 14 x 14 x nk
A2 = reshape(A, nband, []);         % 14 x (14*nk)
W2 = M * A2;                     % 7  x (14*nk)
w_layer = reshape(W2, nlayer, nband, []);% 7 x 14 x nk
% w_layer(L,n,ik) = 第ik个k点第n条band在第L层的权重

%%
[w_layer, M] = get_layer_weight_fromUnk(Unk, 7);
%%
filename="rhg_7s_band/bands.dat";
fid=fopen(filename,'w');
fprintf(fid,'K_distance\t Energy\t w1\t w2\t w3\t w4\t w5\t w6\t w7\n');
fclose(fid);
for iband=1:nbands
   outlist=[kpath',Energy(iband,:)',squeeze(w_layer(1,iband,:)),squeeze(w_layer(2,iband,:)),squeeze(w_layer(3,iband,:)),...
       squeeze(w_layer(4,iband,:)),squeeze(w_layer(5,iband,:)),squeeze(w_layer(6,iband,:)),squeeze(w_layer(7,iband,:))];
   writeoutput(filename,outlist)
end

%%
%%
Efield_list=linspace(0,0.1,251);% 4100
outdir='rhg_7s_band';
for E_idx=1:length(Efield_list)
    g = build_NGra_slab_tb('N',7,'t2',-0.007,'use_offset',true);
    Electric_field_in_evpA=Efield_list(E_idx); %0.001*20.1 for20.1meV
    g=add_elec(g,Electric_field_in_evpA);
    efermi=get_ef(g);
    % [Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
    [Enk,Unk,kpath,kindex]=MTB.ham.get_bulk_bands_psi_add_electric(g.ham,g.hopr,g.wpos,0,nbands,nrpts,hkpoints,nk,g.a,g.b);
    Energy=Enk;
    [w_layer, M] = get_layer_weight_fromUnk(Unk, 7);
    % MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
    ylim([-0.1,0.1])

    fig = figure('Visible','off');
    ax  = axes(fig);
    hold(ax,'on');
    kk=kindex;
    linesize=2;
    xrange=[kpath(1),kpath(end)];
    for i=1:length(Energy(:,1))
        plot(kpath,Energy(i,:)-efermi,'Color','black','LineWidth',linesize);
        hold on
    end
    plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',2)
    for i=1:length(kk)-2
        plot([kk(i+1) kk(i+1)],[min(min(Energy))-efermi-1 max(max(Energy))-efermi+1],'--k','LineWidth',linesize)
    end
    grid off
    box on
    xlim(xrange)
    xticks(kk)
    xticklabels(labels)
    ylim([-2 1])
    ylabel('Energy (eV)','FontSize',24)

    E_meV_int = round(Electric_field_in_evpA * 20.1 * 1000);

    title(ax, sprintf('E_{gap} = %04d meV', E_meV_int));
    box(ax,'on');

    % Transparent backgrounds (works best for svg/pdf; png may depend on version)
    set(fig,'Color','none');
    set(ax,'Color','none');

    % Output filename
    filename = fullfile(outdir, sprintf('Band_%03d_E_%04dmeV.svg', E_idx-1, E_meV_int));
    % Export
    exportgraphics(ax, filename, 'BackgroundColor','none');
    % Close figure to avoid memory buildup
    close(fig);

    filename = fullfile(outdir, sprintf('Band_%03d_E_%04dmeV.dat', E_idx-1, E_meV_int));

    fid=fopen(filename,'w');
    fprintf(fid,'K_distance\t Energy\t w1\t w2\t w3\t w4\t w5\t w6\t w7\n');
    fclose(fid);
    for iband=1:nbands
        outlist=[kpath',Energy(iband,:)'-efermi,squeeze(w_layer(1,iband,:)),squeeze(w_layer(2,iband,:)),squeeze(w_layer(3,iband,:)),...
            squeeze(w_layer(4,iband,:)),squeeze(w_layer(5,iband,:)),squeeze(w_layer(6,iband,:)),squeeze(w_layer(7,iband,:))];
        writeoutput(filename,outlist)
    end    
end


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


function ef=get_ef(g)
    knum=401;%501
    % kxline=[-0.5,0.5];
    % kyline=[-0.5,0.5];
    kxline=[-0.1+1/3,0.1+1/3];
    kyline=[-0.1+2/3,0.1+2/3];
    u=0.5;
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
    ef=calculate_ef(Enk(:), u);
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



function obj=add_elec(obj,Electric_field_in_evpA)
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        obj.wpos(i,3)*Electric_field_in_evpA;
    end 
end

function g = build_NGra_slab_tb(varargin)
%BUILD_NGRA_SLAB_TB  Build geometry + real-space TB for N-layer "NGra" slab.
%
% Usage:
%   g = build_NGra_slab_tb();
%   g = build_NGra_slab_tb('N',7,'a0',2.46,'c',62,'d_inter',3.35, ...
%                          'neighbors',[2 5 6 21], 't2',-0.007, ...
%                          'use_offset',true, 'offset1',0.001, 'offset2',1.003/3);
%
% Notes:
% - Requires MTB.geometry, find_neighbor_data, sk_appro_bab (and/or sk_appro).
% - g.atoms stored in fractional coordinates; g.wpos in Cartesian.
% - g.ham is [nbands x nbands x nrpts], with g.hopr listing fractional hoppings.

% -------------------------
% 0) Parse inputs
% -------------------------
p = inputParser;
p.addParameter('N', 7, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('a0', 2.46, @(x)isnumeric(x)&&isscalar(x)&&x>0);     % in-plane lattice constant
p.addParameter('c', 62, @(x)isnumeric(x)&&isscalar(x)&&x>0);       % out-of-plane lattice parameter (vacuum)
p.addParameter('d_inter', 3.35, @(x)isnumeric(x)&&isscalar(x)&&x>0); % interlayer distance (same unit as a0)
p.addParameter('neighbors', [2 5 6 21], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('t2', -0.007, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('use_offset', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('offset1', 0.001, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('offset2', 1.003/3, @(x)isnumeric(x)&&isscalar(x)); % second sublattice x/y offset base

% choose SK function handle
p.addParameter('sk_fun', @sk_appro_bab, @(f)isa(f,'function_handle'));

p.parse(varargin{:});
opt = p.Results;

Nlay = opt.N;
a0   = opt.a0;
c    = opt.c;
d0_frac = (opt.d_inter / c);     % fractional spacing along z
nei_list = opt.neighbors;        % [2 5 6 21]
t2 = opt.t2;
sk_fun = opt.sk_fun;

% -------------------------
% 1) Geometry: lattice vectors and reciprocal vectors
% -------------------------
g = MTB.geometry("NGra");

a1 = [sqrt(3)/2, -1/2, 0.0] * a0;
a2 = [sqrt(3)/2,  1/2, 0.0] * a0;
a3 = [0.0, 0.0, c];

omega = dot(a1, cross(a2, a3));
b1 = 2*pi * cross(a2, a3) / omega;
b2 = 2*pi * cross(a3, a1) / omega;
b3 = 2*pi * cross(a1, a2) / omega;

g.a = [a1; a2; a3];
g.b = [b1; b2; b3];

% -------------------------
% 2) Build atoms (fractional) and wpos (Cartesian)
% -------------------------
g.atoms = [];

for i = 1:Nlay
    shift = (i-1)/3;
    zfrac = (i-1) * d0_frac;

    if opt.use_offset
        x1 = mod(opt.offset1 + shift, 1);
        y1 = mod(opt.offset1 + shift, 1);
        x2 = mod(opt.offset2 + shift, 1);
        y2 = mod(opt.offset2 + shift, 1);
    else
        x1 = mod(0.0 + shift, 1);
        y1 = mod(0.0 + shift, 1);
        x2 = mod(1/3 + shift, 1);
        y2 = mod(1/3 + shift, 1);
    end

    g.atoms = [g.atoms; ...
        x1, y1, zfrac; ...
        x2, y2, zfrac];
end

% Cartesian positions
g.wpos = g.atoms * g.a;

% center slab in the cell along z
g.wpos(:,3) = g.wpos(:,3) - mean(g.wpos(:,3)) + a3(3)/2;

% recompute fractional atoms from centered wpos (optional but keeps consistency)
g.atoms = g.wpos / g.a;

% -------------------------
% 3) Neighbor list (hopr) + allocate ham
% -------------------------
hopr = [];
result_matrices = cell(1,4);

% you used: (2, 5, 6) with SK; (21) with constant t2
result_matrices{1} = find_neighbor_data(g.wpos(1:end,:), g.a, nei_list(1), 3); % nn
hopr = [hopr; result_matrices{1}(:,3:5)];

result_matrices{2} = find_neighbor_data(g.wpos(1:end,:), g.a, nei_list(2), 3); % t1
hopr = [hopr; result_matrices{2}(:,3:5)];

result_matrices{3} = find_neighbor_data(g.wpos(1:end,:), g.a, nei_list(3), 3); % v3,v4
hopr = [hopr; result_matrices{3}(:,3:5)];

result_matrices{4} = find_neighbor_data(g.wpos(1:end,:), g.a, nei_list(4), 3); % t2
hopr = [hopr; result_matrices{4}(:,3:5)];

hopr = unique(hopr,'rows');

nbands = size(g.wpos, 1);
nrpts  = size(hopr, 1);
ham = zeros(nbands, nbands, nrpts);

% -------------------------
% 4) Fill ham for first three neighbor shells using Slater-Koster (sk_fun)
% -------------------------
for n_idx = 1:3
    result_matrix = result_matrices{n_idx};
    for ii = 1:size(result_matrix,1)
        [~, raw_index] = ismember(result_matrix(ii,3:5), hopr, "rows");
        orbital_1 = result_matrix(ii,1);
        orbital_2 = result_matrix(ii,2);

        % real-space vector from orbital_1 to orbital_2 (include lattice hop)
        R = result_matrix(ii,3:5) * g.a + g.wpos(orbital_2,:) - g.wpos(orbital_1,:);
        tsk = sk_fun(R);

        if ~isnan(tsk)
            ham(orbital_1, orbital_2, raw_index) = ham(orbital_1, orbital_2, raw_index) + tsk;
        end
    end
end

% -------------------------
% 5) Fill t2 (fourth shell) with constant value
% -------------------------
result_matrix = result_matrices{4};
for ii = 1:size(result_matrix,1)
    [~, raw_index] = ismember(result_matrix(ii,3:5), hopr, "rows");
    orbital_1 = result_matrix(ii,1);
    orbital_2 = result_matrix(ii,2);
    ham(orbital_1, orbital_2, raw_index) = ham(orbital_1, orbital_2, raw_index) + t2;
end

% -------------------------
% 6) Pack outputs
% -------------------------
g.hopr = hopr;
g.ham  = ham;

% keep as-is; your original line does nothing effectively
g.wpos = kron(g.wpos, ones(1,1));

end


function [w_layer, M] = get_layer_weight_fromUnk(Unk, nlayer, orbPerLayer)
%GET_LAYER_WEIGHT_FROMUNK  Compute layer weights from eigenvectors Unk.
%
%   [w_layer, M] = get_layer_weight_fromUnk(Unk, nlayer, orbPerLayer)
%
% INPUT
%   Unk          : (nband x nband x nk) eigenvector matrices along k.
%                  Convention: each column is an eigenvector of a band
%                  in the orbital basis.
%   nlayer       : number of layers (e.g., 7)
%   orbPerLayer  : number of orbitals per layer (default: 2)
%
% OUTPUT
%   w_layer      : (nlayer x nband x nk) layer weights
%                  w_layer(L, n, ik) = weight of band n on layer L at k-index ik
%   M            : (nlayer x nband) 0/1 layer mask matrix used in the projection
%
% NOTE
%   Requires that orbitals are ordered layer-by-layer, i.e.
%   layer 1 -> orbitals 1..orbPerLayer,
%   layer 2 -> orbitals orbPerLayer+1..2*orbPerLayer, etc.

    if nargin < 3 || isempty(orbPerLayer)
        orbPerLayer = 2;
    end

    % --- size checks
    [nb1, nb2, nk] = size(Unk);
    if nb1 ~= nb2
        error('Unk must be nband x nband x nk (square in first 2 dims).');
    end
    nband = nb1;

    if nband ~= nlayer * orbPerLayer
        error('Size mismatch: nband=%d but nlayer*orbPerLayer=%d.', ...
              nband, nlayer*orbPerLayer);
    end

    % --- build layer mask matrix M (nlayer x nband)
    M = zeros(nlayer, nband);
    for L = 1:nlayer
        idx = (orbPerLayer*(L-1)+1) : (orbPerLayer*L);
        M(L, idx) = 1;
    end

    % --- fast batch computation using reshape (no loop over k)
    A  = abs(Unk).^2;                 % nband x nband x nk
    A2 = reshape(A, nband, []);       % nband x (nband*nk)
    W2 = M * A2;                      % nlayer x (nband*nk)
    w_layer = reshape(W2, nlayer, nband, nk);

    % --- optional sanity check (comment out if you want absolute speed)
    % Each band should have total weight ~1 across layers
    % err = max(abs(sum(w_layer,1) - 1), [], 'all');
    % if err > 1e-6
    %     warning('Layer weights do not sum to 1 (max err = %.3e). Check Unk convention/normalization.', err);
    % end
end




function writeoutput(filename,list)
    file=fopen(filename,'a+');
    cloumns=size(list,2);
    raws=size(list,1);
    % fprintf(file,'Time on %s\n',datetime('today'));
    % fprintf(file,'raws %d cloums %d\n',raws,cloumns);
    for i=1:raws
        for j=1:cloumns
            fprintf(file,'%12.6f \t',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
    fclose(file);
end
