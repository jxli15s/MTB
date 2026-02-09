clear;
clear all;
%p=parpool(8)
%!!!!!!!!! note that this should be used at tb/matlab directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             Construct Band structure from DFT                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d/for_ming/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/for_ming/wannier90_hr_p1.dat','data/TaIrTe4_2d/for_ming/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;
orbital_num=[18,18,12,12,8,8,8,8,8,8,8,8];
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
        repmat(g.wpos(3,:),[orbital_num(3),1]);...
        repmat(g.wpos(4,:),[orbital_num(4),1]);...
        repmat(g.wpos(5,:),[orbital_num(5),1]);...
        repmat(g.wpos(6,:),[orbital_num(6),1]);...
        repmat(g.wpos(7,:),[orbital_num(7),1]);...
        repmat(g.wpos(8,:),[orbital_num(8),1]);...
        repmat(g.wpos(9,:),[orbital_num(9),1]);...
        repmat(g.wpos(10,:),[orbital_num(10),1]);...
        repmat(g.wpos(11,:),[orbital_num(11),1]);...
        repmat(g.wpos(12,:),[orbital_num(12),1]);...
    ];
g.wpos=g.wpos;

%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
labels={'X','Y','\Gamma','R'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

% labels={'R','Y','X'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.0,0.0]};% hkpoints-high symmetry k points

efermi=-0.4805;
nk=301;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.2,0.2])
xlim([0.6,1.3])

%% calculate Density of States
data=textread("/Volumes/T9/work/dft/TaIrTe4/1s/tairte_vasp/pbe/noivdw/new/wtool/bulkek_plane-matlab-del.dat");
Enk=squeeze(data(:,7:end));
plottap=1;
nk=301;
Enum=1000;
Emin=-0.2;
Emax=0.2;
eps=(Emax-Emin)/Enum*15;
Nband=size(Enk,2);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
Dos=Dos*(Emax-Emin)/Enum;
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16;
%%
figure()
hold on;
plot(Eaxis,Dos*5*10^16)
plot(Eaxis,abs(TDos_new))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Plot Band and Dos together                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = kindex;
linesize = 1;

Eband = Energy - efermi;   % band能量（相对费米能级）

% 你DOS的变量假设：
%   Eaxis: 能量轴（eV，和Energy同单位，建议也减去efermi）
%   Dos  : DOS(Eaxis)  (同长度向量)
Edos = Eaxis;

fig=figure('Color','white')
set(gcf,'Units','normalized','Position',[0.2 0.2 0.6 0.9]);  % 整个窗口大小（可选）

% ====== 你想控制的参数（改这里） ======
left   = 0.10;   % 左边距
right  = 0.03;   % 右边距
bottom = 0.12;   % 下边距
top    = 0.06;   % 上边距
gap    = 0.015;  % 两子图间距（越小越紧）

w1 = 0.8;       % 左图占可用宽度比例（0~1）
% 右图自动 = 1 - w1
% ======================================

W = 1 - left - right - gap;    % 可用总宽
H = 1 - bottom - top;          % 可用总高

ax1 = axes('Units','normalized', ...
           'Position',[left, bottom, w1*W, H]);

ax2 = axes('Units','normalized', ...
           'Position',[left + w1*W + gap, bottom, (1-w1)*W, H]);


% ---- 左：Band
axes(ax1); hold on

for i = 1:size(Eband,1)
    plot(ax1, kpath, Eband(i,:), 'k', 'LineWidth', linesize);
end

% y=0 费米能级
plot(ax1, kpath, zeros(1,length(kpath)), '--k', 'LineWidth', 1);

% 高对称点竖线
for i = 1:length(kk)-2
    x = kk(i+1);
    plot(ax1, [x x], [min(Eband(:))-1, max(Eband(:))+1], '--k', 'LineWidth', linesize);
end

box(ax1,'on'); grid(ax1,'off')
xlim(ax1, [0.6, 1.3])
xticks(ax1, kk)
xticklabels(ax1, labels)
ylim(ax1, [-0.2 0.2])
yticks(-0.2:0.1:0.2)

ax1.LineWidth = 1;
ax1.XAxis.FontSize = 18;
ax1.YAxis.FontSize = 18;
ylabel(ax1,'E - E_F (eV)')

% ---- 右：DOS
axes(ax2); hold on
% Dos_s = smoothdata(Dos,'movmean',100);
plot(ax2, Dos/max(Dos), Edos, 'k', 'LineWidth', 1);   % 注意是 (Dos, E)
yline(ax2, 0, '--k', 'LineWidth', 1);          % 对齐E=0
box(ax2,'on'); grid(ax2,'off')

% 与左图共享能量范围（关键）
ylim(ax2, ylim(ax1))

% 让右图的能量刻度不重复占空间（可选）
ax2.YTickLabel = [];
xlabel(ax2,'DOS ($\times 10^{14} eV^{-1}\cdot cm^{-1}$)','Interpreter','latex','FontSize',17)
xlim([0,1.2])

% ax2 = nexttile(2);
% ax2.YTick = [];                 % 或者 ax2.YTickLabel = [];
% ax2.YColor = 'none';            % 彻底不画右边y轴

% ax2.YTickLabel = [];
% ax2.YLabel.String = '';

ax2.LineWidth = 1;

set(findall(gcf,'-property','FontName'), 'FontName','Arial','FontSize',28)


% exportgraphics(fig, 'tit_QMD/band_dos_fplo.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Plot the 3D band structure                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=textread("/Volumes/T9/work/dft/TaIrTe4/1s/tairte_vasp/pbe/noivdw/new/wtool/bulkek_plane-matlab-del.dat");
kx=reshape(data(:,4),301,301);
ky=reshape(data(:,5),301,301);
band_v=reshape(data(:,8),301,301);
band_v= 0.5 * (band_v + fliplr(band_v));

band_c=reshape(data(:,9),301,301);
band_c= 0.5 * (band_c + fliplr(band_c));

%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx, ky: 分别是 nkx × 1 和 nky × 1 的向量
% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序

       
% 绘图
figure('Color','white');
hold on;

% 1. 三维能带
% s = surf(kx, ky, Ez, 'EdgeColor', 'none');
% colormap turbo
% shading interp
% alpha(0.9);
Ez1 = band_c;
surf(kx, ky, Ez1, 'EdgeColor', 'none');
Ez2 = band_v;
% surf(kx, ky, Ez2, 'EdgeColor', 'none');


colormap(slanCM('RdBu'))
colormap(flipud(colormap));

Ef_list1 =[0.0502 0.0502]  % 0.0137  cb  37.5meV vhS electron side
Ef_list2=[-0.05928 -0.05928] % -0.0093 vt  50meV vhS hole side
% 2. 费米面等高线（Ef）
contour3(kx, ky, Ez1, Ef_list1, 'LineColor', '[0.3,0.5,0.7]', 'LineWidth', 2);
% contour3(kx, ky, Ez2, Ef_list2, 'LineColor', 'blue' ,'LineWidth', 2);
Ef_list = -0.2:0.031:0.2;  % 多个等高值
% contour3(kx, ky, Ez2, Ef_list, 'LineColor', '[0.5,0.5,0.5]' ,'LineWidth', 0.5);
Ef_list = -0.2:0.033:0.2;  % 多个等高值
contour3(kx, ky, Ez1, Ef_list, 'LineColor', '[0.5,0.5,0.5]' ,'LineWidth', 0.5);
% 3. 视图与标签

view(3);
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
% title('Bands and Fermi Surface ');
% colorbar;
axis tight;
% axis equal;
box off
axis off
xlim([-0.3,0.3])
zlim([-0.2,0.2])
clim([-0.3,0.3])
% camlight('headlight')
camlight(45, 30)      % 方位角/俯仰角
% camlight(-60, 20)
% 灯光：至少两个光源会更像"环境光"
camlight('headlight')                      % 跟相机走
camlight('right')                          % 侧面补光
camlight('left')

lighting phong
shading interp    % 表面平滑
material dull
camproj perspective
set(gcf,'Renderer','opengl')               % 很关键：用 OpenGL 才有光照/高光
% view(-6,15)


view(-6,32)

% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
%export svg
% exportgraphics(gca, 'tit_QMD/3D_band_fplo.svg', ...
    % 'BackgroundColor','none', 'ContentType','vector', 'Resolution', 600);
%export pdf
% exportgraphics(gca, 'tit_QMD/3D_band_fplo.pdf', ...
%     'BackgroundColor','none', 'ContentType','vector');
% export png
exportgraphics(gca, 'tit_QMD/3D_band_fplo_cband.png', ...
    'BackgroundColor','none', 'Resolution', 1200);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% load('xinitial_good.mat')
% Step 1: Initialize Geometry and Hamiltonian
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);

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

% Step 5: Neighbor Search
[result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs);
%
% Step 6: Initialize fixed Hartree-Fock States parameters
nec1 = (4 * n1 + 4) / n1 / 2 / 2 / 2;
nec2 = (4 * n1) / n1 / 2 / 2 / 2;
% xinitial_0 = diag(kron(ones(1, n1), [nec1 - 0.4, nec1 - 0.4, nec2 + 0.4, nec2 + 0.4, ...
%                                     nec1 + 0.4, nec1 + 0.4, nec2 - 0.4, nec2 - 0.4]));
xinitial_0 = diag(kron(ones(1, n1*n2), rand(1,8)));

[nbands,~,nrpts]=size(gs.ham);
% labels={'Y','\Gamma','X','R'}; % labels for k
% hkpoints={[0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'R','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Read and Write the wannier90_hr.dat for TI    %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=load("data/tit_hf/15x1/10nmd/U-epsilon-f/"+int2str(12)+"-15x1"+"/result_s1_eps"+int2str(13)+".00_r"+int2str(11)+".mat");
U0=result.U0;
U=result.V;
V=result.V;
xinitial=result.xinitial;
pairsU=result.pairsU;
pairsV=result.pairsV;
efermi=result.efermi;
modifyHam(gs, xinitial, V, V, pairsU, pairsV)
% Add Zeeman and Electric Field
%
Electric_field_in_evpA=0.1;
gs=add_elec(gs,Electric_field_in_evpA);

s3=[1  0
    0  -1];
Zeeman=kron(eye(60),0.001*s3);
gs.add_zeeman(Zeeman)
% plot the band along HSL
Electric_field_in_evpA=0.00*0.529177; nk=251;
%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1+1,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
ylim([-0.1,0.1])
%%
% Calculate the Ham on K-mesh plane
knum  = 200;
kxline = [-0.5,0.5];
kyline = [-0.5,0.5];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]  = MTB.ham.get_bulk_plane_bands(g, Kx,Ky,Kz);
[Hamk,Unk,Enk]  = MTB.ham.get_bulk_plane_bands_with_Ham(gs, Kx,Ky,Kz);
efermi=calculate_ef(Enk(:),0.5);
Enk=Enk-efermi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the DOS and TDOS for the carrier density  %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plottap=2;
nk=knum;
Enum=3000;
Emin=-0.1;
Emax=0.1;
eps=(Emax-Emin)/Enum*10;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
Dos=Dos*(Emax-Emin)/Enum;
TDos_new=TDos/norm(cross(gs.a(1,:),gs.a(2,:)))*10^16;
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Plot Band and Dos together                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = kindex;
linesize = 1;

Eband = Energy - efermi-0.008;   % band能量（相对费米能级）

% 你DOS的变量假设：
%   Eaxis: 能量轴（eV，和Energy同单位，建议也减去efermi）
%   Dos  : DOS(Eaxis)  (同长度向量)
Edos = Eaxis;

fig=figure('Color','white')
set(gcf,'Units','normalized','Position',[0.2 0.2 0.6 0.9]);  % 整个窗口大小（可选）

% ====== 你想控制的参数（改这里） ======
left   = 0.10;   % 左边距
right  = 0.03;   % 右边距
bottom = 0.12;   % 下边距
top    = 0.06;   % 上边距
gap    = 0.015;  % 两子图间距（越小越紧）

w1 = 0.8;       % 左图占可用宽度比例（0~1）
% 右图自动 = 1 - w1
% ======================================

W = 1 - left - right - gap;    % 可用总宽
H = 1 - bottom - top;          % 可用总高

ax1 = axes('Units','normalized', ...
           'Position',[left, bottom, w1*W, H]);

ax2 = axes('Units','normalized', ...
           'Position',[left + w1*W + gap, bottom, (1-w1)*W, H]);


% ---- 左：Band
axes(ax1); hold on

for i = 1:size(Eband,1)
    plot(ax1, kpath, Eband(i,:), 'k', 'LineWidth', linesize);
end

% y=0 费米能级
plot(ax1, kpath, zeros(1,length(kpath)), '--k', 'LineWidth', 1);

% 高对称点竖线
for i = 1:length(kk)-2
    x = kk(i+1);
    plot(ax1, [x x], [min(Eband(:))-1, max(Eband(:))+1], '--k', 'LineWidth', linesize);
end

box(ax1,'on'); grid(ax1,'off')
xlim(ax1, [0, kk(end)])
xticks(ax1, kk)
xticklabels(ax1, labels)
ylim(ax1, [-0.08 0.08])
yticks(-0.08:0.04:0.08)

ax1.LineWidth = 1;
ax1.XAxis.FontSize = 18;
ax1.YAxis.FontSize = 18;
ylabel(ax1,'E - E_F (eV)')

% ---- 右：DOS
axes(ax2); hold on
% Dos_s = smoothdata(Dos,'movmean',100);
plot(ax2, Dos/max(Dos), Edos-0.008, 'k', 'LineWidth', 1);   % 注意是 (Dos, E)
yline(ax2, 0, '--k', 'LineWidth', 1);          % 对齐E=0
box(ax2,'on'); grid(ax2,'off')

% 与左图共享能量范围（关键）
ylim(ax2, ylim(ax1))

% 让右图的能量刻度不重复占空间（可选）
ax2.YTickLabel = [];
xlabel(ax2,'DOS ($\times 10^{14} eV^{-1}\cdot cm^{-1}$)','Interpreter','latex','FontSize',17)
xlim([0,1.2])

% ax2 = nexttile(2);
% ax2.YTick = [];                 % 或者 ax2.YTickLabel = [];
% ax2.YColor = 'none';            % 彻底不画右边y轴

% ax2.YTickLabel = [];
% ax2.YLabel.String = '';

ax2.LineWidth = 1;
axis normal
set(findall(gcf,'-property','FontName'), 'FontName','Arial','FontSize',28)

exportgraphics(fig, 'tit_QMD/band_dos_hfmf.svg', ...
    'BackgroundColor','none', 'ContentType','vector');

%%
% 绘图
band_v=Enk(:,:,60);
band_c=Enk(:,:,61);
figure('Color','white');
hold on;

% 1. 三维能带
% s = surf(kx, ky, Ez, 'EdgeColor', 'none');
% colormap turbo
% shading interp
% alpha(0.9);
Ez1 = band_c;
surf(Kx, Ky, Ez1, 'EdgeColor', 'none');
Ez2 = band_v;
surf(Kx, Ky, Ez2, 'EdgeColor', 'none');

colormap(slanCM('RdBu'))
colormap(flipud(colormap));
view(3);
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
% title('Bands and Fermi Surface ');
colorbar;
% axis tight;
axis equal;
box off
axis off
% xlim([-0.3,0.3])
% zlim([-0.2,0.2])
clim([-0.06,0.06])
% camlight('headlight')
% camlight(45, 30)      % 方位角/俯仰角
% camlight(-60, 20)
% 灯光：至少两个光源会更像"环境光"
camlight('headlight')                      % 跟相机走
% camlight('right')                          % 侧面补光
lighting phong
shading interp    % 表面平滑
material dull
camproj perspective
set(gcf,'Renderer','opengl')               % 很关键：用 OpenGL 才有光照/高光
% view(-25,8)
view(-90,12)

% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% % export png
% exportgraphics(gca, 'tit_QMD/3D_band_hfmf_side.png', ...
%     'BackgroundColor','none', 'Resolution', 600);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       QMD for different B                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siamg_abc_all=zeros(2,2,2,3000,10);
for i =1:9
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_0"+int2str(i)+"meV_30K.mat";
load(filename)
sigma_abc_all(:,:,:,:,i)=sigma_abc;
end
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat";
load(filename)
sigma_abc_all(:,:,:,:,10)=sigma_abc;
%%
fig=figure('Color','white')
abc=[2,1,1];
a=abc(1);b=abc(2);c=abc(3);
peak_1=[]
% figure()
hold on

N = 10;
cmap = [linspace(0,1,N)', zeros(N,1), linspace(1,0,N)'];  % [R,G,B]：红->蓝
plot(Ef_list-0.008,zeros(size(y)),'LineWidth',1.5,'Color',cmap(1,:))
for i=1:N
    y=1000*squeeze(sigma_abc_all(a,b,c,:,i))*pi;
    peak_1=[peak_1,y(1755)]; % 0.009meV
    plot(Ef_list-0.008,y,'LineWidth',1.5,'Color',cmap(i,:))
end

T_list=[0,1,2,3,4,5,6,7,8,9,10]
legend(arrayfun(@(T)sprintf('%g meV',T/10), T_list, 'UniformOutput', false), ...
       'Location','best');

% legend('0 meV',  'Location','best');

box on;
xlim([-0.03,0.03])
% axis normal
% axis equal
ylim([-1.5,1.5])
daspect([1 20 1])   % x:y:z 的数据单位比例（2D时等同于x:y=1:1）
xlabel('E-E$_f$ (eV)','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm/V$^2$','Interpreter','latex','FontSize',18)

exportgraphics(fig, 'tit_QMD/QMD_different_B_baa_longy.svg', ...
    'BackgroundColor','none', 'ContentType','vector');
% exportgraphics(fig, 'tit_QMD/QMD_different_B_abb.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
fig=figure('Color','white');
B_list=0.1:0.1:1;
plot(B_list,peak_1-0.05,'ro-','LineWidth',1)
xlabel('$B_{eff}$ (meV)','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm/V$^2$','Interpreter','latex','FontSize',18)
ylim([0,1.4])
% xlim([0,1.1])

exportgraphics(fig, 'tit_QMD/QMD_different_B_baa_peak.svg', ...
    'BackgroundColor','none', 'ContentType','vector');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              QMD for different B 1meV 30K for baa                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat")
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            % plot(Ef_list-0.008,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            plot(TDos_new*1.8,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            % filename="sigma_abc_"+int2str(x)+int2str(y)+int2str(z)+".dat";
            % outlist=[Ef_list.',reshape(sigma_abc(x,y,z,:),300,1)];
            % writeoutput(filename,outlist)
        end
    end
end
plot(Ef_list-0.008,1000*squeeze(sigma_abc(1,2,2,:))*pi)

% plot(Ef_list,1000*sigma_yxx,'ko')
ylabel('mA/V^2')
xlabel('Energy(eV)')
%%
x = Ef_list - 0.008;
% y = 1000*squeeze(sigma_abc(1,2,2,:))*pi-1.5e-4;
y = 1000*squeeze(sigma_abc(2,1,1,:))*pi;

yp = max(y,0);
yn = min(y,0);

fig=figure('Color','white'); hold on
hp = area(x, yp, 0);  hp.FaceColor = [165/255 48/255 50/255]; hp.EdgeColor='none'; hp.FaceAlpha=0.4;
hn = area(x, yn, 0);  hn.FaceColor = [50/255 111/255 160/255]; hn.EdgeColor='none'; hn.FaceAlpha=0.4;
plot(x,y,'k','LineWidth',1.5); yline(0,'k-');
box on;
axis normal
xlabel('E-E$_f$ (eV)','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm/V$^2$','Interpreter','latex','FontSize',18)
xlim([-0.03,0.03])
% ylim([-1,1.6])
% 
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_baa.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
figure()
hold on;
% plot(Ef_list-0.008,1000*squeeze(sigma_abc(2,1,1,:))*pi) %1496-1741
end1=1450;
% start2=1741;
start2=1451;
TDos_new_del=[TDos_new(1:end1),TDos_new(start2:end)]*1.8/10^12;
sigma_abc_p=1000*squeeze(sigma_abc(2,1,1,:))*pi;
sigma_abc_p=sigma_abc_p';
sigma_abc_del=[sigma_abc_p(1:end1),sigma_abc_p(start2:end)];
% plot(TDos_new,1000*squeeze(sigma_abc(2,1,1,:))*pi)
% plot(TDos_new_del,sigma_abc_del)
plot(TDos_new_del,sigma_abc_del,'k','LineWidth',1.5);
x=TDos_new_del;
y=sigma_abc_del;
yp = max(y,0);
yn = min(y,0);
%%
fig=figure('Color','white'); hold on
hp = area(x, yp, 0);  hp.FaceColor = [165/255 48/255 50/255]; hp.EdgeColor='none'; hp.FaceAlpha=0.4;
hn = area(x, yn, 0);  hn.FaceColor = [50/255 111/255 160/255]; hn.EdgeColor='none'; hn.FaceAlpha=0.4;
plot(TDos_new_del,sigma_abc_del,'k','LineWidth',1);
% yline(0,'k-');
xlabel('n $(10^{12}/cm^{2})$','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm / V$^2$','Interpreter','latex','FontSize',18)
box on;
axis normal
xlim([-8,8])
ylim([-1,1.5])
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_baa_density_n.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              QMD for different B 1meV 30K  for abb                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat")
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            % plot(Ef_list-0.008,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            plot(TDos_new*1.8,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            % filename="sigma_abc_"+int2str(x)+int2str(y)+int2str(z)+".dat";
            % outlist=[Ef_list.',reshape(sigma_abc(x,y,z,:),300,1)];
            % writeoutput(filename,outlist)
        end
    end
end
plot(Ef_list-0.008,1000*squeeze(sigma_abc(1,2,2,:))*pi)

% plot(Ef_list,1000*sigma_yxx,'ko')
ylabel('mA/V^2')
xlabel('Energy(eV)')
%% 
x = Ef_list - 0.008;
y = 1000*squeeze(sigma_abc(1,2,2,:))*pi-1.2e-4;
% y = 1000*squeeze(sigma_abc(2,1,1,:))*pi;

yp = max(y,0);
yn = min(y,0);

fig=figure('Color','white'); hold on
hp = area(x, yp, 0);  hp.FaceColor = [165/255 48/255 50/255]; hp.EdgeColor='none'; hp.FaceAlpha=0.4;
hn = area(x, yn, 0);  hn.FaceColor = [50/255 111/255 160/255]; hn.EdgeColor='none'; hn.FaceAlpha=0.4;
plot(x,y,'k','LineWidth',1.5); yline(0,'k-');
box on;
axis normal
xlabel('E-E$_f$ (eV)','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm/V$^2$','Interpreter','latex','FontSize',18)
xlim([-0.03,0.03])
% ylim([-1,1.6])
% 
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_abb.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
figure()
hold on;
% plot(Ef_list-0.008,1000*squeeze(sigma_abc(2,1,1,:))*pi) %1496-1741
end1=1450;
% start2=1741;
start2=1451;
TDos_new_del=[TDos_new(1:end1),TDos_new(start2:end)]*1/10^12;
sigma_abc_p=1000*squeeze(sigma_abc(1,2,2,:))*pi;
sigma_abc_p=sigma_abc_p';
sigma_abc_del=[sigma_abc_p(1:end1),sigma_abc_p(start2:end)];
% plot(TDos_new,1000*squeeze(sigma_abc(2,1,1,:))*pi)
% plot(TDos_new_del,sigma_abc_del)
plot(TDos_new_del,sigma_abc_del,'k','LineWidth',1.5);
x=TDos_new_del;
y=sigma_abc_del;
yp = max(y,0);
yn = min(y,0);
%%
fig=figure('Color','white'); hold on
hp = area(x, yp, 0);  hp.FaceColor = [165/255 48/255 50/255]; hp.EdgeColor='none'; hp.FaceAlpha=0.4;
hn = area(x, yn, 0);  hn.FaceColor = [50/255 111/255 160/255]; hn.EdgeColor='none'; hn.FaceAlpha=0.4;
plot(TDos_new_del,sigma_abc_del,'k','LineWidth',1);
% yline(0,'k-');
xlabel('n $(10^{12}/cm^{2})$','Interpreter','latex','FontSize',18)
ylabel('mA$\cdot$nm / V$^2$','Interpreter','latex','FontSize',18)
box on;
axis normal
xlim([-8,8])
ylim([-0.03,0.03])
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_abb_density_n.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        QMD for QMD and QM along the high symmetry line            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K_band.mat')
load('/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat')
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            % plot(Ef_list,1000*squeeze(sigma_abc(2,1,1,:)*pi))
            plot(TDos_new,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            % filename="sigma_abc_"+int2str(x)+int2str(y)+int2str(z)+".dat";
            % outlist=[Ef_list.',reshape(sigma_abc(x,y,z,:),300,1)];
            % writeoutput(filename,outlist)
        end
    end
end
% plot(Ef_list,1000*sigma_yxx,'ko')
ylabel('mA/V^2')
% xlabel('Energy(eV)')
xlabel('Carrier Density n (1/cm^2)')
% xlim([-0.1 0.1])
%%
fig=figure('Color','white'); hold on
plot(Ef_list-0.008,1000*squeeze(sigma_abc(2,1,1,:))*pi,'DisplayName','D_{yxx}')
plot(Ef_list-0.008,1000*squeeze(sigma_abc(1,2,1,:))*pi,'DisplayName','D_{xyx}')
plot(Ef_list-0.008,1000*squeeze(sigma_abc(1,2,2,:))*pi*2,'DisplayName','D_{xyy}')
plot(Ef_list-0.008,1000*squeeze(sigma_abc(2,1,2,:))*pi*2,'DisplayName','D_{yxy}')
box on;
legend('Location','best','Interpreter','tex') 
xlim([-0.03,0.03])
ylabel('$\sigma^{2\omega}$(mA$\cdot$nm / $V^2$) ','Interpreter','latex')
% ylabel('mA\cdot nm / V^2')
xlabel('E-Ef (eV)')
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
axis normal
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_abb_baa_sigma_lines.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
% legend(arrayfun(@(T)sprintf('%g meV',T/10), T_list, 'UniformOutput', false), ...
%        'Location','best');

%%

theta = linspace(0,2*pi,721);   % 0..360 deg
E0 = 1;                         % 你也可以取 rms 或 1，反正只差整体系数

s_yxx = 1.1;    % sigma_yxx
s_xyy = -0.03;    % sigma_xyy

Jx = s_xyy * sin(theta).^2.*sin(theta);
Jy = s_yxx * cos(theta).^2.*cos(theta);

fig=figure('Color','white'); hold on
plot(theta*180/pi, Jx+Jy, 'LineWidth', 1.5)
% plot(theta*180/pi, Jy, 'LineWidth', 1.5)
% plot(theta*180/pi, Jpar,'LineWidth', 1.5)
yline(0,'k--','LineWidth',1.2)
xlabel('\theta (deg)'); xlim([0 360])
ylabel('mA\cdot nm / V^2')
box on;
axis normal;
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
clim([-5000,5000])
axis normal

exportgraphics(fig, 'tit_QMD/sigma_theta.svg', ...
    'BackgroundColor','none', 'ContentType','vector');
% figure()
% polarplot(theta,Jx+Jy)

%%
abc=[2,1,1];
a=abc(1);b=abc(2);c=abc(3);
knum=500;
band1=zeros(knum/2,4);
banddk1=zeros(knum/2,4);
bandgxx1=zeros(knum/2,4);
bandgyx1=zeros(knum/2,4);

idx=0;
for i=knum:-1:knum/2+1
    idx=idx+1;
    for j=59:62
        band1(idx,j-58)=Enk(i,i,j);
        banddk1(idx,j-58)=Dk_core(i,i,a,b,c,j-54);
        bandgxx1(idx,j-58)=gk(i,i,1,1,j-54);
        bandgyx1(idx,j-58)=gk(i,i,2,1,j-54);
    end
end

band2=zeros(knum/2,4);
banddk2=zeros(knum/2,4);
bandgxx2=zeros(knum/2,4);
bandgyx2=zeros(knum/2,4);
idx=0
for i=knum/2+1:knum
    idx=idx+1;
    for j=59:62
        band2(idx,j-58)=Enk(knum/2+1,i,j);
        banddk2(idx,j-58)=Dk_core(knum/2+1,i,a,b,c,j-54);
        bandgxx2(idx,j-58)=gk(knum/2+1,i,1,1,j-54);
        bandgyx2(idx,j-58)=gk(knum/2+1,i,2,1,j-54);
    end
end

band3=zeros(knum/2,4);
banddk3=zeros(knum/2,4);
bandgxx3=zeros(knum/2,4);
bandgyx3=zeros(knum/2,4);
idx=0
for i=knum:-1:knum/2+1
    idx=idx+1
    for j=59:62
        band3(idx,j-58)=Enk(knum/2*3+1-i,i,j);
        banddk3(idx,j-58)=Dk_core(knum/2*3+1-i,i,a,b,c,j-54);
        bandgxx3(idx,j-58)=gk(knum/2*3+1-i,i,1,1,j-54);
        bandgyx3(idx,j-58)=gk(knum/2*3+1-i,i,2,1,j-54);
    end
end

band=[band1;band2;band3];
band=band-0.008;
banddk=[banddk1;banddk2;banddk3];
bandgxx=[bandgxx1;bandgxx2;bandgxx3];
bandgyx=[bandgyx1;bandgyx2;bandgyx3];
%%
%for gxx
fig=figure('Color','white'); hold on
kdist=kpath(1:3*knum/2);
msize=20;
Dmax=max(abs(banddk),[],'all');
hold on;
for i=1:4
% scatter(kdist,band(:,i),msize,banddk(:,i)/Dmax,'filled')
scatter(kdist,band(:,i),msize,bandgxx(:,i),'filled')
% plot(squeeze(Enk(100,:,i+54)))
end

kk=kindex;
linesize=1;
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',1)

for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[-0.1 0.1],'k--','LineWidth',linesize)
end

for i=1:4
plot(kdist,band(:,i),'--','LineWidth',1,'Color',[0.65 0.65 0.65])
end


grid off
box on
colorbar
cmax = max(bandgxx(:));             % 或者你自己指定
% cmax=1;


% 白 -> 红 colormap
% n = 256;
% cmap = [ linspace(1,0,n)', linspace(1,0,n)',ones(n,1)];  % [R,G,B]
% colormap(cmap)
% clim([0 cmax])

% colormap(slanCM('RdBu'))
% colormap(slanCM('red'))
% colormap(flipud(colormap));

shading interp

colormap(slanCM('RdBu'))
colormap(flipud(colormap));
xlim([0 kpath(end-1)])
ylim([-0.075,0.06])
xticks(kk)
xticklabels(labels)
% FontSize
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
clim([-9000,9000])
axis normal

% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_band_gxx.pdf', ...
%     'BackgroundColor','none', 'ContentType','vector');
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_band_gxx.png', ...
%     'BackgroundColor','none', 'Resolution', 600);
%%
%for gyx
fig=figure('Color','white'); hold on
kdist=kpath(1:3*knum/2);
msize=20;
Dmax=max(abs(banddk),[],'all');
hold on;
for i=1:4
% scatter(kdist,band(:,i),msize,banddk(:,i)/Dmax,'filled')
scatter(kdist,band(:,i),msize,bandgyx(:,i),'filled')
end

kk=kindex;
linesize=1;
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',1)

for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[-0.1 0.1],'k--','LineWidth',linesize)
end

for i=1:4
plot(kdist,band(:,i),'-','LineWidth',1,'Color',[0.65 0.65 0.65])
end

grid off
box on

shading interp

colormap(slanCM('RdBu'))
colormap(flipud(colormap));
xlim([0 kpath(end-1)])
ylim([-0.075,0.06])
xticks(kk)
xticklabels(labels)
% FontSize
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
colorbar 
clim([-60,60])
axis normal

% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_band_gyx.pdf', ...
%     'BackgroundColor','none', 'ContentType','vector');
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_band_gyx.png', ...
%     'BackgroundColor','none', 'Resolution', 600);
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_band_gyx.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
%for Dcore
fig=figure('Color','white'); hold on
kdist=kpath(1:3*knum/2);
msize=20;
Dmax=max(abs(banddk),[],'all');
hold on;
for i=1:4
scatter(kdist,band(:,i),msize,banddk(:,i)/Dmax,'filled')
% scatter(kdist,band(:,i),msize,bandgxx(:,i),'filled')
% plot(squeeze(Enk(100,:,i+54)))
end

kk=kindex;
linesize=1;
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',1)

for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[-0.1 0.1],'k--','LineWidth',linesize)
end

for i=1:4
plot(kdist,band(:,i),'--','LineWidth',1,'Color',[0.65 0.65 0.65])
end


grid off
box on
colorbar

shading interp
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
xlim([0 kpath(end-1)])
ylim([-0.075,0.06])
xticks(kk)
xticklabels(labels)
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
clim([-0.5,0.5])
axis normal
% set(gcf,'Color','none');      % figure 背景透明
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(fig, 'tit_QMD/QMD_B1meV_band_Dkbaa.svg', ...
%     'BackgroundColor','none', 'ContentType','vector');
%%
dk121=squeeze(Dk_core(:,:,2,1,1,6));
% dk121=sign(dk121).*log(abs(dk121)+100)
s = prctile(abs(dk121(:)), 80);   % 也可 95/99
dk121 = asinh(dk121 / s);
% dk121=dk121(1:71,31:71)
% dk121=squeeze(Dk_core(:,:,2,1,1,8));
 figure()
 surf(Kx,Ky,dk121,'EdgeColor','none')

%  dx = max(Kx(:)) - min(Kx(:));
% 
% Kx3 = [Kx-dx, Kx, Kx+dx];
% Ky3 = [Ky,    Ky, Ky   ];
% Z3  = [dk121, dk121, dk121];

% figure
% surf(Kx3, Ky3, Z3, 'EdgeColor','none');
% shading interp
view(2)
xlim([-0.0554,0.0554])
 shading interp
 colormap(flipud(colormap));
colormap(slanCM('RdBu'))
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
axis normal
view(2)
%%

 Ef=-0.01+0.008;
 eta       = 3e-4;
 band_list=55:66;
 Enk_sel = Enk(:,:,band_list);
 Nb_sel=length(band_list);
 delta_n = (1/pi)*eta ./ ((Enk_sel - Ef).^2 + eta^2); % Nkx x Nky x Nb_sel
 [Nkx,Nky]=size(Kx);
 Dmap = zeros(Nkx,Nky);
 for in=1:Nb_sel
     Dmap = Dmap + squeeze(Dk_core(:,:,2,1,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
     % Dmap = Dmap + squeeze(gk(:,:,2,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
 end


 % clim([-4e18,4e18]) %for  Ef=0.008+0.008;
 % clim([-1.5e19,1.5e19])   %for Ef=-0.01+0.008;
% cmax = 4e18;
cmax = 2e19;


%  thr = 0.0000000000000001*cmax;     % 例如把 2% cmax 内的都当作0
% Dmap2 = Dmap;
% Dmap2(abs(Dmap2) < thr) = 0;
% Dmap=Dmap2;

 D_yxx = sum(Dmap,'all')



 fig=figure('Color','white'); hold on
 surf(Kx,Ky,Dmap,'EdgeColor','none')

 colormap(slanCM('RdBu'))
 colormap(flipud(colormap));
colorbar

clim([-cmax,cmax])
cmap = flipud(slanCM('RdBu'));
n = size(cmap,1);
frac = 0.0006;
mid = (n+1)/2; w = round(frac*n);
cmap(max(1,floor(mid-w)):min(n,ceil(mid+w)),:) = 1;
colormap(cmap);
colorbar


 xlim([min(Kx,[],'all'),max(Kx,[],'all')])
 ylim([min(Ky,[],'all'),max(Ky,[],'all')])


 box on;
 view(2)
 xlabel('$K_x$','Interpreter','latex')
 ylabel('$K_y$','Interpreter','latex')
 ax=gca;
 ax.LineWidth=1;
 shading interp
 ax.YAxis.FontSize=18;
 ax.XAxis.FontSize=18;
 axis normal
 % set(gcf,'Color','none');      % figure 背景透明
 % set(gca,'Color','none');      % axes 背景透明
 % exportgraphics(gca, 'tit_QMD/QMD_B1meV_Dkbaa_plane.png', ...
 %     'BackgroundColor','none', 'Resolution', 600);
 
 %  set(gcf,'Color','none');      % figure 背景透明
 % set(gca,'Color','none');      % axes 背景透明
 % exportgraphics(gca, 'tit_QMD/QMD_B1meV_Dkbaa_plane_white.png', ...
 %     'BackgroundColor','none', 'Resolution', 1200);

% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_Dkbaa_plane_negef01.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_Dkbaa_plane_negef01_white.png', ...
%      'BackgroundColor','none', 'Resolution', 1200);

 %%
 Ef=0.008+0.008;
 eta       = 5e-4;
 band_list=55:66;
 Enk_sel = Enk(:,:,band_list);
 Nb_sel=length(band_list);
 delta_n = (1/pi)*eta ./ ((Enk_sel - Ef).^2 + eta^2); % Nkx x Nky x Nb_sel

 %  kB_eV = 8.617333262e-5;   % eV/K
 %  T_K=30;
 % kBT = kB_eV * T_K;
 % x   = (Enk_sel - Ef) ./ (2*kBT);
 % delta_n = (1./(4*kBT)) ./ cosh(x).^2;



 [Nkx,Nky]=size(Kx);
 Dmap = zeros(Nkx,Nky);
 for in=1:Nb_sel
     Dmap = Dmap + squeeze(gk(:,:,1,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
 end

%  thr = 0.002*cmax;     % 例如把 2% cmax 内的都当作0
% Dmap2 = Dmap;
% Dmap2(abs(Dmap2) < thr) = 0;
% 
% Dmap=Dmap2;

 D_yxx = sum(Dmap,'all')

 fig=figure('Color','white'); hold on
 surf(Kx,Ky,Dmap,'EdgeColor','none')
 colorbar
 colormap(slanCM('RdBu'))
 colormap(flipud(colormap));

% clim([-0.5e6,0.5e6]) %  for Ef=0.008+0.008;
% clim([-3e5,3e5]) %  for Ef=-0.01+0.008;

cmax =3e5;
clim([-cmax,cmax])

cmap = flipud(slanCM('RdBu'));
n = size(cmap,1);
frac = 0.01;
mid = (n+1)/2; w = round(frac*n);
cmap(max(1,floor(mid-w)):min(n,ceil(mid+w)),:) = 1;
colormap(cmap);
colorbar


 xlim([min(Kx,[],'all'),max(Kx,[],'all')])
 ylim([min(Ky,[],'all'),max(Ky,[],'all')])
 % clim([-4e18,4e18])
 box on;
 view(2)
 xlabel('$K_x$','Interpreter','latex')
 ylabel('$K_y$','Interpreter','latex')
 ax=gca;
 ax.LineWidth=1;
 shading interp
 ax.YAxis.FontSize=18;
 ax.XAxis.FontSize=18;

axis normal
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gxx_plane.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gxx_plane_white.png', ...
%      'BackgroundColor','none', 'Resolution', 1200);

% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gxx_plane_negef01.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
% % set(gcf,'Color','none');      % figure 背景透明 
% % set(gca,'Color','none');      % axes 背景透明
% % exportgraphics(gca, 'tit_QMD/QMD_B1meV_gxx_plane_negef01_white.png', ...
% %      'BackgroundColor','none', 'Resolution', 1200);


%%
 Ef=0.008+0.008;
 eta       = 5e-4;
 band_list=55:66;
 Enk_sel = Enk(:,:,band_list);
 Nb_sel=length(band_list);
 delta_n = (1/pi)*eta ./ ((Enk_sel - Ef).^2 + eta^2); % Nkx x Nky x Nb_sel
 [Nkx,Nky]=size(Kx);
Dmap = zeros(Nkx,Nky);
for in=1:Nb_sel
    Dmap = Dmap + squeeze(gk(:,:,2,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
end
D_yxx = sum(Dmap,'all');

figure()
surf(Kx,Ky,Dmap,'EdgeColor','none')

 clim([-400,400]) %  for Ef=0.008+0.008;
 % clim([-1200,1200]) %  for Ef=-0.01+0.008; 

cmax =400;
clim([-cmax,cmax])

cmap = flipud(slanCM('RdBu'));
n = size(cmap,1);
frac = 0.0001;
mid = (n+1)/2; w = round(frac*n);
cmap(max(1,floor(mid-w)):min(n,ceil(mid+w)),:) = 1;
colormap(cmap);
colorbar


 xlim([min(Kx,[],'all'),max(Kx,[],'all')])
 ylim([min(Ky,[],'all'),max(Ky,[],'all')])
 % clim([-4e18,4e18])
 box on;
 view(2)
 xlabel('$K_x$','Interpreter','latex')
 ylabel('$K_y$','Interpreter','latex')
 ax=gca;
 ax.LineWidth=1;
 shading interp
 ax.YAxis.FontSize=18;
 ax.XAxis.FontSize=18;


% 
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gyx_plane.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gyx_plane_white.png', ...
%      'BackgroundColor','none', 'Resolution', 1200);

% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gyx_plane_negef01.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_B1meV_gyx_plane_negef01_white.png', ...
%      'BackgroundColor','none', 'Resolution', 1200);
%%
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_sigma_TDos_all_30K.mat";
load(filename)
Emin=-0.1;
Emax=0.1;
NEF=3000;
Ef_list   = linspace(Emin,Emax,NEF);
%%
z=1000*squeeze(sigma_abc_all(2,1,1,:,:)*pi);
kx=squeeze(TDos_new_all*1.8/1e12);
Blist=-100:100;
ky=kron(ones(3000,1),Blist/100);

figure()
surf(kx,ky,z,'EdgeColor','none')
% pcolor(z')
colorbar
colormap(slanCM('RdBu'))
colormap(flipud(colormap));

shading interp
view(2)
 xlabel('n $(10^{12}$ cm$^{-2})$','Interpreter','latex')
 ylabel('B$_{eff}$ (meV)','Interpreter','latex')
 ax=gca;
 ax.LineWidth=1;
 shading interp
 ax.YAxis.FontSize=18;
 ax.XAxis.FontSize=18;
clim([-0.3,0.3]*pi)
xlim([-8,8])

% set(gcf,'Color','none');      % figure 背景透明 
% set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca, 'tit_QMD/QMD_differentB_sigma_plane.png', ...
%      'BackgroundColor','none', 'Resolution', 600);
%%

% --- 1) 造一个"plane-wave/鼓包"表面 ---
nx = 200; ny = 200;
[x,y] = meshgrid(linspace(-2.2,2.2,nx), linspace(-2.2,2.2,ny));

% 一个平滑鼓包 + 轻微起伏（你可以改成 sin/cos 的plane wave）
% x=cos(x)
% y=cos(y)
z =1.5*exp(-1*(x.^2 + 1.9*y.^2)) ...         % 中央鼓包
    + 0.1*sin(1.2*x).*cos(0.9*y);               % 轻微波纹

% 让颜色随高度变化（淡紫 -> 深紫）
c = z;

figure('Color','white');
ax = axes; hold(ax,'on');

% --- 2) 透明彩色面（边线关掉）---
hs = surf(ax, x, y, z, c, 'EdgeColor','none');
shading interp
colormap(ax, purple_white_colormap(256));  % 自定义淡紫色
colormap(ax,skyblue_white_colormap(256))
% 
% colormap(ax,tealwhite_colormap(256))
hs.FaceAlpha = 0.3;                       % 透明度（0~1）

% --- 3a) 叠加等高线（像你图里那种曲线）---
% nLevels = 10;
% [~,hc] = contour3(ax, x, y, z, nLevels, 'LineColor',[0.15 0.15 0.15], 'LineWidth',1.0);
% hc.ZData = hc.ZData + 1e-3;                % 抬高一点避免z-fighting

% --- 3b) 再叠加"网格线"（可选，更像wireframe）---

% --- 2) 选择画多少条线（少量）---
col = [0.2 0.2 0.2]; 
lw  = 0.9;               % 线宽
alphaLine = 0.35;        % 老版本 plot3 没有 EdgeAlpha，只能用浅色代替
Nx = 12;      % 竖向线（固定x）条数
Ny = 15;      % 横向线（固定y）条数
npts = 600;  % 每条线采样点数（越大越连续）

xmin = min(x(:)); xmax = max(x(:));
ymin = min(y(:)); ymax = max(y(:));

xlines = linspace(xmin, xmax, Nx);
ylines = linspace(ymin, ymax, Ny);

% 竖向线：x 固定，y 扫描
yy = linspace(ymin, ymax, npts);
for xi = xlines
    xx = xi*ones(size(yy));
    zz = interp2(x, y, z, xx, yy, 'linear');   % 'linear' 最稳
    plot3(xx, yy, zz + 1e-6, '-', 'Color', col + (1-alphaLine)*(1-col), 'LineWidth', lw);
end

% 横向线：y 固定，x 扫描
xx = linspace(xmin, xmax, npts);
for yi = ylines
    yy = yi*ones(size(xx));
    zz = interp2(x, y, z, xx, yy, 'linear');
    plot3(xx, yy, zz + 1e-6, '-', 'Color', col + (1-alphaLine)*(1-col), 'LineWidth', lw);
end

view(2); axis tight; box on


% surf(x,y,z,'EdgeColor','none'); 

% step = 11;  % 越大线越少（比如 6~12）
% xm = x(1:step:end, 1:step:end);
% ym = y(1:step:end, 1:step:end);
% zm = z(1:step:end, 1:step:end);
% 
% hm = mesh(ax,xm,ym,zm);
% hm.FaceColor = 'none';
% hm.EdgeColor = [0.15 0.15 0.15];
% hm.LineWidth = 0.6;
% hm.EdgeAlpha = 0.35;   % 线条也可半透明


% hm = mesh(ax, x, y, z);
% hm.FaceColor = 'none';
% hm.EdgeColor = [0.15 0.15 0.15];
% hm.LineWidth = 0.6;
% hm.EdgeAlpha = 0.35;                       % 网格线也透明一点

% --- 4) 视角/光照，让它更像示意图 ---
axis(ax,'off'); axis(ax,'tight'); axis(ax,'vis3d')
view(ax, 25, 25);
camproj(ax,'perspective')
set(gcf,'Renderer','opengl')               % 透明+光照更稳

lighting phong
camlight('headlight')
camlight(90,20)
material dull
daspect([1 1 1])   % z 方向看起来更扁（数字越大越扁）
%
% --- 选波包中心（自己改 x0,y0）---
x0 = 0;
y0 = -0.5;

% 让单极子贴在表面上方一点点
z0 = interp2(x, y, z, x0, y0) -0.5 ;

% --- 单极子颜色/尺寸 ---
% col = [0.65 0.20 0.20];   % 红色
col = [0.10 0.45 0.55];
ms  = 120;                % 中心点大小
L   = 0.15;               % 箭头长度（按你的坐标尺度调）

% 中心"单极子"（红点/小球）
% scatter3(x0, y0, z0, ms, 'filled', ...
%     'MarkerFaceColor', col, 'MarkerEdgeColor','none');
scatter3(x0,y0,z0, 120, 'filled', ...
  'MarkerFaceColor', col, ...
  'MarkerEdgeColor', min(col+0.25,1), ...  % 边缘提亮一点像高光
  'LineWidth', 1.0);

% plot3(x0,y0,z0,'o','MarkerSize',18,'LineWidth',2,'Color',col);

% --- 8 个平面内放射箭头（更像你截图那种）---
theta = linspace(0, 2*pi, 9); theta(end) = [];     % 8 directions
dx = cos(theta(:));
dy = sin(theta(:));
dz = zeros(size(dx));

h = quiver3(x0*ones(size(dx)), y0*ones(size(dy)), z0*ones(size(dz)), ...
            L*dx, L*dy, L*dz, 0, ...
            'Color', col, 'LineWidth', 2.2, 'MaxHeadSize', 0.7);
set(h, 'Clipping','off');  % 防止箭头被坐标轴裁剪

% --- 可选：再加上下两个箭头（像"3D 单极子"）---
Lu = 0.85*L;
hu = quiver3(x0, y0, z0, 0, 0, +Lu, 0, 'Color', col, 'LineWidth',2.2, 'MaxHeadSize',0.7);
hd = quiver3(x0, y0, z0, 0, 0, -Lu, 0, 'Color', col, 'LineWidth',2.2, 'MaxHeadSize',0.7);
set([hu,hd], 'Clipping','off');

% box on
% axis on
set(gcf,'Color','none');      % figure 背景透明 
set(gca,'Color','none');      % axes 背景透明
exportgraphics(gca, 'tit_QMD/Fig1_QM_BC.png', ...
     'BackgroundColor','none', 'Resolution', 1200);

% % --- 5) 加箭头和文字（可选，按你图里的 |u_k>、J^DC）---
% % 箭头位置选在鼓包附近
% x0 = -0.25; y0 = 0.0;
% z0 = interp2(x,y,z,x0,y0);
% 
% quiver3(ax, x0, y0, z0+0.02, -0.25, 0.15, 0.20, 0, ...
%     'LineWidth',2, 'Color',[0.2 0.2 0.2], 'MaxHeadSize',0.7);
% text(ax, x0-0.45, y0+0.25, z0+0.30, '$|u_k\rangle$', ...
%     'Interpreter','latex','FontSize',16, 'Color',[0.2 0.2 0.2]);
% 
% x1 = 0.35; y1 = 0.1;
% z1 = interp2(x,y,z,x1,y1);
% quiver3(ax, x1, y1, z1+0.02, 0.25, 0.10, 0.18, 0, ...
%     'LineWidth',2, 'Color',[0.7 0.7 0.7], 'MaxHeadSize',0.7);
% text(ax, x1+0.20, y1+0.18, z1+0.28, '$|u_{k^\prime}\rangle$', ...
%     'Interpreter','latex','FontSize',16, 'Color',[0.65 0.65 0.65]);
% 
% % J^DC 小箭头
% quiver3(ax, 0, 0.15, max(z(:))+0.08, 0, 0.0, 0.12, 0, ...
%     'LineWidth',3, 'Color',[0.65 0.15 0.15], 'MaxHeadSize',0.9);
% text(ax, 0.05, 0.18, max(z(:))+0.22, '$J^{\mathrm{DC}}$', ...
%     'Interpreter','latex','FontSize',18, 'Color',[0.65 0.15 0.15]);

% 导出（透明背景）
% exportgraphics(gcf,'schematic_surface.svg','BackgroundColor','none','ContentType','vector','Resolution',600);
%%
% --------- 内部函数：淡紫 colormap ----------
function cmap = purple_white_colormap(n)
    % white -> light purple -> deeper purple
    t = linspace(0,1,n)';
    c1 = [1 1 1];
    c2 = [0.92 0.80 0.92];
    c3 = [0.70 0.45 0.75];
    cmap = (t<=0.6).* (c1 + (c2-c1).*(t/0.6)) + ...
           (t>0.6).*  (c2 + (c3-c2).*((t-0.6)/0.4));
end

function cmap = skyblue_white_colormap(n)
% white -> light sky blue -> deeper sky blue
    t  = linspace(0,1,n)';

    c1 = [1.00 1.00 1.00];   % white
    c2 = [0.56 0.84 0.99];   % very light sky blue
    c3 = [0.35 0.65 0.90];   % deeper sky blue

    cmap = (t<=0.6).* (c1 + (c2-c1).*(t/0.6)) + ...
           (t>0.6).*  (c2 + (c3-c2).*((t-0.6)/0.4));
end

function cmap = tealwhite_colormap(n)
% white -> pale teal -> deep teal
    t  = linspace(0,1,n)';

    c1 = [1.00 1.00 1.00];   % white
    c2 = [0.86 0.96 0.96];   % pale teal
    c3 = [0.10 0.55 0.65];   % deep teal

    cmap = (t<=0.6).* (c1 + (c2-c1).*(t/0.6)) + ...
           (t>0.6).*  (c2 + (c3-c2).*((t-0.6)/0.4));
end



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
    % % 
    % % T = [1,0,0,0,0,0,0,0; ...
    % %      0,0,0,0,0,1,0,0; ...
    % %      0,0,1,0,0,0,0,0; ...
    % %      0,0,0,0,0,0,0,1; ...
    % %      0,1,0,0,0,0,0,0; ...
    % %      0,0,0,0,1,0,0,0; ...
    % %      0,0,0,1,0,0,0,0; ...
    % %      0,0,0,0,0,0,1,0];

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
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
    pairsV = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
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
            xinitial_U = [xinitial_U, correlation];
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
            xinitial_V = [xinitial_V, correlation];
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
            for j = 1:size(pairs, 1)
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
        correlation_temp(col_up_i, row_up_j) = -upup_avg; % 上自旋部分
        correlation_temp(col_dn_i, row_dn_j) = -dndn_avg; % 下自旋部分
        correlation_temp(col_dn_i, row_up_j) = -updn_avg; % 上自旋部分
        correlation_temp(col_up_i, row_dn_j) = -dnup_avg; % 下自旋部分

        % 将当前原子对的矩阵存储到输出列表中
        offsite_V(pair_idx, :, :) = correlation_temp;
    end
end

function projector = get_projector(gs, kpoints,Unk, whichbands)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % gs: geometry
    % kpoints: 面内mesh 维度(knum^2,3)
    % whichbands: 选择投影的band 维度(1,nbands)
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    %
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)
    [nbands,~,knum2]=size(Unk);
    projectors_k = zeros(nbands,nbands,knum2);
    % 构造 projectors_k
    for i = 1:knum2
        projectors_k(:, :, i) = diag(which_bands); % 按对角阵填充
    end

    % 应用投影算符
    vs_conj = conj(Unk);
    vs_transpose = permute(Unk, [1, 2, 1]); % 转置第一和第二维
    projectors_k = pagemtimes(pagemtimes(vs_conj, projectors_k), vs_transpose);


    rvectors=gs.hopr;
    expr=exp(1j*(kpoints*(rvectors*gs.a).'));
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

function check_hermitian_elements(obj)
    % H: Hamiltonian 矩阵，大小为 (nbands, nbands, nhopr)
    % g.hopr: hopping 位置矩阵，大小为 (nhopr, 3)

    [nbands, ~, nhopr] = size(obj.ham);
    
    % 遍历所有 (i, j)
    for i = 1:nbands
        for j = 1:nbands
            % 遍历所有 hopping 位置 R
            for k = 1:nhopr
                R = obj.hopr(k, :);  % 当前 hopping 位置
                Hij = obj.ham(i, j, k);  % H_ij(R)
                
                % 找到对应的 -R 位置
                idx_negR = find(all(obj.hopr == -R, 2));
                
                if isempty(idx_negR)
                    warning('No matching -R found for R = (%d, %d, %d)', R(1), R(2), R(3));
                    continue;
                end
                
                Hji_conj = conj(obj.ham(j, i, idx_negR)); % H^*_ji(-R)

                % 检查是否满足厄米性
                if abs(Hij - Hji_conj) > 1e-5
                    fprintf('Non-Hermitian element found at (i=%d, j=%d) for R = (%d, %d, %d)\n', ...
                            i, j, R(1), R(2), R(3));
                    fprintf('H(%d,%d, [%d %d %d]) = %f + %fi\n', i, j, R(1), R(2), R(3), real(Hij), imag(Hij));
                    fprintf('H*(%d,%d, [%d %d %d]) = %f + %fi\n', j, i, -R(1), -R(2), -R(3), real(Hji_conj), -imag(Hji_conj));
                    fprintf('Difference: %f + %fi\n\n', real(Hij - Hji_conj), imag(Hij - Hji_conj));
                    break;
                end
            end
        end
    end
    
    disp('Hamiltonian elements checked.');
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
end


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