clear;
clear all;
g = MTB.geometry("Haldane");
g = MTB.read_poscar(g,"data/Haldane/POSCAR");
[g.ham,g.hopr] = MTB.read_hr('data/Haldane/Haldane_hr.dat');
 e=1.6*10^-19;
h=6.626*10^-34;
G=25812;%h/e^2(Om*Cm)^-1
coef=1/G*10^8;
%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'K','\Gamma','M'}; % labels for k
hkpoints={[0.333333,0.333333,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=0.0;
%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)
%% Calculate plane bands
knum=201;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);

%% Calculate Dos
nk=knum;
Enum=501;
Emin=-10;
Emax=10;
eps=(Emax-Emin)/Enum;
plottap=0;
efermi=0.0;

[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);
plot(Eaxis,Dos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
plot(Eaxis,TDos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
%plot(Eaxis,TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16,'Linestyle','-','Color','#4DA1D7','LineWidth',2)


%% Calculate Berry Curvature by LOOP method
plottap=1;
bandindex=1;
[Omega_k,KX,KY] = MTB.ham.get_Berry_curvature(bandindex,Unk,Kx,Ky,plottap);

%% Calculate Berry Curvature by Quantum metric
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,1:1,0.01,1);
%% Calculate Quantum Geometry at a single point
% k 是当前 k 点（1x2），Uk, Ek 已经从 H(k) 对角化得到
dk_list = [g.b(1,:)/knum;  % 方向 1
           g.b(2,:)/knum]; % 方向 2

delta = 1e-6;      % 规避近简并
band  = 5;         % 目标能带号举例
k=[0.1,0.1,0.0];
Uk=squeeze(Unk(:,:,1));
Ek=squeeze(Enk(1,1,:));
band=1;
delta=1e-6;
[g, Q, F] = MTB.ham.quantum_geometry_general_k(g, k, dk_list, Uk, Ek, band, delta);

gxx = g(1,1);
gyy = g(2,2);
gxy = g(1,2);      % = g(2,1)
Omega_z = F(1,2);  % F_12 = -F_21，即 2D 的 Berry 曲率
%%
%% Calculate Quantum Geometry on a mesh
% k 是当前 k 点（1x2），Uk, Ek 已经从 H(k) 对角化得到
% ===== 已有的部分 =====
knum  = 401;
kxline = [0,1];
kyline = [0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]  = MTB.ham.get_bulk_plane_bands(g, Kx,Ky,Kz);
[Hamk,Unk,Enk]  = MTB.ham.get_bulk_plane_bands_with_Ham(g, Kx,Ky,Kz);
%
% ===== 定义 dk_list (2D 情况) =====
b1 = g.b(1,:);    % 你的 reciprocal vectors，如果是别的字段，就替换
b2 = g.b(2,:);
dk_list = [b1/knum;   % 对应 "kx 方向"的单位向量
           b2/knum];  % 对应 "ky 方向"的单位向量



% ===== 选择要算的能带 & delta =====
band_list = 1:2;          % 比如 valence top + conduction bottom…
delta     = 1e-6;             % 能隙平方阈值，根据能量单位稍微调一下

% ===== 在整张 mesh 上算 g, Q, F =====
[gk, Qk, Fk,~] = MTB.ham.quantum_geometry_general_plane(g, Kx,Ky,Kz, Unk,Enk, band_list, dk_list, delta);

% 例如：取 band_list 中第 1 条带的 g_xx, g_yy, Omega_z
iband = 1;
gxx = squeeze(gk(:,:,1,1,iband));
gyy = squeeze(gk(:,:,2,2,iband));
gxy = squeeze(gk(:,:,1,2,iband));
Omega_z = squeeze(Fk(:,:,1,2,iband));  % F_{xy}
%%
% 1. 先算量子几何
[gk, Qk, Fk, vk] = MTB.ham.quantum_geometry_general_plane( ...
    g, Kx,Ky,Kz, Unk,Enk, band_list, dk_list, delta);

% 2. BZ 面积 & EF 列表
% AreaBZ  = norm(cross(g.b(1,:), g.b(2,:)));   % 2D BZ 面积
AreaBZ = norm(g.b(1,:))*norm(g.b(2,:));
Emin=-2;
Emax=2;
NEF=30;
Ef_list = linspace(Emin, Emax, NEF);
eta     = 1e-2;
weights = [];    % 均匀网格就传空

% 3. 得到 k-resolved 和 band-resolved 的 D_abc^(n)
[Dk_abc_n, D_abc_n] = MTB.ham.get_D_metric_kresolved_from_gv( ...
    gk, vk, Enk, band_list, Ef_list, AreaBZ, eta, weights);

% 例如：第 iband 条带的 D_{y;xx}^{(n)}(E_F)
a = 2; b = 1; c = 1;      % 2->y, 1->x (按你的 dk_list 约定)
iband = 1;
D_yxx_band1_EF = squeeze(D_abc_n(a,b,c,iband,:));  % 长度 N_EF

% 某个 E_F 下的 k-resolved D_{y;xx}^{(n)}(k;E_F)
ief = 3;
D_yxx_k = squeeze(Dk_abc_n(:,:,a,b,c,iband,ief));
figure; pcolor(Kx,Ky,D_yxx_k); shading interp; axis equal tight; colorbar;
shading interp;          % 平滑一下，不要格子线
colorbar;                % 右边加个颜色条
colormap("hot");           % 颜色图随意换：parula, jet, etc.
title(sprintf('D^{(n)}_{y;xx}(k; E_F=%.3f eV)', Ef_list(ief)));
%%
figure()
pcolor(Kx,Ky,squeeze(Dk_metric(:,:,1,2,2,1,3)))

shading interp;          % 平滑一下，不要格子线
colorbar;                % 右边加个颜色条
colormap("hot");           % 颜色图随意换：parula, jet, etc.
%% plot quantum metric
figure();
% 如果 Omega_z 是按 (ix,iy) 存的，大小和 Kx,Ky 一样：
pcolor(Kx, Ky, gxx);

shading interp;          % 平滑一下，不要格子线
colorbar;                % 右边加个颜色条
colormap("hot");           % 颜色图随意换：parula, jet, etc.

axis equal tight;        % 等比例 + 紧凑
xlabel('k_x');
ylabel('k_y');
title('\Omega_z(k_x,k_y)');  % LaTeX 风格标题
b1 = g.b(1,:);   % 1x3
b2 = g.b(2,:);   % 1x3
%
% Brillouin zone area
AreaBZ = norm(g.b(1,:))*norm(g.b(2,:));
% AreaBZ = norm(cross(b1, b2));    % |b1 x b2|

[Nkx, Nky] = size(Kx);
DeltaS = AreaBZ / (Nkx * Nky);   % 每个 k 点对应的面积

% 选一条能带的 Berry curvature 网格
g_all = squeeze(gxx+gyy);  % F_{xy}(k_x,k_y) for band iband

% 数值积分得到 Chern number
g_trace = (DeltaS) * sum(g_all, 'all');

fprintf('Band %d: gxx+gyy ≈ %.6f\n', band_list(iband), g_trace);

%% plot Berry curvature
figure();
% 如果 Omega_z 是按 (ix,iy) 存的，大小和 Kx,Ky 一样：
pcolor(Kx, Ky, Omega_z);

shading interp;          % 平滑一下，不要格子线
colorbar;                % 右边加个颜色条
colormap("hot");           % 颜色图随意换：parula, jet, etc.

axis equal tight;        % 等比例 + 紧凑
xlabel('k_x');
ylabel('k_y');
title('\Omega_z(k_x,k_y)');  % LaTeX 风格标题
%%
% reciprocal primitive vectors
b1 = g.b(1,:);   % 1x3
b2 = g.b(2,:);   % 1x3

% Brillouin zone area
AreaBZ = norm(g.b(1,:))*norm(g.b(2,:));
% AreaBZ = 
% AreaBZ = norm(cross(b1, b2));    % |b1 x b2|

[Nkx, Nky] = size(Kx);
DeltaS = AreaBZ / (Nkx * Nky);   % 每个 k 点对应的面积
DeltaS = norm(dk_list(1,:))*norm(dk_list(2,:));
% 选一条能带的 Berry curvature 网格
Omega_z = squeeze(Fk(:,:,1,2,iband));  % F_{xy}(k_x,k_y) for band iband

% 数值积分得到 Chern number
Chern = (DeltaS / (2*pi)) * sum(Omega_z, 'all');

% 由于数值误差，一般会离整数有一点点偏离，取最近的整数
Chern_rounded = round(real(Chern));

fprintf('Band %d: Chern ≈ %.6f, rounded = %d\n', band_list(iband), Chern, Chern_rounded);

%% calcualte the nonlinear Hall sigma_abc by quantum metric dipole 
% 1. 准备好 Hk, Unk, Enk
%    Hk: nb x nb x Nkx x Nky  (你可以在构建 bands 的时候一并存下来)

g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/Graphene/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');
%%
% ===== 已有的部分 =====
knum  = 301;
kxline = [1/3-0.05,1/3+0.05];
kyline = [1/3-0.05,1/3+0.05];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]  = MTB.ham.get_bulk_plane_bands(g, Kx,Ky,Kz);
[Hamk,Unk,Enk]  = MTB.ham.get_bulk_plane_bands_with_Ham(g, Kx,Ky,Kz);
%


% 2. dk_vecs: 对应两个方向的 k 步长, 比如沿 b1, b2
% ===== 定义 dk_list (2D 情况) =====
b1 = g.b(1,:);    % 你的 reciprocal vectors，如果是别的字段，就替换
b2 = g.b(2,:);
dk_vecs = [b1/knum;   % 对应 "kx 方向"的单位向量
           b2/knum];  % 对应 "ky 方向"的单位向量


% 3. 其它参数
band_list = 1:2;               % 比如 valence top + conduction bottom…
Emin=-1.35;
Emax=-1.15;
NEF=500;
Ef_list   = linspace(Emin,Emax,NEF);
% AreaBZ    = norm(cross(b1,b2));
AreaBZ = norm(g.b(1,:))*norm(g.b(2,:))*0.1^2;
eta       = 1e-3;
weights   = [];               % 让函数内部设置均匀 AreaBZ/(Nkx*Nky)
deltaE_reg = 1e-5;

sigma_abc = MTB.ham.get_sigma_quantum_metric_dipole( ...
    Hamk, Unk, Enk, band_list, dk_vecs, Ef_list, ...
    AreaBZ, eta, weights, deltaE_reg);
%%
% 例如 σ_{yxx}(E_F)
a = 2; b = 1; c = 1;
sigma_yxx = squeeze(sigma_abc(a,b,c,:));
%%
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            plot(Ef_list,squeeze(sigma_abc(x,y,z,:)))
        end
    end
end
% plot(sigma_yxx)

%% Calculate Hall Conductivity
Enum=101;
Emin=-10;
Emax=10;
T=0.001
Eaxis=linspace(Emin,Emax,Enum);
sigma=zeros(1,Enum);
for i=1:Enum
    E=Eaxis(i);
    %sigma(i)=sum(Omega_k(Enk<E));
    [Eaxis,sigma]=MTB.ham.get_ahc(Omega_k,Enk,Enum,Emin,Emax,T);
end
%%
figure('Color','white')
hold on;
plot(Eaxis,sigma,'LineWidth',2)
ylabel('$\sigma$ ($\Omega cm)^{-1}$','Interpreter','latex','FontSize',24)
xlabel('Electric field(V/nm)','FontSize',24)
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
%ylim([0.0,70])
box on

%% Calculate anc 
Enum=101;
Emin=-2;
Emax=2;
T=60
Eaxis=linspace(Emin,Emax,Enum);
sigma=zeros(1,Enum);
for i=1:Enum
    E=Eaxis(i);
    %sigma(i)=sum(Omega_k(Enk<E));
    [Eaxis,sigma]=MTB.ham.get_anc(Omega_k,Enk,Enum,Emin,Emax,T);

end

figure('Color','white')
hold on;
plot(Eaxis,sigma,'LineWidth',2)
ylabel('$\sigma$ ($\Omega cm)^{-1}$','Interpreter','latex','FontSize',24)
xlabel('Electric field(V/nm)','FontSize',24)
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
%ylim([0.0,70])
box on


%% Calculate bulk bands
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
efermi=0.0258; %% set Fermi Level 0.0258
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"MnBiTe-bulk",0)
%% Calculate bulk bands add electric field
g.readwpos("data/MnBiTe-xu/wpos_MnBiTe");
g.wpos=g.wpos*27
Electric_field_in_evpA=0.0024; %eV/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b)
E_gate=0.0024
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig10",E_gate*10000)
%%
E_gate=0.0024
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig",E_gate*10000)
%% Calculate bulk bands add electric fields
g.readwpos("data/MnBiTe-xu/wpos_MnBiTe")
g.wpos=g.wpos*27
Electric_field_list=linspace(-0.003,0.003,11);
nk=101
for i=1:11
    Electric_field_in_evpA=Electric_field_list(i)
    [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b)
    MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig"+num2str(i),Electric_field_in_evpA*10000)
%    Egap(i)=(min(Energy(172,:))-max(Energy(171,:)))*1000
    Egap(i)=(min(Energy(175,:))-max(Energy(174,:)))*1000
    if Egap(i)<1
        Egap(i)=0
    end
end
%%
figure('Color','white')
hold on;
plot(Electric_field_list*10,Egap,'--o','LineWidth',2)
ylabel('E_{gap} (meV)','FontSize',24)
xlabel('Electric field(V/nm)','FontSize',24)
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
ylim([0.0,70])
box on
print("EvsGap",'-dpng','-r600')

%% DOS cals
nk=101;
Electric_field_in_evpA=0;
Energy=MTB.ham.get_bulk_plane_bands_add_electric(g.ham, g.hopr, g.wpos, Electric_field_in_evpA, nbands,nrpts,nk,g.a,g.b);

%%
nk=101;
Enum=101;
Emin=-0.2;
Emax=0.2;
eps=0.0025;
Nband=276;
plottap=0;
Electric_field_list=linspace(-0.003,0.003,11);
Dos_all=zeros(length(Electric_field_list),Enum);
TDos_all=zeros(length(Electric_field_list),Enum);
for i=1:11
    Electric_field_in_evpA=Electric_field_list(i);
    Energy=MTB.ham.get_bulk_plane_bands_add_electric(g.ham, g.hopr, g.wpos, Electric_field_in_evpA, nbands,nrpts,nk,g.a,g.b);
    [Eaxis,Dos,TDos]=MTB.ham.get_dos(Energy-efermi,eps,Enum,Emin,Emax,Nband,nk,plottap);
    Dos_all(i,:)=Dos;
    TDos_all(i,:)=TDos;
end
%%
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_all(6,:),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('Dos')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%% Plot Dos
figure('Color','white')
hold on;
TDos_all_new=TDos_all/16.5012*10^16
for i=1:1:11
    plot(Eaxis,TDos_all_new(i,:),'Linestyle','-','LineWidth',2,'Color',[(i-1)/10,100/255 255/255])
end
labels=Electric_field_list
%Eaxis_all=repmat(Eaxis,11,1)
%plot(Eaxis_all',TDos_all','LineStyle','-')
box on;
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
ylabel('n (cm^{-2})','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
legend('-3.0meV','-2.4meV','-1.8meV','-1.2meV','-0.6meV','0meV','0.6meV','1.2meV','1.8meV','2.4meV','3.0meV','Location','best','NumColumns',2)
legend('boxon')
legend('FontSize',13)
xlim([-0.1,0.1])
ylim([0.0,3*10^13])
print('MnBiTe-Density','-dpng','-r300')

%%

%%
figure('Color','white')
hold on;
for i=1:11
    plot(Eaxis,Dos_all(i,:)+i*0.00,'Linestyle','-','LineWidth',2,'Color',[i/11,30/255 255/255])
end
xline(0,'LineStyle','--','LineWidth',2)
legend('-3.0meV','-2.4meV','-1.8meV','-1.2meV','-0.6meV','0meV','0.6meV','1.2meV','1.8meV','2.4meV','3.0meV','Location','northeast','NumColumns',2)
% legend('boxoff')
legend('FontSize',13)
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xlim([-0.1,0.1])
%ylim([0.0,0.08])
yticklabels('')
ylabel('DOS','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
box on
print('MnBiTe-DOS','-dpng','-r300')