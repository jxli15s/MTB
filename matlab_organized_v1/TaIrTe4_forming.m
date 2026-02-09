clear;
clear all;
%p=parpool(8)
%!!!!!!!!! note that this should be used at tb/matlab directory
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
efermi=-0.4805;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-1,1])
%
set(gcf,'Color','none');      % figure 背景透明
set(gca,'Color','none');      % axes 背景透明
exportgraphics(gca,'for_arpes/band_XGX.svg','BackgroundColor','none');
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'Y','\Gamma','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.4805;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-1,1])
%
set(gcf,'Color','none');      % figure 背景透明
set(gca,'Color','none');      % axes 背景透明
exportgraphics(gca,'for_arpes/band_YGY.svg','BackgroundColor','none');
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.4805;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-1,1])
%
set(gcf,'Color','none');      % figure 背景透明
set(gca,'Color','none');      % axes 背景透明
exportgraphics(gca,'for_arpes/band_RYR.svg','BackgroundColor','none');

%% calculate 3d band structure
data=textread("/Volumes/T9/work/dft/TaIrTe4/1s/tairte_vasp/pbe/noivdw/new/wtool/bulkek_plane-matlab-del.dat");
kx=reshape(data(:,4),301,301);
ky=reshape(data(:,5),301,301);
band_v=reshape(data(:,8),301,301);
band_v= 0.5 * (band_v + fliplr(band_v));

band_c=reshape(data(:,9),301,301);
band_c= 0.5 * (band_c + fliplr(band_c));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx, ky: 分别是 nkx × 1 和 nky × 1 的向量
% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序

       
% 绘图
figure;
hold on;

% 1. 三维能带
% s = surf(kx, ky, Ez, 'EdgeColor', 'none');
% colormap turbo
% shading interp
% alpha(0.9);
Ez1 = band_c;
surf(kx, ky, Ez1, 'EdgeColor', 'none');
Ez2 = band_v;
surf(kx, ky, Ez2, 'EdgeColor', 'none');
lighting phong    % 或 gouraud，phong 更平滑
shading interp    % 表面平滑
colormap(slanCM('RdBu'))
colormap(flipud(colormap));


% surf(kx, ky, Ez1, 'EdgeColor', 'none', 'FaceAlpha', 0.6,'FaceColor', [1,0,0]);
% surf(kx, ky, Ez2, 'EdgeColor', 'none', 'FaceAlpha', 0.6,'FaceColor', [0,0,1]);

% 2. 费米面等高线（Ef）
% contour3(kx, ky, Ez, [Ef Ef], 'k', 'LineWidth', 2);  % Fermi contour

Ef_list = -0.2:0.2:0.6;  % 多个等高值
Ef_list1 =[0.0502 0.0502]  % 0.0137  cb  37.5meV vhS electron side
Ef_list2=[-0.05928 -0.05928] % -0.0093 vt  50meV vhS hole side
% 2. 费米面等高线（Ef）
contour3(kx, ky, Ez1, Ef_list1, 'LineColor', 'red', 'LineWidth', 1.5);
contour3(kx, ky, Ez2, Ef_list2, 'LineColor', 'blue' ,'LineWidth', 1.5);

% 3. 视图与标签
view(3);
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
title(['Bands and Fermi Surface ']);
colorbar;
axis tight;
% axis equal;
box on
xlim([-0.3,0.3])
zlim([-0.2,0.2])
clim([-0.3,0.3])

view(-5,15)

% print(gcf,'for_arpes/band_3d.png','-dpng','-r300');   % 透明是否生效取决于版本/查看器
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx, ky: 分别是 nkx × 1 和 nky × 1 的向量
% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序       
% 绘图
figure('Color','none');
hold on;
Ez1 = band_c;
Ez2 = band_v;
ef_val=[0.0502 0.0502]
% 2. 费米面等高线（Ef）
contour(kx, ky, Ez1, ef_val, 'LineColor', 'red', 'LineWidth', 1.5);
contour(kx, ky, Ez2, ef_val, 'LineColor', 'blue' ,'LineWidth', 1.5);
% 3. 视图与标签
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
title(['Fermi Surface for vhs at electron side']);
% colorbar;
axis tight;
axis equal;
box on

set(gcf,'Color','none');      % figure 背景透明
set(gca,'Color','none');      % axes 背景透明
exportgraphics(gca,'for_arpes/fermi_surface_vhs_c.svg','BackgroundColor','none');
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx, ky: 分别是 nkx × 1 和 nky × 1 的向量
% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序       
% 绘图
figure('Color','none');
hold on;
Ez1 = band_c;
Ez2 = band_v;
ef_val=[-0.05928 -0.05928]
% 2. 费米面等高线（Ef）
contour(kx, ky, Ez1, ef_val, 'LineColor', 'red', 'LineWidth', 1.5);
contour(kx, ky, Ez2, ef_val, 'LineColor', 'blue' ,'LineWidth', 1.5);
% 3. 视图与标签
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
title(['Fermi Surface for vhs at hole side']);
% colorbar;
axis tight;
axis equal;
box on

set(gcf,'Color','none');      % figure 背景透明
set(gca,'Color','none');      % axes 背景透明
% exportgraphics(gca,'for_arpes/fermi_surface_vhs_v.svg','BackgroundColor','none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx, ky: 分别是 nkx × 1 和 nky × 1 的向量
% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序       
% 绘图
ef_list = linspace(-1, 0.2, 1201);

Ez1 = band_c;
Ez2 = band_v;

outdir = 'for_arpes';
if ~exist(outdir, 'dir'); mkdir(outdir); end

for idx = 1:length(ef_list)
    ef_val = ef_list(idx);

    % Create an invisible figure (no window popup)
    fig = figure('Visible','off');
    ax  = axes(fig);
    hold(ax,'on');

    % Contour levels: single level is enough
    contour(ax, kx, ky, Ez1, [ef_val ef_val], 'LineColor','red',  'LineWidth',1.5);
    contour(ax, kx, ky, Ez2, [ef_val ef_val], 'LineColor','blue', 'LineWidth',1.5);

    xlabel(ax,'k_x'); ylabel(ax,'k_y');
    title(ax, sprintf('Fermi Surface (E_f = %.3f)', ef_val));
    axis(ax,'tight'); axis(ax,'equal'); box(ax,'on');

    % Transparent backgrounds (works best for svg/pdf; png may depend on version)
    set(fig,'Color','none');
    set(ax,'Color','none');

    % Output filename
    filename = fullfile(outdir, sprintf('fermi_surface_%03d_Ef_%+.3f.svg', idx, ef_val));

    % Export
    exportgraphics(ax, filename, 'BackgroundColor','none');

    % Close figure to avoid memory buildup
    close(fig);
end
