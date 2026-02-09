%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the geomtery information                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;
clear;
%%
g=get_g_from_wannier_encut20("Rhg-15s-wannier");
Electric_field_in_evpA=0.000/46.8968; %0.0544 for zero indirect gap ef=0.0272 
g=add_elec(g,Electric_field_in_evpA);
efermi=get_ef(g);%0.0544eV gap ef=0.0272; 0eV gap ef=0.0238;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[1/3,1/3,0.0]*0.9,...
          [1/3,1/3,0.0],...
          [1/3,1/3,0.0]+([0.5,0.0,0.0]-[1/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;

% efermi=0.0272;
% Electric_field_in_evpA=0.00; %0.08-0.12 V/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
% MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
MTB.plot.plot_bands(Energy*1000,nbands,efermi*1000,kpath,labels,kindex,"bulk",0)
ylim([-200,200])
ylabel('Energy (meV)')
%%
tic;
knum=301;
kxline=[1/3-0.04,1/3+0.04];
kyline=[1/3-0.04,1/3+0.04];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
Enk=Enk-efermi;
toc;
%%
figure()
surf(Enk(:,:,14),'EdgeColor','none')
surf(Enk(:,:,15),'EdgeColor','none')
hold on;
surf(Enk(:,:,16),'EdgeColor','none')
surf(Enk(:,:,17),'EdgeColor','none')
zlim([-0.05,0.05])
%% Calculate the Dos 
plottap=1;
nk=knum;
Enum=801;
Emin=-0.2;
Emax=0.2;
eps=(Emax-Emin)/Enum;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
% Dos=sum(Dos)*(Emax-Emin)/Enum;
Dos=Dos*(Emax-Emin)/Enum;

%% Plot the Dos and TDos
Dos_new=Dos;
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16*(0.1)^2*6;
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos*10^15,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
% plot(TDos_new,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,TDos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(TDos_new,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('States/eV')
% xlim([-2,2])
%yticks([0 5 10])
% ylim([0,100])
%yticklabels({})
legend('DOS')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
print('data/TaIrTe4/dos/TaIrTe4-DOS','-dpng','-r300')

%%
figure('Color','white')
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16;
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
% xlim([-2,2])
legend('Density')
%ylim([0,200])
print('data/TaIrTe4/dos/TaIrTe4-Density','-dpng','-r300')


%%
load("data/Graphene/15s/e_n_dos/rhg_dos_e_n_large_scal.mat")
Eaxis=linspace(-0.05,0.05,501);
density=TDos_all/norm(cross(g.a(1,:),g.a(2,:)))*10^16*(0.08)^2*6;
E_y=kron(ones(1,size(density,2)),E_list');
E_x=kron(ones(size(density,1),1),Eaxis);
%%
% load("data/Graphene/15s/e_n_dos/rhg_dos_e_n.mat")
% load("data/Graphene/15s/fplo/encut_25/dos/rhg_dos_e_n_fplo_encut25.mat")
% E_list=linspace(0,3e-3,601)*46.8968;

load("data/Graphene/15s/fplo/encut_25/dos/rhg_dos_e_n_fplo_encut25_801.mat")
E_list=linspace(0,3e-3,1201)*46.8968;
Eaxis=linspace(-0.12,0.12,1001);
density=TDos_all/norm(cross(g.a(1,:),g.a(2,:)))*10^16*(0.08)^2;
E_y=kron(ones(1,size(density,2)),E_list');
E_x=kron(ones(size(density,1),1),Eaxis);

% desity=flipud(density)
E_y=flipud(E_y);
% Dos_all=flipud(Dos_all)
figure()
surf(density/10^12,E_y*1000,Dos_all,'EdgeColor','none')
set(gca, 'YDir', 'normal'); % 确保 y 轴从小到大显示
% surf(E_x,E_y,Dos_all,'EdgeColor','none')
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp
colorbar;
% clim([0,0.005])
axis tight;
% axis equal;
box on
view(0,90)
% xlim([-1*10^13,1*1*10^13])
xlim([-6,3])
xlabel('Carrier density (10^{12}/cm^2)')
% ylabel('E field (meV)')
ylabel('E_{gap}^K (meV)')
yticks([0,1,2,3,4,5,6,7]*20)
yticklabels({'140','120','100','80','60','40','20','0'})
% pcolor(Dos_all)
%%
knum=801;
kxline=[1/3-0.04,1/3+0.04];
kyline=[1/3-0.04,1/3+0.04];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
%%
figure()
E_y=flipud(E_y);
% surf(density/10^12,E_y*1000,Dos_all,'EdgeColor','none')
surf(E_x,E_y,Dos_all,'EdgeColor','none')
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp
colorbar;
% clim([0,0.005])
axis tight;
% axis equal;
box on
view(0,90)
% xlim([-1*10^13,1*1*10^13])
xlim([-7,2])
xlabel('Carrier density (10^{12}/cm^2)')
ylabel('E field (meV)')
% pcolor(Dos_all)

%%
%% for 0.0meV Efield
figure()
plot(Eaxis,Dos_all(1,:))
figure()
plot(Dos_all(1,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by GreenFunction     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=0
load("data/Graphene/15s/fplo/encut_25/dos/Egap_0/Ez1_mesh_801_Elec_1.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_0/Ez2_mesh_801_Elec_1.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_0/Ez3_mesh_801_Elec_1.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_0/Ez4_mesh_801_Elec_1.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);

% mu=-0.02;
mu=0.002;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%%
mu=-0.01;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.405,1.55])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_x/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=-0.03;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.405,1.55])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_x/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=-0.0025;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.54])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_x/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=0.025;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.54])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_x/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%% for 9.26meV Efield
figure()
plot(Eaxis,Dos_all(80,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_9/Ez1_mesh_801_Elec_80.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_9/Ez2_mesh_801_Elec_80.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_9/Ez3_mesh_801_Elec_80.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_9/Ez4_mesh_801_Elec_80.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=0.00168;
mu=0.0014;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');


%% for 12.427meV Efield
figure()
plot(Eaxis,Dos_all(107,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_12/Ez1_mesh_801_Elec_107.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_12/Ez2_mesh_801_Elec_107.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_12/Ez3_mesh_801_Elec_107.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_12/Ez4_mesh_801_Elec_107.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
mu=0.0007;

eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%% for 15.71meV Efield
figure()
plot(Eaxis,Dos_all(135,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_15/Ez1_mesh_801_Elec_135.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_15/Ez2_mesh_801_Elec_135.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_15/Ez3_mesh_801_Elec_135.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_15/Ez4_mesh_801_Elec_135.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=0.0005;
mu=0
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%% for 18.99meV Efield
figure()
plot(Eaxis,Dos_all(163,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_18/Ez1_mesh_801_Elec_163.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_18/Ez2_mesh_801_Elec_163.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_18/Ez3_mesh_801_Elec_163.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_18/Ez4_mesh_801_Elec_163.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=0;
mu=-0.00048
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%% for 22.39meV Efield
figure()
plot(Eaxis,Dos_all(192,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_22/Ez1_mesh_801_Elec_192.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_22/Ez2_mesh_801_Elec_192.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_22/Ez3_mesh_801_Elec_192.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_22/Ez4_mesh_801_Elec_192.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=-0.00072;
mu=-0.001;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');


%% for 25.67meV Efield
figure()
plot(Eaxis,Dos_all(220,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Plot the fermi surface by GreenFunction         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=9.26meV E_list(80)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_25/Ez1_mesh_801_Elec_220.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_25/Ez2_mesh_801_Elec_220.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_25/Ez3_mesh_801_Elec_220.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_25/Ez4_mesh_801_Elec_220.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=-0.00144;
mu=-0.0019;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%% for 30meV Efield
figure()
plot(Eaxis,Dos_all(257,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by GreenFunction     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=30meV E_list(257)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_30/Ez1_mesh_801_Elec_257.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_30/Ez2_mesh_801_Elec_257.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_30/Ez3_mesh_801_Elec_257.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_30/Ez4_mesh_801_Elec_257.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
% mu=-0.006;
mu=-0.00312;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=0.02;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%%
mu=0.0;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=0.03;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%% for 80meV Efield
figure()
plot(Eaxis,Dos_all(683,:))
% xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by GreenFunction     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=80meV E_list(683)
load("data/Graphene/15s/fplo/encut_25/dos/Egap_80/Ez1_mesh_801_Elec_683.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_80/Ez2_mesh_801_Elec_683.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_80/Ez3_mesh_801_Elec_683.mat");
load("data/Graphene/15s/fplo/encut_25/dos/Egap_80/Ez4_mesh_801_Elec_683.mat");
Enk=cat(3,Ez1,Ez2,Ez3,Ez4);
%%
mu=-0.01128;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%%
mu=-0.005;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=0.03;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
mu=0.042;
eta=2*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);

% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
xlim([1.41,1.545])
ylim([0.78,0.93])
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

%%
figure()
plot(density(1,:),Dos_all(1,:))
xlabel('Carrier density (n)')
ylabel('DOS(1/cell)')
%
efermi=Ef_array(1);
g=get_g_from_wannier_encut20("Rhg-15s-wannier");
Electric_field_in_evpA=0.00;
g=add_elec(g,Electric_field_in_evpA);

%% for 0.2meV Efield
figure()
plot(Eaxis,Dos_all(81,:))
xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
figure()
plot(density(81,:),Dos_all(81,:))
xlabel('Carrier density (n)')
ylabel('DOS(1/cell)')
%% for 0.3meV Efield
figure()
plot(Eaxis,Dos_all(121,:))
xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
figure()
plot(density(121,:),Dos_all(121,:))
xlabel('Carrier density (n)')
ylabel('DOS(1/cell)')
%% for 0.5meV Efield
figure()
plot(Eaxis,Dos_all(201,:))
xlim([-0.12,0.12])
xlabel('Chemical potential(eV)')
ylabel('DOS(1/cell)')
figure()
plot(density(201,:),Dos_all(201,:))
xlabel('Carrier density (n)')
ylabel('DOS(1/cell)')
%% for 1meV Efield
figure()
plot(Eaxis,Dos_all(401,:))
figure()
plot(density(401,:),Dos_all(401,:))
%% for 2meV Efield
figure()
plot(Eaxis,Dos_all(801,:))
figure()
plot(density(801,:),Dos_all(801,:))

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi Surface-1             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate plane bands
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Construct the geomtery information                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%set the parameters
g=get_g_from_wannier_encut20("Rhg-15s-wannier");
Electric_field_in_evpA=0.000;
g=add_elec(g,Electric_field_in_evpA);
%efermi=get_ef(g);
%%
knum=301;
kxline=[1/3-0.04,1/3+0.04];
kyline=[1/3-0.04,1/3+0.04];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
Enk=Enk-efermi;
%%
% save('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_0meV.mat',"Enk",'efermi','Kx','Ky')
% load('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_0meV.mat')
% save('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_02meV.mat',"Enk",'efermi','Kx','Ky')
% load('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_02meV.mat')
% save('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_03meV.mat',"Enk",'efermi','Kx','Ky')
% load('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_03meV.mat')
% save('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_05meV.mat',"Enk",'efermi','Kx','Ky')
% load('data/Graphene/15s/e_n_dos/Enk_mesh_801_e_05meV.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kx=Kx;
ky=Ky; 
%%
% 绘图
%for E=0
mu=0.002 % peak1
mu=-0.05;
%for E=0.2meV
% mu=-0.001638;  %for E=0 peak1
% mu=-0.003;  %for E=0 peak1 left
% mu=-0.008;  %for E=0 peak1 left
% mu=-0.016;  %for E=0 peak1 left
% mu=-0.10705;
% mu=0.002;
%for E=0.3meV
% mu=-0.0008;  %for E=0 peak1
% mu=-0.003237;  %for E=0 peak2
% mu=-0.015;
% mu=-0.03;
% mu=-0.10690;
% mu=-0.10705;
%for E=0.5meV
% mu=-0.004805;  %for E=0 peak1
% mu=-0.017;  %for E=0 peak1
% mu=-0.10496;  %for E=0 peak1
% mu=0.01

figure('Color', 'w');
hold on;
Ez1 = squeeze(Enk(:,:,15))';
Ez2 = squeeze(Enk(:,:,16))';
Ez3 = squeeze(Enk(:,:,14))';
Ez4 = squeeze(Enk(:,:,17))';
%%
figure()
Ef_list=[0 0]; % for E=0 peak1

% 1. 费米面等高线（Ef）
contour3(kx, ky, Ez1, Ef_list, 'LineColor', 'blue', 'LineWidth', 1.5);
contour3(kx, ky, Ez2, Ef_list, 'LineColor', 'red' ,'LineWidth', 1.5);
contour3(kx, ky, Ez3, Ef_list, 'LineColor', 'cyan' ,'LineWidth', 1.5);
contour3(kx, ky, Ez4, Ef_list, 'LineColor', 'black' ,'LineWidth', 1.5);
% 2. 视图与标签
view(3);
% xlabel('k_x/A'); ylabel('k_y/A'); zlabel('E(k)');
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');

title('Fermi Surface');
box on
axis equal;
% xlim([1.4,1.55])
% ylim([-0.91,-0.78])
% 
% xlim([1.38,1.57])
% ylim([-0.94,-0.76])
view(0,90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by GreenFunction     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for E=0
mu=0.002;
mu=-0.05
% mu=0.00048;  %for E=0 peak1
% mu=-0.0091;  %for E=0 peak1 left
% mu=-0.02;  %for E=0 peak1 left2
% mu=-0.092;  %for E=0 peak1 left3
% mu=-0.10768;  %for E=0 peak1 left3
%for E=0.2meV
% mu=-0.001638;  %for E=0.2meV peak1
% mu=-0.003;  %for E=0.2meV peak1 left
% mu=-0.008;  %for E=0 peak1 left
% mu=-0.016;  %for E=0 peak1 left
% mu=-0.10705;
%for E=0.3meV
% mu=-0.0008;  %for E=0 peak1
% mu=-0.003237;  %for E=0 peak2
% mu=-0.015;
% mu=-0.03;
%for E=0.5meV
% mu=-0.004805;  %for E=0 peak1
% mu=-0.017;  %for E=0 peak1
% mu=-0.06;  %for E=0 peak1
% mu=-0.10496;  %for E=0 
% mu=0.01;  %for E=0 peak1
 % mu=0.0175;  %for E=0 peak1
 % mu=-0.001;
eta=10*10^-4;
A = compute_spectral_function_kplane(Enk, mu, eta);
%% plot by pcolor
 figure()
 imagesc(A'); axis equal; axis off; colormap hot; colorbar;
%% plot by imagesc
figure('Color', 'w');
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
% colormap hot;        % 色图
shading interp
colorbar;
% axis equal tight;
axis equal;
% xlim([1.4,1.55])
% ylim([-0.91,-0.78])

% xlim([1.38,1.57])
% ylim([-0.94,-0.76])

title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu),'eV']);
xlabel('$k_y/\AA$', 'Interpreter', 'latex');
ylabel('$k_y/\AA$', 'Interpreter', 'latex');
%% plot by surf
figure('Color', 'w');
surf(Kx, Ky, A', 'EdgeColor', 'none');
% colormap turbo;
colorbar;
view([0, 90]);       % 视角角度
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('$k_y$', 'Interpreter', 'latex');
zlabel('$A(k, \mu)$', 'Interpreter', 'latex');
title(['Spectral Function A(k, \mu = ', num2str(mu), ')']);
axis equal;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[1/3,1/3,0.0]*0.9,...
          [1/3,1/3,0.0],...
          [1/3,1/3,0.0]+([0.5,0.0,0.0]-[1/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;

% efermi=0
% Electric_field_in_evpA=0.00; %0.08-0.12 V/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
ylim([-0.2,0.2])
%
%%n=3.1*10^12 n=6.2*10^12
%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% efermi=min(Energy(16,:));
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)
g.iniham=g.ham+0;
ylim([-0.1,0.1])
% bands_file = 'data/Graphene/dft-bands/dft-bands-3s.dat';
bands_file = 'data/Graphene/dft-bands/dft-bands-15s-small.dat';
bands_data = readmatrix(bands_file);  % [nq × nmode]
bands_data(:,2)=bands_data(:,2)-3.0030+0.016;
kpath_vasp=bands_data(:,1);
bands_vasp=bands_data(:,2:end);
kpath_vasp=reshape(kpath_vasp,[],192);
bands_vasp=reshape(bands_vasp,[],192);
bands_vasp=bands_vasp(:,46:75);
% bands_vasp=bands_vasp(:,60:61);
kpath_vasp(end/2,:)=[];
bands_vasp(end/2,:)=[];
hold on;
% plot(bands_data(:,1),bands_data(:,2),'r.',LineWidth=2)
for i=1:size(bands_vasp,2)
    plot(kpath_vasp(:,i),bands_vasp(:,i),'r-',LineWidth=2)
end
ylim([-1,1])

%%
function g=get_g_from_wannier_encut20(name)
%Create geometry object
g = MTB.geometry(name);
% Read geometry from POSCAR file
g = MTB.read_poscar(g, fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "POSCAR"));
% Read Wannier Hamiltonian data
[g.ham, g.hopr] = MTB.wannier.read_hr(...
    fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "wannier90_hr_p1.dat"), ...
    fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "wannier90_hr_p2.dat"));
g.wpos=g.atoms*g.a;
end

function ef=get_ef(g)
    knum=401;%501
    % kxline=[-0.5,0.5];
    % kyline=[-0.5,0.5];
    kxline=[-0.1+1/3,0.1+1/3];
    kyline=[-0.1+1/3,0.1+1/3];
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

function A = compute_spectral_function_kplane(Enk, mu, eta)
% 输入参数
% Enk: [nkx, nky, nbands] 的能带数组
% mu: 能量点
% eta: Lorentzian 展宽参数，例如 eta = 0.01

    [nkx, nky, nbands] = size(Enk);
    A = zeros(nkx, nky);  % 初始化谱函数

    for ix = 1:nkx
        for iy = 1:nky
            for n = 1:nbands
                E = Enk(ix, iy, n);
                A(ix, iy) = A(ix, iy) + eta / ( (mu - E)^2 + eta^2 );
            end
        end
    end

    A = A / pi;  % 加上 prefactor 1/pi
end

