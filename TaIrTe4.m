clear;
clear all;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.read_hr('data/TaIrTe4/wannier90_hr.dat');
%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=7.3742;
%% plot bulk bands
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-bulk",0)

%% plot slab bands
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
%%
nk=51;
nslab=3;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)


%% calculate slab plane bands

knum=81;
nslab=3;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky] = get_2Dkmesh(g,kxline,kyline,knum);
%% 
[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);
%% Calculate the Dos 
plottap=1
nk=knum;
Enum=101;
Emin=-10;
Emax=10;
eps=(Emax-Emin)/Enum;
efermi=7.3742;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,Nband,nk,plottap);
%% Plot the Dos and TDos
Dos_new=Dos
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('DOS')
xlim([-2,2])
%yticks([0 5 10])
yticklabels({})
legend('DOS')
print('TaIrTe4-DOS','-dpng','-r300')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%
figure('Color','white')
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
xlim([-2,2])
legend('Density')
%ylim([0,200])
print('TaIrTe4-Density','-dpng','-r300')
%%
myCluster = parcluster('Processes')
delete(myCluster.Jobs)
