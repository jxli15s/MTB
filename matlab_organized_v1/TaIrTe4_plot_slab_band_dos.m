clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/wannier90_hr_p1.dat','data/TaIrTe4/wannier90_hr_p2.dat');

MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
knum=400;
nslab=1;
Occ=120*nslab;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky] = get_Slab2Dkmesh(g,kxline,kyline,knum);

%% Calculate slab bands
[nbands,~,nrpts]=size(g.ham);
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=51;
efermi=7.3742;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%% Load or Calculate slab plane Enk
% load("data/TaIrTe4/dos/TaIrTe4_Enk-400x400.mat")
load("data/TaIrTe4/dos/TaIrTe4_2sEnk-400x400.mat")
%[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);
%% Calculate the Dos 
plottap=1;
nk=100;
Enum=501;
Emin=-0.1;
Emax=0.1;
eps=0.01;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
Dos=sum(Dos)*(Emax-Emin)/Enum;

%% Plot the Dos and TDos
Dos_new=Dos;
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
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
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
xlim([-2,2])
legend('Density')
%ylim([0,200])
print('data/TaIrTe4/dos/TaIrTe4-Density','-dpng','-r300')

%%
myCluster = parcluster('Processes')
delete(myCluster.Jobs)
