clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/wannier90_hr_p1.dat','data/TaIrTe4_2d/wannier90_hr_p2.dat');
% g = MTB.read_poscar(g,"data/TaIrTe4/qe/gamma_low/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
% Volumes/T9/work/tb/matlab/data/TaIrTe4/qe/gamma_low
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.4423;
nk=101;

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

Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

%% Berry Curvature
nk=31;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0;
% nsband=1:88;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
%%
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,1);

%%
figure('Color','White')
contourf(Kx,Ky,Omega_k(:,:,88))
%% Berry Curvature Dipole
nk=31;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0.1;
nsband=1:100;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,0);
%%
Enum=101;
Emin=-1;
Emax=1;
T=0.01;
[Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband)-efermi,Enum,Emin,Emax,T);
figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz} (\AA)$','Interpreter','latex','FontSize',24)
%yrange=[min(sigma)-20,max(sigma)+20]
%xrange=[-1,1]
%xlim(xrange)
%ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%
[Omega_dk2,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,1);
%%
figure('Color','White')
contourf(Kx,Ky,Omega_dk2(:,:,88))
%%%BCD along Eaxis

[Eaxis,bcd]=get_bcd(Omega_dk,Enk,Enum,Emin,Emax,tem)
%%

MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
knum=400;
nslab=3;
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
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%% Load or Calculate slab plane Enk
% load("data/TaIrTe4/dos/TaIrTe4_Enk-400x400.mat")
%[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);
%% Calculate the Dos 




plottap=1;
nk=knum;
Enum=401;
Emin=-2;
Emax=2;
eps=0.05;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
sum(Dos)*(Emax-Emin)/Enum

%% Plot the Dos and TDos
Dos_new=Dos
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('States/eV')
xlim([-2,2])
%yticks([0 5 10])
ylim([0,100])
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
