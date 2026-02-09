clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

gs = MTB.geometry("TaIrTe4_s");
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y','R'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;

g.wpos=[];
g.wpos=g.atoms*g.a;
orbital_num=[4,4];
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1])
    ]
g.wpos=g.wpos

Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

%% Berry Curvature
nk=31;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0.0;
nsband=1:8;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
%%
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,1);

%%
figure('Color','White')
contourf(Kx,Ky,Omega_k(:,:,4))
%% Berry Curvature Dipole
nk=101;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=1.5;
nsband=1:8;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,0);
%%
Enum=201;
Emin=-1;
Emax=1;
T=100;
[Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband),Enum,Emin,Emax,T);
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
nk=81;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_list=[0.0,0.02,0.04,0.06,0.08,0.1];
%Electric_field_list=[0.0,0.05]
Enum=501;
Emin=-1;
Emax=1;
T=100;
nsband=1:8;
Enk_all=zeros(nk,nk,nbands,length(Electric_field_list));
Omega_dk_all=zeros(nk,nk,size(nsband,2),length(Electric_field_list));
bcd_all=zeros(length(Electric_field_list),Enum);
for i=1:length(Electric_field_list)
    tic;
	Electric_field_in_evpA=Electric_field_list(i);
	[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
	Enk_all(:,:,:,i)=Enk-efermi;
	[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,0);
	Omega_dk_all(:,:,:,i)=Omega_dk;
	[Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband)-efermi,Enum,Emin,Emax,T);
	bcd_all(i,:)=bcd;
    toc;
end



figure('Color','White')
% plot(Eaxis,bcd_all(:,:),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,bcd_all(:,:),'Linestyle','-','LineWidth',2)
legend('E=0.00','E=0.02','E=0.04','E=0.06','E=0.08','E=0.10','Location','northwest','NumColumns',1)
xlabel('E-E_f(eV)')
ylabel('$D_{xz} (\AA)$','Interpreter','latex','FontSize',24)
% yrange=[-6,6]
xrange=[-0.2,0.2]
xlim(xrange)
% ylim(yrange)
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
load("data/TaIrTe4/dos/TaIrTe4_Enk-400x400.mat")
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
