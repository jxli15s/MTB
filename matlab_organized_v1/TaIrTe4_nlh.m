clear;
clear all;
%parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/wannier90_hr_p1.dat','data/TaIrTe4/wannier90_hr_p2.dat');


%% Creat slab ham
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
efermi=7.3742;
%% Berry Curvature Dipole
knum=1
nslab=1;
kxline=[0,1];
kyline=[0,1];
% nsband=1:nbands*nslab;
nsband=1:2;
[Kx,Ky] = g.get_Slab2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);
Enk=Enk-efermi;

[Kx_d,Ky_d] = g.get_Slab2Dkmesh(kxline,kyline,knum);

dkx=g.b2(1,:)./knum;
dkx=norm(dkx);
Kx_d=Kx_d+dkx
[~,Enk_d]=MTB.ham.get_slab_plane_bands(g,Kx_d,Ky_d,nslab);
Enk=Enk-efermi;

[Omega_dk,Kx,Ky]=MTB.ham.get_slab_Berrycurvature_dip2(g,Kx,Ky,Enk,Unk,nsband,nslab,0.001,0);

% Omega_dk_bac=Omega_dk;
%%
load("data/TaIrTe4/nlh/201-201/Dip2-2s-201.mat")
Enum=1000;
Emin=-0.2;
Emax=0.2;
nsband=110:130;
T=100;
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
load("data/TaIrTe4/nlh/101-101/Dip2-2s-101.mat")
Enum=1000;
Emin=-0.2;
Emax=0.2;
nsband=1:176;
T=100;
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
load("data/TaIrTe4/nlh/3s/Dip2-3s-201.mat")
Enum=1000;
Emin=-0.3;
Emax=0.3;
nsband=1:200;
T=60;
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-10,max(bcd)+10];
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%


