clear;
clear all;
%p=parpool(8)
%delete(gcp('nocreate'))
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
%%
g.wpos=[];
g.wpos=g.atoms*g.a;

[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;


Electric_field_in_evpA=0.1;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
%% Berry Curvature Dipole
nk=301;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.1;
nsband=1:8;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);

knum=nk;
dkx=g.b(1,:)./knum;
dkx=norm(dkx);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
% Kx_d=Kx_d+10^-3;
Kx_d=Kx_d+dkx;
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);

[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);

%%


%%
Enum=1000;
Emin=-0.5;
Emax=0.5;
T=100;
[Eaxis,bcd]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
% [Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)
% yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%
%% Here You could also load("Dip2-V00-E001.mat")
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.4;
gs=moire_potential(g,gs,Vamp);
nk=301;
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R','\Gamma'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;

Electric_field_in_evpA=0.10;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

%%
load("data/TaIrTe4_2d_tb/omegadata/Dip2-V04-E001-301.mat")

Enum=1000;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=50;
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

%%
% bcd00=bcd;
bcd04=bcd;
Eaxis00=Eaxis;
figure('Color','White')
plot(Eaxis00,bcd00,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
plot(Eaxis,bcd04,'Linestyle','-','Color','red','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,zeros(1,1000),'--k','LineWidth',2)
x=ones(1,11)*0.16
y=-5:5;
plot(x,y,'--k','LineWidth',2)
legend('V=0.0','V=0.4','Location','northwest','NumColumns',1)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-0.8,max(bcd)+0.8]
xrange=[0,0.3]
xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%
function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(1,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,1),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));      
 end
 end