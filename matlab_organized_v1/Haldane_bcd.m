clear;
clear all;
g = MTB.geometry("Haldane");
g = MTB.read_poscar(g,"data/Haldane/POSCAR");
[g.ham,g.hopr] = MTB.read_hr('data/Haldane/Haldane_hr.dat');
 e=1.6*10^-19;
h=6.626*10^-34;
G=25812;%h/e^2(Om*Cm)^-1
coef=1/G*10^8;
g.wpos=[];
g.wpos=g.atoms*g.a;
%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'K','\Gamma','M'}; % labels for k
hkpoints={[0.333333,0.333333,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)

%% Berry Curvature Dipole
nk=500;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.1;
nsband=1:2;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Kx_d=Kx+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
% [Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%%
Enum=100;
Emin=-0.5;
Emax=0.5;
T=30;
[Eaxis,bcd]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
% [Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband),Enum,Emin,Emax,T);
figure('Color','White')
% plot(Eaxis,bcd./norm(g.a(1,:)),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(a_0)$','Interpreter','latex','FontSize',24)
% yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


Eaxis_1=Eaxis;
bcd_1=bcd;
%%
%% Calculate Berry Curvature by Quantum metric
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,1:2,0.01,1);
%% Berry Curvature Dipole
nk=101;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.1;
nsband=1:2;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Kx_d=Kx+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%%
Enum=300;
Emin=-0.5;
Emax=0.5;
T=30;
[Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband),Enum,Emin,Emax,T);
figure('Color','White')
plot(Eaxis,bcd./norm(g.a(1,:)),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(a_0)$','Interpreter','latex','FontSize',24)
% yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


Eaxis_1=Eaxis;
bcd_1=bcd;

%%
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
%% Create supercell
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;

Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
% ylim([-0.5,0.5])

%% Berry Curvature Dipole
nk=201;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.1;
nsband=1:30;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
[Kx_d,Ky_d,Kz_d] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
Kx_d=Kx+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(gs,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
% [Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(gs,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%%
Enum=100;
Emin=-0.5;
Emax=0.5;
T=100;
[Eaxis,bcd]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
% [Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband),Enum,Emin,Emax,T);
figure('Color','White')
% plot(Eaxis,bcd./norm(g.a(1,:)),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(a_0)$','Interpreter','latex','FontSize',24)
% yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


Eaxis_1=Eaxis;
bcd_1=bcd;
%% Calculate Berry Curvature by Quantum metric
nsbands=1:30;
nk=101;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(gs,Kx,Ky,Kz);
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(gs,Kx,Ky,Kz,Enk,Unk,nsbands,0.01,1);
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
 %%
 