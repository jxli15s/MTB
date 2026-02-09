clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
g = MTB.read_poscar(g,"data/TaIrTe4/terminal2/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/terminal2/wannier90_hr_p1.dat','data/TaIrTe4/terminal2/wannier90_hr_p2.dat');
%%
g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=-0.48;
nk=101;
%%
Electric_field_in_evpA=0.00;
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.5,0.3,0.0]
[Energy,Psik1,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.5,-0.3,-0.0]
[Energy,Psik2,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(i*s2,eye(4));
% T=kron(eye(4),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")
%%
b=eig(Psik1(:,1:8)'*T*Psik1(:,1:8));
%%
[v,e]=eig(T);
[e1,ind1]=sort(diag(imag(e)));
v=v(:,ind1);

d1=inv(T)*conj(hk1)*T;
d2=inv(v)*T*v.*1i;
% d2=v*T*inv(v).*1i;
d3=inv(v)*hk1*v;
%%

knum=301;
band1=1;
band2=2;

[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
% [wx,unk]=MTB.ham.get_wilsonloop_time(g,T,knum,band1,band2);

%% Calculate the Dos
nk=301;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0;
% nsband=1:88;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);

plottap=1;

Enum=301;
Emin=-0.25;
Emax=0.25;
eps=0.01;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
sum(Dos)*(Emax-Emin)/Enum
xlim([-0.25,0.25])
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
% ylim([0,100])
%yticklabels({})
legend('DOS')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
print('data/TaIrTe4/dos/TaIrTe4-DOS','-dpng','-r300')


%% Berry Curvature Dipole
nk=81;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.01;
nsband=1:8;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Kx_d=Kx+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%%
Enum=300;
Emin=-0.5;
Emax=0.5;
T=10;
[Eaxis,bcd,~]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband),Enum,Emin,Emax,T);
figure('Color','White')
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
knum=101;
band1=1;
band2=4;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
%% edge states
gslab = MTB.geometry("TaIrTe4");
gslab = MTB.read_poscar(gslab,"data/TaIrTe4_2d_tb/POSCAR");
[gslab.ham,gslab.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

MillerIndices=[0,1,0];
Umatrix = gslab.MillerIndicestoumatrix(MillerIndices);
Urot = gslab.surfab;

[nbands,~,nrpts]=size(gslab.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=101;
nslab=50;
efermi=0.0; %% set Fermi Level
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gslab.ham,gslab.hopr2,nslab,nbands,nrpts,hkpoints,nk,gslab.a2,gslab.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")

%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
g = MTB.read_poscar(g,"data/TaIrTe4/terminal2/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/terminal2/wannier90_hr_p1.dat','data/TaIrTe4/terminal2/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;
g.wpos=zeros(size(g.ham,1),3);

n1=1;
n2=2;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R','\Gamma'}; % labels for k
% hkpoints={[0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points

hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points

labels={'R','Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'-X','\Gamma','X'}; % labels for k
hkpoints={[-0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points

% efermi=0.0141;
efermi=-0.52;
nk=101;

Electric_field_in_evpA=0;
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
knum=51;
band1=1;
band2=60;

[wx1,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
[wx2,unk]=MTB.ham.get_wilsonloop_ky(gs,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Write the wannier90_hr.dat       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename="data/TaIrTe4_2d_tb/wannier90_cdw_hr_0910.dat";
MTB.write_hr(gs,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% filename='band_cdw_V01_E001.dat'
% for iband=1:nbands
%     outlist=[kpath',Energy(iband,:)']
%     writeoutput(filename,outlist)
% end

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.0,0.0,0.0]
[Energy,Psik1,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.0,-0.0,-0.0]
[Energy,Psik2,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
s2=[0  -1i
    1i  0];
s2=[1 0
    0 -1];
T=kron(s2,eye(4));
T=kron(eye(15),T);
% T=kron(eye(4),i*s2);
h1=T*hk1*inv(T)-hk2;
max(h1,[],"all")
%%
b=eig(Psik1(:,1:62)'*T*Psik1(:,1:62));
%%
[v,e]=eig(T);

[e1,ind1]=sort(diag(imag(e)));
v=v(:,ind1);

d1=T*conj(hk1)*inv(T);
d2=inv(v)*T*v.*1j;
d3=inv(v)*hk1*v;
%%

knum=101;
band1=1;
band2=62;

[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%

knum=101;
band1=29;
band2=30;
[wx,unk]=MTB.ham.get_wilsonloop_time(gs,v,knum,band1,band2);

%%
%% Calculate Berry Curvature by Quantum metric
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,1:2,0.01,1);
%% edge states

MillerIndices=[0,1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;

[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
%%
nk=31;
nslab=50;
efermi=0.0; %% set Fermi Level
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")

%% Calculate the Dos
nk=301;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0;
% nsband=1:88;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);

plottap=1;

Enum=201;
Emin=-0.25;
Emax=0.25;
eps=0.01;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
sum(Dos)*(Emax-Emin)/Enum
%%
knum=81;
band1=1;
band2=60;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%% Nonlinear Hall Effect Berry Curvature Dipole
nk=31;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.0;
nsband=1:120;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
[Kx_d,Ky_d,Kz_d] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
Kx_d=Kx+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
%%
[Omega_k,Kx,Ky,Kz]=MTB.ham.get_Berrycurvature_cop(gs,Kx,Ky,Kz,Enk,Unk,1:120,0.001,1)
% [Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(gs,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(gs,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
%%
Enum=100;
Emin=-0.5;
Emax=0.5;
T=100;
[Eaxis,bcd]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
figure('Color','White')
% plot(Eaxis,bcd./norm(gs.a(1,:)),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(a_0)$','Interpreter','latex','FontSize',24)
% yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

Eaxis_2=Eaxis;
bcd_2=bcd;
%%
figure('Color','White')
plot(Eaxis_1,bcd_1./norm(g.b(1,:)),'Linestyle','-','LineWidth',2)
hold on
plot(Eaxis_2,bcd_2./norm(gs.b(1,:)),'Linestyle','-','Color','red','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}$','Interpreter','latex','FontSize',24)
yrange=[min(bcd)-10,max(bcd)+10]
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;

g.wpos=[];
g.wpos=g.atoms*g.a;
orbital_num=[4,4];
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1])
    ]
g.wpos=g.wpos*5

Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
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
Electric_field_list=[0.0,0.01,0.02,0.03,0.04,0.05];
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


function writeoutput(filename,list)
    file=fopen(filename,'a+');
    cloumns=size(list,2);
    raws=size(list,1);
    % fprintf(file,'Time on %s\n',datetime('today'));
    % fprintf(file,'raws %d cloums %d\n',raws,cloumns);
    for i=1:raws
        for j=1:cloumns
            fprintf(file,'%12.6f',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
 end