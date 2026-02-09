clc;
clear;
g = MTB.geometry("WS2_3s");
g = MTB.read_poscar(g,"data/WS2/3-slab/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/WS2/3-slab/wannier90_hr_p1.dat','data/WS2/3-slab/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;
%%
[nbands,~,nrpts]=size(g.ham);
labels={'N','\Gamma','N1'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]
          };% hkpoints-high symmetry k points

nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"2M-WS2-fplo")
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Check Inversion Symmetry              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1=[0  1
    1  0];
kpoint1=[0.5,0.0,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);%[Energy,Psik,hk]
kpoint2=-kpoint1;
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
P=blkdiag(d1,p1,p2);

d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(3))
P=kron(orbital_site,p_bulk)

c=P*hk1*inv(P)-hk2; % P^{-1}H(k)P=H(-k)
max(c,[],'all')
a=eig(Psik'*P*Psik); %% P\phi=e\phi


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("WS2_3s");
gs = MTB.read_poscar(gs,"data/WS2/3-slab/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/WS2/3-slab/wannier90_hr_p1.dat','data/WS2/3-slab/wannier90_hr_p2.dat');
gs.wpos=[];
gs.wpos=gs.atoms*gs.a;
%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=10;
%
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'};
hkpoints={[-0.5,0.0],...
          [0.0,0.0],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0);
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("WS2_3s");
gs = MTB.read_poscar(gs,"data/WS2/3-slab/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/WS2/3-slab/wannier90_hr_p1.dat','data/WS2/3-slab/wannier90_hr_p2.dat');
gs.wpos=[];
gs.wpos=gs.atoms*gs.a;
%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=10;

[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'};
hkpoints={[0.5,0.0],...
          [0.0,0.0],...
          [0.5,0.0]};% hkpoints-high symmetry k points

Np=3;
omegamin=-1;
omegamax=1;
omeganum=100;
omegas=linspace(omegamin,omegamax,omeganum);
nk=51;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);
%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('heat')); %magma plasma inferno cividis inferno hot heat

shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_r)
colormap(slanCM('ice'))
shading interp
% caxis([1, 50])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_bulk)
colormap(slanCM('heat'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Calculate bands along HSL for slab DW BdG      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("WS2_3s");
gs = MTB.read_poscar(gs,"data/WS2/3-slab/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/WS2/3-slab/wannier90_hr_p1.dat','data/WS2/3-slab/wannier90_hr_p2.dat');
gs.wpos=[];
gs.wpos=gs.atoms*gs.a;
[nbands,~,nrpts]=size(gs.ham);
%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points

nk=31;
nslab=20;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03

[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2,mu,delta)
%%
% save("data/NbSe2/slab_bands_bdg_dw.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
load("data/NbSe2/slab_bands_bdg_dw.mat")
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Calculate bands along HSL for slab bands BdG      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");
MillerIndices=[1,0,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points

nk=51;
nslab=50;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG_v2(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
%%
save("data/NbSe2/slab_bands_bdg.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
% load("data/WTe2/slab_bands_bdg.mat")
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")