clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("WTe2");
% g=roate_geometry([1,0,0],g);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Write the wannier90_hr.dat from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename="data/WTe2/fplo/wannier90_hr.dat";
% MTB.write_hr(g,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the 2pi periodictivity        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
kpoint=[0.0,0.0,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[1.0,-0.0,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

c=hk1-hk2;
max(c,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y','M','\Gamma'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5],...
          [0.0,0.5,0.5],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(35,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(36,:),"Color",'red','LineWidth',2);

% filename='./data/SrSnO/data_all/SrSnO_bulk_band.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Check Time Reversal Symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(26),i*s2);

h1=inv(T)*conj(hk1)*T-hk2;
max(h1,[],'all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the Mirrorx symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
%%
kpoint1=[0.0,0.5,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[0.0,-0.5,0.3]
[Energy,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
sy=[0,-1j;1j,0];
s0=[1,0;0,1];
mx_eig=diag([1,-1,1,1,1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1]);
% mx_eig=diag(ones(1,26));
% Here we need consider some atoms has been moved to the outside of the
% cell, a phase should be considered by exp(-1j*k*R)
a=kron(diag([repmat(1,1,4)*exp(1j*kpoint1(2)*2*pi),ones(1,4),repmat(1,1,4)*exp(1j*kpoint1(2)*2*pi),ones(1,4),ones(1,5),repmat(1,1,5)*exp(1j*kpoint1(2)*2*pi)]),ones(2))*exp(-1j*kpoint1(2)*pi);
% a2=kron(diag([repmat(-1,1,4),ones(1,4),repmat(-1,1,4),ones(1,4),ones(1,5),repmat(-1,1,5)]),ones(2))*exp(-1j*kpoint1(2)*pi);
mx=kron(mx_eig,1j*sy).*a;

% kpoint=kpoint1*g.b;
% m1=[1,0,0,0,0;...
%     0,-1,0,0,0;...
%     0,0,-1,0,0;...
%     0,0,0,1,0;...
%     0,0,0,0,-1];
% % m1=eye(5)
% m1=kron(m1,1j*[0,1;1,0]);
% m2=m1*1;
% m3=m1*1;
% m4=[-1,0,0;...
%     0,-1,0;...
%     0,0,1;];
% % m4=eye(3)
% m4=kron(m4,1j*[0,1;1,0]);
% % m=[m1,zeros(size(m1)),zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
% %    zeros(size(m1)),m2,zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
% %    zeros(size(m1)),zeros(size(m1)),m3,zeros(size(m1,1),size(m4,2));...
% %    zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),m4];
% m=blkdiag(m1,m2,m3,m4);
% % 
% for i =1:size(g.wpos,1)
%     m(i,:)=exp(-2j*(g.wpos(i,1))*kpoint(1)).*m(i,:);
%      % m(i,:)=exp(-2j*(g.wpos(i,:))*kpoint').*m(i,:);
% end

c=mx*hk1*inv(mx)-hk1;
max(c,[],'all')
%
% a=Psik'*P*Psik; %% P \psi=e \psi    e=\psi.T.conj P \psi
% a=eig(Psik(:,1:36)'*mx*Psik(:,1:36))
% sum(a)
%%

[v,e]=eig(mx);

d1=inv(v)*hk1*v;

d2=inv(v)*mx*v.*1j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Inversion Symmetry         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);

kpoint=[0.1,0.2,0.3];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.1,-0.2,-0.3];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

space_Te=[0,1,0,0;...
          1,0,0,0;...
          0,0,0,1;...
          0,0,1,0;...
    ];
orbita_Te=diag([1,1,repmat(-1,1,6)]);
P_Te=kron(space_Te,orbita_Te);
space_W=[0,1;...
        1,0];
orbita_W=diag(ones(1,10));
P_W=kron(space_W,orbita_W);
P=blkdiag(P_Te,P_W);
%%
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
a=Psik'*P*Psik; %% P \psi=e \psi    e=\psi.T.conj P \psi
c=eig(a);
%%
a=eig(Psik(:,1:36)'*P*Psik(:,1:36));
sum(a)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Calculate Z2 and Z4                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoints=[0.0,0.0,0.0;0.0,0.5,0.0;0.0,0.0,0.5;0.0,0.5,0.5];
Unk = zeros(nbands,nbands,size(kpoints,1));
nocc=36;
Pnk = zeros(nocc,size(kpoints,1));
Pnh0 = zeros(nocc,size(kpoints,1));
% Pnk = zeros(nbands/2+2,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
Z4_2  = zeros(size(kpoints,1),1);
% PnkZ2 = zeros(1,size(kpoints,1));
Z2 = zeros(size(kpoints,1),1);

space_Te=[0,1,0,0;...
          1,0,0,0;...
          0,0,0,1;...
          0,0,1,0;...
    ];
orbita_Te=diag([1,1,repmat(-1,1,6)]);
P_Te=kron(space_Te,orbita_Te);
space_W=[0,1;...
        1,0];
orbita_W=diag(ones(1,10));
P_W=kron(space_W,orbita_W);
P=blkdiag(P_Te,P_W);

for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
    Unk(:,:,i)=Psik;
    % Pnk(:,i)=eig(Psik(:,1:nbands/2+2)'*P*Psik(:,1:nbands/2+2));
    Pnk(:,i)=eig(Psik(:,1:nocc)'*P*Psik(:,1:nocc));
    % Pnk(:,i)=eig(Psik(:,57:nocc)'*P*Psik(:,57:nocc));
    Z4(i)=sum(Pnk(:,i));
    Z4_2(i)=sum(abs(Pnk(:,i) + 1) <= 0.1)/2;
    Z2(i)=(-1)^(sum(abs(Pnk(:,i) + 1) <= 0.1)/2);
    % PnkZ2=eig(Psik(:,1:2)'*P*Psik(:,1:2))
    % h0 = diag([ones(1, size(P,1)/2), -1*ones(1, size(P,1)/2)]);
    % [Psia,Eh0]=eig(h0);
    % Pnh0(:,i)=eig(Psia(:,1:nocc)'*P*Psia(:,1:nocc));
end
Z4_m2=mod(sum(Z4_2),4);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=read_fplo("WTe2");
%%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=50;
%%
[nbands,~,nrpts]=size(gs.ham);
% labels={'\Gamma','X','\Gamma'};
% hkpoints={[0.0,0.0],...
%           [0.0,0.5],...
%           [0.0,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','Y'};
hkpoints={[0.0,0.5],...
          [0.0,0.0],...
          [0.0,0.5]};% hkpoints-high symmetry k points

nk=31;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("WTe2");
%%
MillerIndices=[0,0,1];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=50;
%%
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','Y'};
hkpoints={[0.0,-0.5],...
          [0.0,0.0],...
          [0.0,0.5]};% hkpoints-high symmetry k points
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
gs=read_fplo("WTe2");
MillerIndices=[0,0,1];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

[nbands,~,nrpts]=size(gs.ham);
labels={'\Gamma','X','\Gamma'};
hkpoints={[0.0,0.0],...
          [0.0,0.5],...
          [0.0,0.0]};% hkpoints-high symmetry k points

Np=1;
omegamin=-1;
omegamax=1;
omeganum=200;
omegas=linspace(omegamin,omegamax,omeganum);
nk=101;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

%%
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','Y'};
hkpoints={[0.0,-0.5],...
          [0.0,0.0],...
          [0.0,0.5]};% hkpoints-high symmetry k points


Np=1;
omegamin=-0.6;
omegamax=0.6;
omeganum=300;
omegas=linspace(omegamin,omegamax,omeganum);
nk=201;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);
%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('heat')); %magma plasma inferno cividis inferno hot heat
caxis([0, 100])
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Check Inversion Symmetry for slab         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=2;
s1=[0  1
    1  0];

kpoint1=[0.5,0.0]
kpoint1=kpoint1*g.b2;
hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=[-0.5,0.0]
kpoint2=kpoint2*g.b2
hk2=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint2,g.a2);
[V,D]=eig(hk2);
[Energy2,ind]=sort(diag(D));
Psik2=V(:,ind);

space_Te=[0,1,0,0;...
          1,0,0,0;...
          0,0,0,1;...
          0,0,1,0;...
    ];
orbita_Te=diag([1,1,repmat(-1,1,6)]);
P_Te=kron(space_Te,orbita_Te);
space_W=[0,1;...
        1,0];
orbita_W=diag(ones(1,10));
P_W=kron(space_W,orbita_W);
p_bulk=blkdiag(P_Te,P_W);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab))
P=kron(orbital_site,p_bulk)

c=P*hk1*inv(P)-hk2;
max(c,[],'all')

%%
a=eig(Psik1'*P*Psik1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Check Inversion Symmetry for slab DW Bdg      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=20
mu=0.43
delta=0.03
s1=[0  1
    1  0];

kpoint1=[0.0,0.0];
tic;
[Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2,g.b2,mu,delta);
toc;

space_Te=[0,1,0,0;...
          1,0,0,0;...
          0,0,0,1;...
          0,0,1,0;...
    ];
orbita_Te=diag([1,1,repmat(-1,1,6)]);
P_Te=kron(space_Te,orbita_Te);
space_W=[0,1;...
        1,0];
orbita_W=diag(ones(1,10));
P_W=kron(space_W,orbita_W);
p_bulk=blkdiag(P_Te,P_W);


%p_bulk=eye(44)
orbital_site=flip(eye(nslab));
P=kron(orbital_site,p_bulk);
P=[P,zeros(size(P));zeros(size(P)),-P];

a4=eig(Psik(:,1:nslab*nbands)'*P*Psik(:,1:nslab*nbands));
sum(a4)
%% for four TRI points
clc;
clear;
g=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=60;
mu=-0.2;
delta=0.03;

space_Te=[0,1,0,0;...
          1,0,0,0;...
          0,0,0,1;...
          0,0,1,0;...
    ];
orbita_Te=diag([1,1,repmat(-1,1,6)]);
P_Te=kron(space_Te,orbita_Te);
space_W=[0,1;...
        1,0];
orbita_W=diag(ones(1,10));
P_W=kron(space_W,orbita_W);
p_bulk=blkdiag(P_Te,P_W);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab));
P=kron(orbital_site,p_bulk);
P=[P,zeros(size(P));zeros(size(P)),-P]; %construct P operator

kpoints=[0.0,0.0;0.5,0.0;0.0,0.5;0.5,0.5];
% Unk = zeros(nslab*nbands*2,nslab*nbands*2,size(kpoints,1));
Pnk = zeros(nslab*nbands,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    tic;
    [~,Psik]=MTB.ham.get_slab_bands_BdG_at_q(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
    toc;
    % Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:nslab*nbands)'*P*Psik(:,1:nslab*nbands));
    Z4(i)=sum(Pnk(:,i));
end
Z4
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Calculate bands along HSL for slab DW BdG      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5000000000],...
          [0.0000000000,0.0000000000],...
          [0.0,0.5]};% hkpoints-high symmetry k points

nk=51;
nslab=50;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
%%
% save("data/WTe2/slab_bands_bdg_dw.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
load("data/WTe2/slab_bands_bdg_dw.mat")
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Calculate bands along HSL for slab bands BdG      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("WTe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5000000000],...
          [0.0000000000,0.0000000000],...
          [0.0,0.5]};% hkpoints-high symmetry k points

nk=51;
nslab=50;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG_v2(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
%%
% save("data/WTe2/slab_bands_bdg.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
load("data/WTe2/slab_bands_bdg.mat")
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")

 % Energy_G=Energy;
  % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-M-201mus.mat")
% load("data/WS2/data/E-u/viper/Energy-mus-500s-001-X.mat")
load("data/WTe2/E-u/Energy-mus-G-001-300s.mat")
Energy_G=Energy;

load("data/WTe2/E-u/Energy-mus-X-001-300s.mat")
Energy_Y=Energy;
figure()

% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
% mus=linspace(-0.5,0.3,501)
% mus=linspace(-0.5,0.3,201)
x=repmat(mus,50,1)
plot(x',Energy_G','or',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
hold on;
plot(x',Energy_Y','ob',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
ylabel('E$_{gap}$','FontSize',24,'Interpreter','latex')
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
% ylim([-0.1 0.1])
% ylim([-0.06,0.06])
ylim([-0.005,0.005])
xlim([-1.5,1.5])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
%x=repmat(mus,30,1)
%load("Energy-mus.mat")
%plot(x',Energy','*-','Color','red')

% 
% filename="data/WS2/data/final_data/E-mus/Gamma.dat";
% for i=1:size(Energy_G,1)
%     outlist=[mus',Energy_G(i,:)'];
%     writeoutput(filename,outlist);
% end
% filename="data/WS2/data/final_data/E-mus/Y.dat";
% for i=1:size(Energy_Y,1)
%     outlist=[mus',Energy_Y(i,:)'];
%     writeoutput(filename,outlist);
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")

 % Energy_G=Energy;
  % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-M-201mus.mat")
% load("data/WS2/data/E-u/viper/Energy-mus-500s-001-X.mat")
% load("data/WTe2/E-u/Energy-mus-G-010-500s.mat")
load("data/WTe2/E-u/Energy-mus-G-010-500s-05ev.mat")

Energy_G=Energy;

% load("data/WTe2/E-u/Energy-mus-Y-010-500s.mat")
load("data/WTe2/E-u/Energy-mus-Y-010-500s-05ev.mat")
Energy_Y=Energy;
figure()

% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
% mus=linspace(-0.5,0.3,501)
% mus=linspace(-0.5,0.3,201)
x=repmat(mus,100,1)
% x=repmat(mus,100,1)
plot(x',Energy_G','or',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
hold on;
plot(x',Energy_Y','ob',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
ylabel('E$_{gap}$','FontSize',24,'Interpreter','latex')
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
ylim([-0.001 0.001])
xlim([-0.8,0.8])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
%x=repmat(mus,30,1)
%load("Energy-mus.mat")
%plot(x',Energy','*-','Color','red')

% 
% filename="data/WS2/data/final_data/E-mus/Gamma.dat";
% for i=1:size(Energy_G,1)
%     outlist=[mus',Energy_G(i,:)'];
%     writeoutput(filename,outlist);
% end
% filename="data/WS2/data/final_data/E-mus/Y.dat";
% for i=1:size(Energy_Y,1)
%     outlist=[mus',Energy_Y(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Functions                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/WTe2/fplo/POSCAR");
        pos=textread("data/WTe2/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/WTe2/fplo/mydata-p1");
        orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
        orbital_index=find(orbital_index_logical);
        orbital=ham(orbital_index_logical,1:2)
        orbital_num=sqrt(size(orbital_index,1))

        ham2=ham;
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
            end
        end
        i=size(orbital_index,1);
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)/g.a;
        toc;

                %
        tic;
        dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
        hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
        [x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                         round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                         round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
        hopr=round([x(:),y(:),z(:)]);
        toc;
        %
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                for j=orbital_index(i)+1:orbital_index(i+1)-1
                    index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
                    hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
                end
            end
        end


        i=size(orbital_index,1)-1
        for j=orbital_index(i+1)+1:size(ham2,1)
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
        end

        toc;
        g.ham=hambac;
        g.hopr=hopr;
        fplo=g;
end