clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.5,0.5],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;
%%
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(6,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(7,:),"Color",'red','LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the Mirrorx symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint1=[0.0,0.2,0.4]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.0,0.2,0.4]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);

kpoint=kpoint1*g.b;
m1=[1,0,0,0,0;...
    0,-1,0,0,0;...
    0,0,-1,0,0;...
    0,0,0,1,0;...
    0,0,0,0,-1]
% m1=eye(5)
m1=kron(m1,1j*[0,1;1,0])
m2=m1*1;
m3=m1*1;
m4=[-1,0,0;...
    0,-1,0;...
    0,0,1;]
% m4=eye(3)
m4=kron(m4,1j*[0,1;1,0])
m=[m1,zeros(size(m1)),zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
   zeros(size(m1)),m2,zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
   zeros(size(m1)),zeros(size(m1)),m3,zeros(size(m1,1),size(m4,2));...
   zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),m4];
m=blkdiag(m1,m2,m3,m4);
% 
for i =1:size(g.wpos,1)
    m(i,:)=exp(-2j*(g.wpos(i,1))*kpoint(1)).*m(i,:);
     % m(i,:)=exp(-2j*(g.wpos(i,:))*kpoint').*m(i,:);
end

c=m*hk1*inv(m)-hk2;
max(c,[],'all')


[v,e]=eig(m);

d1=inv(v)*hk1*v;

d2=inv(v)*m*v.*1j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=601;
band1=1;
band2=3;
[wx,unk]=MTB.ham.get_wilsonloop_mirror(g,v,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=301;
band1=1;
band2=6;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check the Mirrorxy symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint1=[0.2,0.2,0.4]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
% [Energy,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b)
kpoint2=[0.2,0.2,0.4]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [Energy,Psik,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b)

kpoint=kpoint1*g.b;


m12=[0,0,0,0,0,-1,0,0,0,0;...
     0,0,0,0,0,0,0,0,-1,0;...
     0,0,0,0,0,0,0,-1,0,0;...
     0,0,0,0,0,0,-1,0,0,0;...
     0,0,0,0,0,0,0,0,0,1;...
     -1,0,0,0,0,0,0,0,0,0;...
     0,0,0,-1,0,0,0,0,0,0;...
     0,0,-1,0,0,0,0,0,0,0;...
     0,-1,0,0,0,0,0,0,0,0;...
     0,0,0,0,1,0,0,0,0,0];
m12=kron(m12,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);

m3=[-1,0,0,0,0;...
    0,0,0,-1,0;...
    0,0,-1,0,0;...
    0,-1,0,0,0;...
    0,0,0,0,1];
m3=kron(m3,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);
% m3(1:10,:)=m3(1:10,:).*exp(1j*(g.wpos(21,:)-[g.wpos(21,2),g.wpos(21,1),g.wpos(21,3)])*kpoint');

m4=[0,0,-1;...
    0,-1,0;...
    -1,0,0;];
m4=kron(m4,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0])%.*exp(2j*(g.wpos(32,:))*kpoint');
m4(1:6,:)=m4(1:6,:).*exp(1j*(g.wpos(31,:)-[g.wpos(31,2),g.wpos(31,1),g.wpos(31,3)])*kpoint');

m=blkdiag(m12,m3,m4);


c=inv(m)*hk1*m-hk2;
max(c,[],'all')

%%
[v,e]=eig(m);

d1=inv(v)*hk1*v;

d2=inv(v)*m*v.*1j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=301;
band1=1;
band2=3;
[wx,unk]=MTB.ham.get_wilsonloop_mirror_xy(g,v,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=301;
band1=1;
band2=6;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check the Mirrorxy symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint1=[0.2,0.2,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[0.2,0.2,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);

kpoint=kpoint1*g.b;
m12=[0,0,0,0,0,1,0,0,0,0;...
     0,0,0,0,0,0,0,0,1,0;...
     0,0,0,0,0,0,0,1,0,0;...
     0,0,0,0,0,0,1,0,0,0;...
     0,0,0,0,0,0,0,0,0,-1;...
     1,0,0,0,0,0,0,0,0,0;...
     0,0,0,1,0,0,0,0,0,0;...
     0,0,1,0,0,0,0,0,0,0;...
     0,1,0,0,0,0,0,0,0,0;...
     0,0,0,0,-1,0,0,0,0,0];

% m12=[0,0,0,0,0,1,0,0,0,0;...
%      0,0,0,0,0,0,0,0,-1,0;...
%      0,0,0,0,0,0,0,-1,0,0;...
%      0,0,0,0,0,0,-1,0,0,0;...
%      0,0,0,0,0,0,0,0,0,-1;...
%      1,0,0,0,0,0,0,0,0,0;...
%      0,0,0,-1,0,0,0,0,0,0;...
%      0,0,-1,0,0,0,0,0,0,0;...
%      0,-1,0,0,0,0,0,0,0,0;...
%      0,0,0,0,-1,0,0,0,0,0];
m12=kron(m12,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);
m12(1:10,:)=m12(1:10,:).*exp(-1j*(g.wpos(1,:)-[g.wpos(1,2),g.wpos(1,1),-g.wpos(1,3)])*kpoint');
m12(11:20,:)=m12(11:20,:).*exp(-1j*(g.wpos(11,:)-[g.wpos(11,2),g.wpos(11,1),-g.wpos(11,3)])*kpoint');
m3=[1,0,0,0,0;...
    0,0,0,1,0;...
    0,0,1,0,0;...
    0,1,0,0,0;...
    0,0,0,0,-1];
m3=kron(m3,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);
m3(1:10,:)=m3(1:10,:).*exp(-1j*(g.wpos(21,:)-[g.wpos(21,2),g.wpos(21,1),-g.wpos(11,3)])*kpoint');
m4=[0,0,1;...
    0,1,0;...
    1,0,0;];
m4=kron(m4,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0])%.*exp(2j*(g.wpos(32,:))*kpoint');
m4(1:6,:)=m4(1:6,:).*exp(-1j*(g.wpos(31,:)-[g.wpos(31,2),g.wpos(31,1),-g.wpos(31,3)])*kpoint');
% m4(1:6,:)=m4(1:6,:).*exp(-1j*(g.wpos(31,:)-[-g.wpos(31,2),-g.wpos(31,1),g.wpos(31,3)])*kpoint');
% m=[m1,zeros(size(m1)),zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
%    zeros(size(m1)),m2,zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
%    zeros(size(m1)),zeros(size(m1)),m3,zeros(size(m1,1),size(m4,2));...
%    zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),m4];
m=blkdiag(-m12,-m3,m4);
% 
% for i =1:size(g.wpos,1)
%     m(i,:)=exp(-2j*(g.wpos(i,1))*kpoint(1)).*m(i,:);
%      % m(i,:)=exp(-2j*(g.wpos(i,:))*kpoint').*m(i,:);
% end

c=inv(m)*hk1*m-hk2;
max(c,[],'all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(18),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the 2pi periodictivity        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint=[0.0,0.0,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[1.0,-0.0,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

c=hk1-hk2;
max(c,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check the Inversion symmetry at single point     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint1=[0.1,0.3,0.2]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.1,-0.3,-0.2]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
P=eye(36);
kpoint=kpoint1*g.b;
for i =1:size(g.wpos)
    P(i,i)=exp(-2j*(g.wpos(i,:))*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
% b=eig(Psik(:,1:6)'*P*Psik(:,1:6))
% sum(b(1:6))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get the Parity of bulk TRIP                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint1=[0.5,0.0,0.0]
[Energy,Psik1,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.5,-0.0,-0.0]
[Energy,Psik2,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
P=eye(36);
kpoint=kpoint1*g.b;
for i =1:size(g.wpos)
    P(i,i)=exp(-2j*(g.wpos(i,:))*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
b=eig(Psik1(:,1:6)'*P*Psik1(:,1:6))
sum(b(1:6))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get the Parity of bulk TRIP                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint1=[0.0,0.0,0.0]
[Energy,Psik1,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.0,-0.0,-0.0]
[Energy,Psik2,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
P=eye(36);
kpoint=kpoint1*g.b;
for i =1:size(g.wpos)
    P(i,i)=exp(-2j*(g.wpos(i,:))*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
b=eig(Psik1(:,1:6)'*P*Psik1(:,1:6))
sum(b(1:6))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get the Parity of eight bulk TRIPs              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0;0.0,0.0,0.5;0.5,0.0,0.5;0.0,0.5,0.5;0.5,0.5,0.5];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(6,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    kpoint=kpoint*g.b;
    tic;
    P=eye(36);
    for j =1:size(g.wpos,1)
    % fprintf('%.6f\n',P(i,i))
    P(j,j)=exp(-2j*(g.wpos(j,:))*kpoint');
    if j>30
        P(j,j)=-P(j,j);
    end
    end
    kpoint=kpoints(i,:)
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
    toc;
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:6)'*P*Psik(:,1:6));
    Z4(i)=sum(Pnk(:,i));
end


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                Calculate Bulk Band BdG            %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [nbands,~,nrpts]=size(g.ham);
% labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
% hkpoints={[0.5,0.5,0.5],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.5,0.5,0.0],...
%           [0.0,0.0,0.0]...
%           };% hkpoints-high symmetry k points
% nk=51;
% efermi=0;
% mu=0.0;
% delta=0.03;
%  [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b,mu,delta)
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b,mu,delta);
% 
% MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
% hold on;
% plot(kpath,Energy(6,:),'Color','magenta','LineWidth',2);
% plot(kpath,Energy(7,:),"Color",'red','LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=31;
nslab=50;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
labels={'X','\Gamma','X'};
Np=1;
omegamin=-1;
omegamax=1;
omeganum=100;
omegas=linspace(omegamin,omegamax,omeganum);
nk=31;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(g.ham,g.hopr2,nbands,nrpts,hkpoints,nk,Np,g.a2,g.b2,omegamax,omegamin,omeganum);

%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_l)
colormap(slanCM('inferno'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_r)
colormap(slanCM('hot'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
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
kpoint=[0.0,0.0]
[h00,h01]=MTB.ham.get_slab_h00(g.ham,g.hopr2,nbands,nrpts,kpoint,Np,g.a2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the Inversion symmetry        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.5,0.5,0.5]
[Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

P=eye(36);

kpoint=kpoint*g.b;
for i =1:size(g.wpos)
    % P(i,i)=exp(2j*(g.wpos(i,:)-[-4.8960,4.8960,-4.8960])*kpoint');
    % P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,-4.8960])*kpoint');
    % fprintf('%.6f\n',P(i,i))
    P(i,i)=exp(2j*(g.wpos(i,:))*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end
% P=eye(36);
% P(31:end,:)=-P(31:end,:);
c=diag(Psik'*P*Psik)
sum(c(1:6))
b=eig(Psik(:,1:6)'*P*Psik(:,1:6))
sum(b(1:6))
%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0;0.0,0.0,0.5;0.5,0.0,0.5;0.0,0.5,0.5;0.5,0.5,0.5];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(6,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    kpoint=kpoint*g.b;
    tic;
    for j =1:size(g.wpos)
    % P(i,i)=exp(2j*(g.wpos(i,:)-[-4.8960,4.8960,-4.8960])*kpoint');
    % P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,-4.8960])*kpoint');
    % fprintf('%.6f\n',P(i,i))
    P(j,j)=exp(2j*(g.wpos(i,:))*kpoint');
    if j>30
        P(j,j)=-P(j,j);
    end
    end
    kpoint=kpoints(i,:)
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
    toc;
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:6)'*P*Psik(:,1:6));
    Z4(i)=sum(Pnk(:,i));
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Check the Inversion symmetry for slab       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=21;
s1=[0  1
    1  0];

del_index=find(g.wpos(:,3)<-1);
kpoint1=[0.0,0.5];
kpoint1=kpoint1*g.b2;
hk1=MTB.ham.get_slab_hk_v2(g.ham,g.hopr2,del_index,nslab,nbands,nrpts,kpoint1,g.a2);
% hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=[0.0,-0.5];
kpoint2=kpoint2*g.b2;
hk2=MTB.ham.get_slab_hk_v2(g.ham,g.hopr2,del_index,nslab,nbands,nrpts,kpoint2,g.a2);
% hk2=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk2);
[Energy2,ind]=sort(diag(D));
Psik2=V(:,ind);


n=nslab-2;
kpoint=kpoint1;
phase1=exp(2j*(g.wpos(10,1:2))*kpoint');
phase2=exp(2j*(g.wpos(20,1:2))*kpoint');
phase3=exp(2j*(g.wpos(30,1:2))*kpoint');
phase4=exp(2j*(g.wpos(36,1:2))*kpoint');
% phase1=exp(2j*(g.wpos(10,1:2)-g.wpos(20,1:2))*kpoint');
% phase2=exp(2j*(g.wpos(20,1:2)-g.wpos(20,1:2))*kpoint');
% phase3=exp(2j*(g.wpos(30,1:2)-g.wpos(20,1:2))*kpoint');
% phase4=exp(2j*(g.wpos(36,1:2)-g.wpos(20,1:2))*kpoint');

block1_1=[eye(10)*phase1,zeros(10,nbands-10)];
block1_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block1_3=zeros(20,nbands);
block1_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block1_5=[zeros(6,30),eye(6)*-1*phase4];
block1=[block1_1;block1_2;block1_3;block1_4;block1_5];
%
block2_1=[eye(10)*phase1,zeros(10,nbands-10)];
block2_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block2_3=zeros(nbands,nbands);
block2_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block2_5=[zeros(6,30),eye(6)*-1*phase4];
block2=[block2_1;block2_2;block2_3;block2_4;block2_5];
%
P=zeros(nslab*nbands-16);
P(1:56,size(P,2)-36+1:end)=block1;
%
step1=nbands-size(del_index); %%20
step2=step1+nbands; %%20+36*n
P(20+1:20+72,size(P,2)-36*2+1:size(P,2)-36*2+36)=block2;

for i =1:(nslab-1)/2-2
    P(20+i*36+1:20+i*36+72,size(P,2)-36*(2+i)+1:size(P,2)-36*(2+i)+36)=block2;
end
P=P+P.';

center1=eye(10,10)*phase1;
center2=eye(10,10)*phase2;
center=blkdiag(center1,center2);
P(20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20,20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20)=center;

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
% a=eig(Psik1'*P*Psik1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Check Inversion Symmetry for slab Bdg         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for single TRI point
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

%%
nslab=21;
mu=0.0;
delta=0.03;
s1=[0  1
    1  0];
kpoint=[0.5,0.0];
del_index=find(g.wpos(:,3)<-1);

%
tic;
[Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q_odd(g.ham,g.hopr2,del_index,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
toc;
%

% get the Parity Operator
%%
n=nslab-2;
kpoint=kpoint*g.b2;
phase1=exp(-2j*(g.wpos(10,1:2))*kpoint');
phase2=exp(-2j*(g.wpos(20,1:2))*kpoint');
phase3=exp(-2j*(g.wpos(30,1:2))*kpoint');
phase4=exp(-2j*(g.wpos(36,1:2))*kpoint');
% phase1=exp(2j*(g.wpos(10,1:2)-g.wpos(20,1:2))*kpoint');
% phase2=exp(2j*(g.wpos(20,1:2)-g.wpos(20,1:2))*kpoint');
% phase3=exp(2j*(g.wpos(30,1:2)-g.wpos(20,1:2))*kpoint');
% phase4=exp(2j*(g.wpos(36,1:2)-g.wpos(20,1:2))*kpoint');

block1_1=[eye(10)*phase1,zeros(10,nbands-10)];
block1_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block1_3=zeros(20,nbands);
block1_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block1_5=[zeros(6,30),eye(6)*-1*phase4];
block1=[block1_1;block1_2;block1_3;block1_4;block1_5];
%
block2_1=[eye(10)*phase1,zeros(10,nbands-10)];
block2_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block2_3=zeros(nbands,nbands);
block2_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block2_5=[zeros(6,30),eye(6)*-1*phase4];
block2=[block2_1;block2_2;block2_3;block2_4;block2_5];
%
P=zeros(nslab*nbands-16);
P(1:56,size(P,2)-36+1:end)=block1;
%
step1=nbands-size(del_index); %%20
step2=step1+nbands; %%20+36*n
P(20+1:20+72,size(P,2)-36*2+1:size(P,2)-36*2+36)=block2;

for i =1:(nslab-1)/2-2
    P(20+i*36+1:20+i*36+72,size(P,2)-36*(2+i)+1:size(P,2)-36*(2+i)+36)=block2;
end
P=P+P.';

center1=eye(10,10)*phase1;
center2=eye(10,10)*phase2;
center=blkdiag(center1,center2);
P(20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20,20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20)=center;


P=[P,zeros(size(P));zeros(size(P)),-P];
%
a4=eig(Psik(:,1:size(P)/2)'*P*Psik(:,1:size(P)/2));
sum(a4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Check Inversion Symmetry for slab Bdg         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for single TRI point
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

%%
nslab=101;
mu=0.1;
delta=0.03;
s1=[0  1
    1  0];
kpoint=[0.5,0.0];
del_index=find(g.wpos(:,3)<-1);

kpoints=[0.0,0.0;0.5,0.0;0.0,0.5;0.5,0.5];
Unk = zeros((nslab*nbands-size(del_index,1))*2,nslab*nbands*2-2*size(del_index,1),size(kpoints,1));
Pnk = zeros(nslab*nbands-size(del_index,1),size(kpoints,1));
Pnh0 = zeros(nslab*nbands-size(del_index,1),size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);

parfor i=1:size(kpoints,1)
        kpoint=kpoints(i,:);
        fprintf("Processing on %d K",i);
        tic;
        [Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q_odd(g.ham,g.hopr2,del_index,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
        toc;
        Unk(:,:,i)=Psik;


% get the Parity Operator

kpoint=kpoint*g.b2;
phase1=exp(-2j*(g.wpos(10,1:2))*kpoint');
phase2=exp(-2j*(g.wpos(20,1:2))*kpoint');
phase3=exp(-2j*(g.wpos(30,1:2))*kpoint');
phase4=exp(-2j*(g.wpos(36,1:2))*kpoint');
% phase1=exp(2j*(g.wpos(10,1:2)-g.wpos(20,1:2))*kpoint');
% phase2=exp(2j*(g.wpos(20,1:2)-g.wpos(20,1:2))*kpoint');
% phase3=exp(2j*(g.wpos(30,1:2)-g.wpos(20,1:2))*kpoint');
% phase4=exp(2j*(g.wpos(36,1:2)-g.wpos(20,1:2))*kpoint');

block1_1=[eye(10)*phase1,zeros(10,nbands-10)];
block1_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block1_3=zeros(20,nbands);
block1_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block1_5=[zeros(6,30),eye(6)*-1*phase4];
block1=[block1_1;block1_2;block1_3;block1_4;block1_5];
%
block2_1=[eye(10)*phase1,zeros(10,nbands-10)];
block2_2=[zeros(10,10),eye(10)*phase2,zeros(10,nbands-20)];
block2_3=zeros(nbands,nbands);
block2_4=[zeros(10,20),eye(10)*phase3,zeros(10,6)];
block2_5=[zeros(6,30),eye(6)*-1*phase4];
block2=[block2_1;block2_2;block2_3;block2_4;block2_5];
%
P=zeros(nslab*nbands-16);
P(1:56,size(P,2)-36+1:end)=block1;
%
step1=nbands-size(del_index); %%20
step2=step1+nbands; %%20+36*n
P(20+1:20+72,size(P,2)-36*2+1:size(P,2)-36*2+36)=block2;

for j =1:(nslab-1)/2-2
    P(20+j*36+1:20+j*36+72,size(P,2)-36*(2+j)+1:size(P,2)-36*(2+j)+36)=block2;
end
P=P+P.';

center1=eye(10,10)*phase1;
center2=eye(10,10)*phase2;
center=blkdiag(center1,center2);
P(20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20,20+((nslab-1)/2-1)*nbands+1:20+((nslab-1)/2-1)*nbands+20)=center;


P=[P,zeros(size(P));zeros(size(P)),-P];

Pnk(:,i)=eig(Psik(:,1:size(P)/2)'*P*Psik(:,1:size(P)/2));
Z4(i)=sum(Pnk(:,i));

h0 = diag([ones(1, size(P,1)/2), -1*ones(1, size(P,1)/2)]);
[Psia,Eh0]=eig(h0)
Pnh0(:,i)=eig(Psia(:,1:size(P)/2)'*P*Psia(:,1:size(P)/2));
end

Z4
tic;
% save("parity_50s_001_mu00.mat","Unk","Pnk","Z4","-v7.3");
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     get wave loacalization for DomainWall         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("data/SrSnO/Unk/301_001/Unk_gamma_10821.mat")
%%
%%
figure()
x=kron(diag(ones(1,300)),ones(1,36));
y=a(10841:end,1)
y=x*abs(y.^2);
x=1:300;
plot(x,y)

%%
x=kron(diag(ones(1,300)),ones(1,36));
%%

load("data/SrSnO/Unk/300_001/Unk_10801.mat")
%%

x=kron(diag(ones(1,600)),ones(1,36));
y=x*abs(a.^2);
x=1:600;
plot(x,y)
%%
x=1:21600;
y=abs(a.^2)
plot(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     get wilsonloop for DomainWall                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO")
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham);
%%
nslab=21;
efermi=0.0; %% set Fermi Level
del_index=find(g.wpos(:,3)<-1);
k_mesh_dim=[31,21];
band1=1;
band2=nslab*nbands-size(del_index,1);
mu=0;
delta=0.03;

%%
tic;
%wx=MTB.ham.get_wilsonloop_slab_BdG_v2(g,knum,nslab,band1,band2,mu,delta)
%wx=MTB.ham.get_wilsonloop_slab_BdG_v3(g,k_mesh_dim,nslab,band1,band2,mu,delta)
wx=MTB.ham.get_wilsonloop_slab_BdG_v3_del(g,k_mesh_dim,nslab,band1,band2,del_index,mu,delta);
toc;
%kx=linspace(0,1,k_mesh_dim(1));
%%
% kx=linspace(0,0.12,k_mesh_dim(1));
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx.mat")
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del.mat")
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del_s1.mat")
load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del_s2.mat")
figure('Color','white')
% kx=linspace(0,1,101);
% kx=linspace(0,0.12,201);
plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
% xticks([0,1/2,1])
% xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%
% save("wilsonloop_60s_ky.mat","wx","kx","-v7.3");
%%

% load("data/SrSnO/wilsonloop/60s/Energy-mus-fplo-node2-M-10.mat")
load("data/SrSnO/wilsonloop/60s/Energy-mus-fplo-node2-M-20.mat")
figure()
mus=linspace(-2,2,501);
x=repmat(mus,40,1);
% plot(x',Energy','*-','Color','red')
plot(x',Energy','*-')

%% Supercell calculations
%%%
n1=10;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','X','\Gamma','Y','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=21;
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
MillerIndices=[1,0,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
%%

[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
nk=51;
nslab=1;
efermi=0; %% set Fermi Level

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")
%%


%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/SrSnO/fplo/POSCAR");
        pos=textread("data/SrSnO/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/SrSnO/fplo/mydata-p1");
        orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
        orbital_index=find(orbital_index_logical);
        orbital=ham(orbital_index_logical,1:2);
        orbital_num=sqrt(size(orbital_index,1));

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

%%
function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(2,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,2),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=pi/2+0.01;
       V=Vamp.*(cos(2*pi/7/a*x+phi));      
 end
end




