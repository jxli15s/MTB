clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplot  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("WS2");
%%
% % filename="data/WS2/fplo3/write_hr/wannier90_test.dat";
% % write_hr(g,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Inversion Symmetry         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];

kpoint=[0.0,0.0,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[0.5,0.0,0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
P=blkdiag(d1,p1,p2);
pos=eye(size(g.wpos,1));
kpoint=kpoint*g.b;
for i=1:size(g.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    fprintf('%.6f',pos(i,i));
end
P=P*pos;
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
a=Psik'*P*Psik; %% P \psi=e \psi    e=\psi.T.conj P \psi
c=eig(a);
%%
a=eig(Psik(:,1:28)'*P*Psik(:,1:28));
sum(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate the band structure           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'N','\Gamma','N1','Y','M','L','N'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.5,0.0,0.5],...
          [0.5,0.5,0.5],...
          [0.5,0.0,0.0]
          };% hkpoints-high symmetry k points
% % labels={'C','\Gamma','C1'}; % labels for k
% % hkpoints={[-0.26856,-0.26856,0.0],...
% %           [0.0,0.0,0.0],...
% %           [-0.26856,-0.26856,0.0],...
% %           [0.33091,0.33091,0.6146],...
% %           [-0.66909,0.66909,0.3854],...
% %           [-0.43765,0.90053,0.3854],...
% %           [0.26856,0.73144,0.0],...
% %           };% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"2M-WS2-fplo")
hold on;
plot(kpath,Energy(29,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(28,:),"Color",'red','LineWidth',2);

% filename="data/WS2/data/final_data/band/WS2_bulk_band.dat";
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Check Time Reversal Symmetry            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(22),i*s2);

h1=T*conj(hk1)*inv(T)-hk2; %% TH=HT T=i*sigma_y*K  THT^-1=H 
max(h1,[],'all')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Check Inversion Symmetry              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1=[0  1
    1  0];
kpoint=[0.5,0.0,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.5,-0.0,-0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
P=blkdiag(d1,p1,p2);
pos=eye(size(g.wpos,1));
kpoint=kpoint*g.b;
for i=1:size(g.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    fprintf('%.6f',pos(i,i))
end
P=P*pos;
c=P*hk1*inv(P)-hk2; % P^{-1}H(k)P=H(-k)
max(c,[],'all')
a=eig(Psik'*P*Psik); %% P\phi=e\phi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Check Inversion Symmetry for slab         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("WS2");
MillerIndices=[1,-1,0];
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

d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab))
P=kron(orbital_site,p_bulk)

c=P*hk1*inv(P)-hk2;
max(c,[],'all')

%%
a=eig(Psik1'*P*Psik1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Check Inversion Symmetry for slab Bdg         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for single TRI point
g=read_fplo("WS2");
%%
MillerIndices=[1,-1,0];
% MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=22
mu=0.0
delta=0.03
s1=[0  1
    1  0];
kpoint=[0.0,0.5];
tic;
[Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
toc;
d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab));
P=kron(orbital_site,p_bulk);
P=[P,zeros(size(P));zeros(size(P)),-P];

a4=eig(Psik(:,1:nslab*nbands)'*P*Psik(:,1:nslab*nbands));
sum(a4)

%% for four TRI points
clc;
clear;
g=read_fplo("WS2");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

nslab=20;
mu=0.0;
delta=0.03;

s1=[0  1
    1  0];

d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab));
P=kron(orbital_site,p_bulk);
P=[P,zeros(size(P));zeros(size(P)),-P]; %construct P operator

kpoints=[0.0,0.0;0.5,0.0;0.0,0.5;0.5,0.5];
Unk = zeros(nslab*nbands*2,nslab*nbands*2,size(kpoints,1));
Pnk = zeros(nslab*nbands,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    tic;
    [~,Psik]=MTB.ham.get_slab_bands_BdG_at_q(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
    toc;
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:nslab*nbands)'*P*Psik(:,1:nslab*nbands));
    Z4(i)=sum(Pnk(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Calculate the wilsonloop of domainwall       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load("./data/WS2/data/wloop/issac/wilsonloop_50s_kx_scal.mat") %1-10
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky_scal1.mat") %1-10
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky.mat") %1-10
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky.mat")
% load("./data/WS2/data/wloop/wilsonloop_50s_ky.mat")
% load("./data/WS2/data/wloop/001/wilsonloop_001_50s_kx.mat")
% load("./data/WS2/data/wloop/issac/new/wilsonloop_50s_kx_scal2.mat")
% load("./data/WS2/data/wloop/issac/new/wilsonloop_50s_kx_scal.mat")
figure('Color','white')
% plot(kx,wx)


% kx=linspace(0,1,knum);

plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
% xticks([0,1/2,1])
% xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Calculate the surface states            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("WS2");


% for i=1:size(g.wpos,1)
%     if g.wpos(i,1)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(1,:);
%     end
%     if g.wpos(i,2)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(2,:);
%     end
%     if g.wpos(i,3)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(3,:);
%     end
% end

MillerIndices=[1,-1,0]; %for a axis
% MillerIndices=[1,1,0]; %for b axis 
% MillerIndices=[0,0,1]; %for c axis 
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=51;
nslab=20;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Surface states by Greenfunction           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=read_fplo("WS2");

for i=1:size(gs.wpos,1)
    if gs.wpos(i,1)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(1,:);
    end
    if gs.wpos(i,2)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(2,:);
    end
    if gs.wpos(i,3)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(3,:);
    end
end
%%
MillerIndices=[0,0,1];
% MillerIndices=[1,-1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
efermi=0; %% set Fermi Level
labels={'C','\Gamma','C'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
%%
Np=1;
omegamin=-0.5;
omegamax=0.5;
omeganum=1000;
omegas=linspace(omegamin,omegamax,omeganum);
nk=101;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);

%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l);
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)
%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_r)
colormap(slanCM('vik'))
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
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")
% load("data/WS2/data/E-u/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/Energy-mus-200s-001_v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001-v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-1000s-001.mat")
 % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Gamma.mat")
  load("data/WS2/data/E-u/viper/Energy-mus-300s-1-10.mat")
   % load("data/WS2/data/E-u/viper/Energy-mus-500s-1-10.mat")
  %%
figure()
% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
mus=linspace(-0.5,0.3,301)
% mus=linspace(-0.5,0.3,201)
x=repmat(mus,50,1)
% plot(x',Energy','or',...
%     'LineWidth',1,...
%     'MarkerSize',4,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r')
plot(x',Energy','o')
hold on;
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
% ylim([-0.1 0.1])
% xlim([])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
%x=repmat(mus,30,1)
%load("Energy-mus.mat")
%plot(x',Energy','*-','Color','red')

filename="data/WS2/data/final_data/E-mus/1-10-300s.dat";
for i=1:size(Energy,1)
    outlist=[mus',Energy(i,:)'];
    writeoutput(filename,outlist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")
% load("data/WS2/data/E-u/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/Energy-mus-200s-001_v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001-v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-1000s-001.mat")
 load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Gamma.mat")
 Energy_G=Energy;
  % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-M-201mus.mat")
% load("data/WS2/data/E-u/viper/Energy-mus-500s-001-X.mat")
load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Y.mat")
Energy_Y=Energy;
figure()
% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
% mus=linspace(-0.5,0.3,501)
mus=linspace(-0.5,0.3,201)
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
%%              Energy-mus-k_{x,y}                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky_scal1.mat")
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
load("data/WS2/data/finite-ky/Energy-mus-500s-001_kx.mat")
% load("data/WS2/data/finite-ky/Energy-mus-500s-001_ky.mat")
absEnergys=abs(Energys)
minEnergys=min(absEnergys)
% minEnergys=reshape(minEnergys(:,:,:),101,101)
 minEnergys=reshape(minEnergys(:,:,:),201,201)
[KX,KY]=meshgrid(ky,mus)
figure
% scatter(KX,KY,50,log(2*minEnergys))
surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
colormap(slanCM('RdBu'))
% clim([-10,-7])
% ylim([-0.2,0.1])
% xlim([0.058,0.1334])
% ylim([-0.2,-0.15])
% xlim([0.043,0.128])
% ylim([-0.049,0.076])
ylabel('$\mu$','FontSize',20,'Interpreter','latex')
xlabel('$k_x$','FontSize',20,'Interpreter','latex')
% colormap(slanCM('heat'))
shading interp
% Energy-mus-500s-0015_ky_scal1.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Wilson loop                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% load("./data/WS2/data/wloop/wilsonloop_50s_kx.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_kx_scal.mat")
% plot(x',Energy','*-','Color','red')
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky_scal1.mat")
load("./data/WS2/data/wloop/001/wilsonloop_001_80s_ky.mat")
figure('Color','white')
% plot(kx,wx)
knum=101;
band1=1;
band2=50*44;

% kx=linspace(0,1,knum);

plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
% xticks([0,1/2,1])
% xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

% filename="data/WS2/data/wloop/001/wilsonloop_001_80s_ky.dat";
% for i=1:size(wx,2)
%     outlist=[kx',wx(:,i)];
%     writeoutput(filename,outlist);
% end



%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/WS2/fplo3/POSCAR");
        pos=textread("data/WS2/fplo3/wpos");
        g.wpos=pos;
        ham=textread("data/WS2/fplo3/mydata-p1");
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


function writeoutput(filename,list)
    file=fopen(filename,'a+');
    cloumns=size(list,2);
    raws=size(list,1);
    % fprintf(file,'Time on %s\n',datetime('today'));
    % fprintf(file,'raws %d cloums %d\n',raws,cloumns);
    for i=1:raws
        for j=1:cloumns
            fprintf(file,'%12.6f \t',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
 end


 function write_hr(g,filename)
    % open a file to write wannier hoppings
    fileID=fopen(filename,'w');
    
    [numBands,~,numRpts]=size(g.ham);

    % Get the current date and time
    currentTime = datetime('now');
    formattedTime = datestr(currentTime, 'mm/dd/yyyy at HH:MM:SS');

    % Write the header with the current data and time to the file
    fprintf(fileID,' write on %s\n', formattedTime);
    % Write the orbital numbers to the file
    fprintf(fileID,'\t %d\n', numBands);
    % Write the sites numbers to the file
    fprintf(fileID,'\t %d\n', numRpts);
    
    % numbers per line
    numsPerLine = 15;
    
    % get the site numbers of hoppings
    len=size(g.hopr,1);

    degeneracy = ones(1,len);

    % Write the data into the file with 15 numbers per line
    for i = 1:numsPerLine:len
        % Determine the index range of data for the current line
        endIdx = min(i+numsPerLine-1,len);
        fprintf(fileID, '%5d', degeneracy(i:endIdx));
        fprintf(fileID, '\n');
    end

    for i = 1:numRpts
        for j = 1:numBands
            for k = 1:numBands
                fprintf(fileID, '%5d %5d %5d %5d %5d %12.6f %12.6f\n',g.hopr(i,:),j,k,real(g.ham(j,k,i)),imag(g.ham(j,k,i)));
            end
        end
    end

    
 end

    
