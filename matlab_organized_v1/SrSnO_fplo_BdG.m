clc;
clear;
g = MTB.geometry("SrSnO");
g = MTB.read_poscar(g,"data/SrSnO/fplo/POSCAR");
%%
ham=[];
% fid=fopen('data/SrSnO/fplo/mydata-p1');
% tline=fgetl(fid);
% while ischar(tline)
%     disp(tline);
%     tline=fgetl(fid);
% end
pos=textread("data/SrSnO/fplo/wpos");
g.wpos=pos;
ham=textread("data/SrSnO/fplo/mydata-p1");
orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
orbital_index=find(orbital_index_logical);
orbital=ham(orbital_index_logical,1:2);
orbital_num=sqrt(size(orbital_index,1));
%%
ham2=ham;

for i=1:size(orbital_index,1)-1
    if (orbital_index(i+1)-orbital_index(i))>2
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)-(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
    elseif i==size(orbital_index,1)-1
        ham2(orbital_index(i+1)+1:end,1:3)=ham2(orbital_index(i+1)+1:end,1:3)-(g.wpos(orbital(i+1,1),:)-g.wpos(orbital(i+1,2),:));
        ham2(orbital_index(i+1)+1:end,1:3)=ham2(orbital_index(i+1)+1:end,1:3)/g.a;
    end
end


dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
[x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                 round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                 round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
hopr=round([x(:),y(:),z(:)]);

for i=1:size(orbital_index,1)-1
    if (orbital_index(i+1)-orbital_index(i))>2
        for j=orbital_index(i)+1:orbital_index(i+1)-1
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
        end
    elseif i==size(orbital_index,1)-1
        for j=orbital_index(i+1)+1:size(ham2,1)
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
        end
    end
end

g.ham=hambac;
g.hopr=hopr;

%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.5,0.5],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0; %% set Fermi Level

%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")

%%
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=100;
delta=0.03;% 0.03~0.05
numEigs=30;
% efermi=3.6616
mus=linspace(-0.3,0.3,101);%mu for E_f-mu to E_f+mu of 100 points
kpoint=[0.0,0.0];% Gamma Point
%%
tic;
Energy=MTB.ham.get_slab_mu_E_sparse_BdG(g.ham,g.hopr2,nslab,nbands,numEigs,nrpts,kpoint,g.a2,mus,delta);
save("Energy-mus.mat","Energy");
toc
%%
mus=linspace(-0.3,0.3,101);
x=repmat(mus,30,1);
load("Energy-mus.mat")
% plot(x',Energy','*-','Color','red')
plot(x',Energy','*')
%%
load("data/SrSnO/fplo/800s/Energy-mus-fplo-node2.mat")
%%
figure()
mus=linspace(-0.3,0.3,501);
x=repmat(mus,30,1);
plot(x',Energy','o')
%%
load("data/SrSnO/1000s/isaac/001/Energy-mus-1000s-001.mat")
% load("data/SrSnO/1000s/isaac/001/Energy-mus-500s-001.mat")
%%
figure()
mus=linspace(-0.5,0.5,501); 
x=repmat(mus,50,1);
plot(x',Energy','.','Color','#007EC9')
ylim([-0.0104,0.0104])
ylabel('Energy(eV)')
xlabel('mu')
title('mu=0.01 layers=1000')
%%
load("data/SrSnO/1000s/isaac/003/Energy-mus-500s-001.mat")
%%
figure()
mus=linspace(-0.3,0.3,201);
x=repmat(mus,50,1);
plot(x',Energy','o','Color','#007EC9')
ylim([-0.031,0.031])
ylabel('Energy(eV)')
xlabel('mu')
title('mu=0.01 layers=500')
%%

[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.333333,0.333333]};% hkpoints-high symmetry k points
labels={'X','\Gamma'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
nk=51;
nslab=100;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03

[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG_v2(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")


%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
nk=11;
nslab=100;
efermi=0.0; %% set Fermi Level
mu=-0.0888
delta=0.03

[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")
ylim([-0.4,0.4])

%%
%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
nk=51;
nslab=200;
numEigs=50;
efermi=0.0; %% set Fermi Level
mu=0;
delta=0.03;
[Energy,kpath,kk]=MTB.ham.get_slab_bands_sparse_BdG(g.ham,g.hopr2,nslab,nbands,numEigs,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta);

%%
MTB.plot.plot_bands(Energy,numEigs,efermi,kpath,labels,kk,"SrSnO-100slab")
ylim([-0.4,0.4])
%%
load("data/SrSnO/fplo/800s/band/Energy-band-0088.mat")
labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.333333,0.333333]};% hkpoints-high symmetry k points
nk=301;
efermi=0;
numEigs=50;
MTB.plot.plot_bands(Energy,numEigs,efermi,kpath,labels,kindex,"SrSnO-100slab")
ylim([-0.4,0.4])

%%
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'o','Color','black');
    hold on
end
mus=[-0.16,-0.13,-0.11,-0.06,-0.03,0.0];
%%
load("data/SrSnO/fplo/800s/band/Energy-band-6.mat")
% load("data/SrSnO/fplo/800s/band/Energy-band-0088.mat")
labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
%%
nk=401;
efermi=0;
numEigs=50;
MTB.plot.plot_bands(Energy,numEigs,efermi,kpath,labels,kindex,"SrSnO-100slab")
ylim([-0.03,0.03])
xlim([0.2,0.4])
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'o','Color','black');
    hold on
end
%%
load("data/SrSnO/fplo/gamma/Energy-q.mat")
figure;
nslab=100:50:1000
plot(nslab,Energy(26,:)-Energy(25,:))

%%
% zshift=g.a*inv(Urot);

zvalue=zeros(nbands*nslab*2,3);
zpos=1:nslab*2;
zpos=zpos-nslab;
zpos=kron(zpos,ones(1,nbands));
wpos=repmat(g.wpos,nslab*2,1);

wpos(:,3)=wpos(:,3)+zshift(3,3)*zpos';
g.wpos=wpos;
%%

knum=11;
nslab=20;
band1=1;
band2=nslab*nbands;
mu=0;
delta=0.03;
%%
tic;
[wx,~]=MTB.ham.get_wilsonloop_slab_BdG(g,knum,nslab,band1,band2,mu,delta);
toc;
%%
plot(wx)
%%
%%Willson loop
wx=MTB.ham.get_wilsonloop_slab_BdG_v2(g,knum,nslab,band1,band2,mu,delta)
%%
load("data/SrSnO/wilsonloop/80s/wilsonloop_80s_kx.mat")
% load("data/SrSnO/wilsonloop/80s/wilsonloop_80s_kx.mat")
figure('Color','white')
plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
xticks([0,1/2,1])
xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)