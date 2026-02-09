clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("NbSe2");
% g=roate_geometry([1,0,0],g);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Write the wannier90_hr.dat from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename="data/WTe2/fplo/wannier90_hr.dat";
% MTB.write_hr(g,filename)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','M','K','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.333333,0.333333,0.0],...
          [0.0,0.0,0.0]
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(1,:),'Color','blue','LineWidth',2);
plot(kpath,Energy(2,:),"Color",'red','LineWidth',2);

% filename='./data/SrSnO/data_all/SrSnO_bulk_band.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Check Mirror M_y  for slab DW BdG           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.25,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end
%
[nbands,~,nrpts]=size(gs.ham);
MillerIndices=[0,1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
%
nslab=3;
delta=0.03;% 0.03~0.05
mu=0;
del_index=find(gs.wpos(:,2)>3);
%%
kpoint1=[0.0,0.0];
kpoint1b=kpoint1*gs.b2;
% hk1=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1b,gs.a2);
% hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
hk1=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1b,gs.a2,gs.b2,mu,delta)

[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=kpoint1;
kpoint2b=kpoint2*gs.b2;
% hk2=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2,gs.a2);
hk2=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2b,gs.a2,gs.b2,mu,delta)
% hk2=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk2);
[Energy2,ind]=sort(diag(D));
Psik2=V(:,ind);
%%
C=[0,  1, 0, 0, 0, 0;...
   -1, 0, 0, 0, 0, 0;...
   0, 0, 0, 1, 0, 0;...
   0, 0, -1, 0, 0, 0;...
   0, 0, 0, 0, 0,1;...
   0, 0, 0, 0, -1, 0;]
% C=[0, 1, 0, 0, 0, 0;...
%    -1, 0, 0, 0, 0, 0;...
%    0, 0, 0,exp(1j*2*pi/3), 0, 0;...
%    0, 0, exp(1j*1*pi/3), 0, 0, 0;...
%    0, 0, 0, 0, 0,exp(1j*4*pi/3);...
%    0, 0, 0, 0, exp(1j*5*pi/3), 0;]
n=(nbands*nslab-size(del_index,1))/6;
P=flip(eye(n));
P=kron(P,C);
P=[P,zeros(size(P));zeros(size(P)),-conj(P)];
h=P'*hk1*P-hk2
max(h,[],'all')
%%
[v,e]=eig(P);
[e,ind]=sort(diag(e));
v=v(:,ind);
d1=inv(v)*hk1*v;
d2=inv(v)*P*v.*1j;
hh=inv(v)*hk1*v.*1j
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Get the Zak phase for BdG DW          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.25,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end
%
[nbands,~,nrpts]=size(gs.ham);
MillerIndices=[0,1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
%
nslab=101;
delta=0.03;% 0.03~0.05
mu=0;
del_index=find(gs.wpos(:,2)>3);
%%
kpoint1=[0.5,0.0];
kpoint1b=kpoint1*gs.b2;
% hk1=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1b,gs.a2);
% hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
hk1=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1b,gs.a2,gs.b2,mu,delta);

tic;
[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);
toc;

kpoint2=kpoint1;
kpoint2b=kpoint2*gs.b2;
% hk2=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2,gs.a2);
hk2=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2b,gs.a2,gs.b2,mu,delta);
% hk2=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk2);
[Energy2,ind]=sort(diag(D));
Psik2=V(:,ind);
%%
C=[0,  1, 0, 0, 0, 0;...
   -1, 0, 0, 0, 0, 0;...
   0, 0, 0, 1, 0, 0;...
   0, 0, -1, 0, 0, 0;...
   0, 0, 0, 0, 0,1;...
   0, 0, 0, 0, -1, 0;]
n=(nbands*nslab-size(del_index,1))/6;
P=flip(eye(n));
P=kron(P,C);
My=[P,zeros(size(P));zeros(size(P)),-conj(P)];

h=My'*hk1*My-hk2;
max(h,[],'all')
%%
[v,e]=eig(My);
[e,ind]=sort(diag(e));
v=v(:,ind);
d1=inv(v)*hk1*v;
d2=inv(v)*My*v;
hh=inv(v)*hk1*v;
%% Calculate the Zak Phase
knum=51;
kx=linspace(0,1,knum);
ky=linspace(0,0,knum);
kpoints=[kx(1:knum)',ky(1:knum)'];
band1=1;
band2=(nbands*nslab-size(del_index,1))/2;
tem1=eye(band2-band1+1);
tem2=eye(band2-band1+1);
psik1=zeros(band2*2,band2*2,knum);
psik2=zeros(band2*2,band2*2,knum);
%%
parfor k_idx=1:knum
fprintf("Procesing on k_idx = %d\n",k_idx);

kpoint=kpoints(k_idx,:);
kpointb=kpoint*gs.b2;
hk0=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpointb,gs.a2,gs.b2,mu,delta);
hk=inv(v)*hk0*v;

hk1=hk(1:size(hk,1)/2,1:size(hk,2)/2);
[V1,D1]=eig(hk1);
[e1,ind1]=sort(diag(real(D1)));
psik1(:,:,k_idx)=V1(:,ind1);

hk2=hk(size(hk,1)/2+1:end,size(hk,2)/2+1:end);
[V2,D2]=eig(hk2);
[e2,ind2]=sort(diag(real(D2)));
psik2(:,:,k_idx)=V2(:,ind2);
end

psik1(:,:,knum)=psik1(:,:,1);
psik2(:,:,knum)=psik2(:,:,1);
for k_idx=1:knum-1
    tem1=tem1*psik1(:,band1:band2,k_idx)'*psik1(:,band1:band2,k_idx+1);
    tem2=tem2*psik2(:,band1:band2,k_idx)'*psik2(:,band1:band2,k_idx+1);
end

wx1(:)=sort(angle(eig(tem1))) / pi;
wx2(:)=sort(angle(eig(tem2))) / pi;
% Psik=blkdiag(psik1,psik2);
% E=[e1;e2];
% Psik=blkdiag(psik1,psik2);
% for j=1:knum-1
%     tem=tem*unk1(:,band1:band2,j,i)'*unk1(:,band1:band2,j+1,i); % along ky give us kx evolution
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Get the corner for open xy          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end
%%
%===Construct xy opened model======
%
n1=31;
n2=31;
tic;
gs = MTB.ham.get_xy_open_wannier(gs,n1,n2);
toc;
%
gs.ham=gs.ham{1};
%%
del_index=find(gs.wpos(:,2)<0.0001);
%
gs.ham(del_index,:)=[];
gs.ham(:,del_index)=[];
gs.wpos(del_index,:)=[];
gs.atoms=gs.wpos(1:6:end,:)*inv(gs.a);
%%
del_index=find(gs.wpos(:,1)<0.0001);
%
gs.ham(del_index,:)=[];
gs.ham(:,del_index)=[];
gs.wpos(del_index,:)=[];
gs.atoms=gs.wpos(1:6:end,:)*inv(gs.a);
%%
hk1=gs.ham;
[V, D] = eig(full(1/2*hk1+1/2*hk1'));
[vals,ind]=sort(diag(D));
% sgn = sign(vals);
E=vals;
Psik=V(:,ind);
%%
figure()
plot(E,'ro')
%%
figure()
x=gs.wpos(1:end,1);
y=gs.wpos(1:end,2);
p=abs(Psik).^2;
scatter3(x, y, sum(p(:, 3721:3724), 2), 100, sum(p(:, 3721:3724), 2), 'filled');
% scatter3(x, y, sum(p(:, 582:584), 2), 100, sum(p(:,582:584), 2), 'filled');
% view(2);
%%
figure()
x=gs.wpos(1:end,1);
y=gs.wpos(1:end,2);
p=abs(Psik).^2;
band=3873:3876
scatter3(x, y, sum(p(:, band), 2), 100, sum(p(:,band), 2), 'filled');
%%
for i=3600:4000
    fig = figure('Visible', 'off');
    % p(1:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i);
    % scatter(x(:), y(:), 100, abs(p(:,i)).^2, 'filled');
    scatter3(x(:), y(:),p(:,i), 100, p(:,i), 'filled');
    caxis([0, max(p(:,i))]);
    view(2)
    wxname="data/NbSe2/states/"+int2str(i)+'.png';
    saveas(fig, wxname)
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Get the BdG DW for open xy          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[];
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end

%
n1=81;
n2=151;
gs_xy = MTB.ham.get_xy_open_wannier(gs,n1,n2);
gs_xy.ham=gs_xy.ham{1};

%
del_index=find(gs_xy.wpos(:,2)<0.0001);
gs_xy.ham(del_index,:)=[];
gs_xy.ham(:,del_index)=[];
gs_xy.wpos(del_index,:)=[];
%
del_index=find(gs_xy.wpos(:,1)<0.0001);
gs_xy.ham(del_index,:)=[];
gs_xy.ham(:,del_index)=[];
gs_xy.wpos(del_index,:)=[];
gs_xy.atoms=gs_xy.wpos(1:6:end,:)*inv(gs_xy.a);


figure;
plot(gs_xy.atoms(:,1),gs_xy.atoms(:,2),'o')
%
ypos = gs_xy.wpos(:,2);  % 取出 y 坐标
ymean = mean(ypos);   % 计算平均值
val = zeros(size(ypos));     % 初始化
val(ypos < ymean-1) = 1;       % 小于平均值 → 赋值为 1
val(ypos > ymean+1) = -1;      % 大于平均值 → 赋值为 -1
val=val(1:2:end);
% val = ones(size(ypos,1)/2,1)
%%
% % yfrac = gs_xy.atoms(:,2);  % 取出 y 坐标
% % ymean = mean(yfrac);   % 计算平均值
% % val = zeros(size(yfrac));     % 初始化
% % val(yfrac < ymean-1e-6) = 1;       % 小于平均值 → 赋值为 1
% % val(yfrac > ymean+1e-6) = -1;      % 大于平均值 → 赋值为 -1
%%
mu=0;
delta=0.03;
numEigs=50;
sigma_y=sparse([0, -1j; 1j, 0]);
orbital=spdiags(val, 0, length(val), length(val));
% orbital=kron(speye(3),orbital);
h_delta=kron(orbital,1j*sigma_y*delta);

nbands=size(gs_xy.ham,1);
h_onsite=speye(nbands)*mu;
hk_e=gs_xy.ham-h_onsite;
hk_h=h_onsite-gs_xy.ham.';
hk=[hk_e,h_delta;h_delta', hk_h];
hk=(hk+hk')/2;
% ishermitian(hk)
%%
fprintf("Calculate begin %d\n",1)
tic;
% numEigs=50;
% hk=gs_xy.ham{1}/2+gs_xy.ham{1}'/2;
% vals=sort(eigs(hk,numEigs,'smallestabs'));
[V, D] = eigs(hk, numEigs, 'smallestabs');
vals = diag(D);
[vals, idx] = sort(vals, 'ascend');  % 升序排序本征值
V = V(:, idx);                       % 对应地重新排序本征矢量
toc;
fprintf("Calculate end %d\n",1)
%%
figure
plot(real(vals),'ro')
%%
figure()
x=gs_xy.wpos(1:end,1);
y=gs_xy.wpos(1:end,2);
p=abs(V).^2;
p=p(1:end/2,:)+p(end/2+1:end,:)
band=24:27
scatter3(x, y, sum(p(:, band), 2), 100, sum(p(:,band), 2), 'filled');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Calculate bands along HSL for slab DW BdG      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");
MillerIndices=[0,1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points

nk=51;
nslab=100;
efermi=0.0; %% set Fermi Level
mu=-0.0
delta=0.03

[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2,mu,delta)
%%
save("data/NbSe2/slab_bands_bdg_dw.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
% load("data/NbSe2/slab_bands_bdg_dw.mat")
%%
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
% save("data/NbSe2/slab_bands_bdg.mat","Energy","nk","nslab","nbands","efermi","kpath","labels","kindex","mu","delta","-v7.3");
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
load("data/NbSe2/E-u/Energy-mus-G-100-300s.mat")
Energy_G=Energy;

load("data/NbSe2/E-u/Energy-mus-X-100-300s.mat")
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
load("data/WTe2/E-u/Energy-mus-G-010-500s.mat")
Energy_G=Energy;

load("data/WTe2/E-u/Energy-mus-Y-010-500s.mat")
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
%%

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','M','K','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.333333,0.333333,0.0],...
          [0.0,0.0,0.0]
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")

%%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

%%
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','Y','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.0,0.000000],...
          [0.000000000,0.0000000000,0.000000],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]...
          [0.000000000,0.0000000000,0.000000]};% hkpoints-high symmetry k points
% labels={'\Gamma','M','K','\Gamma'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.333333,0.333333,0.0],...
%           [0.0,0.0,0.0]
%           };% hkpoints-high symmetry k points
nk=201;
% efermi=-1.2533;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)
% plot_geometry_sub(gs,g)

%%


%%
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
% shift=[0.0,0,0] %%
shift=[0.25,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end

MillerIndices=[0,1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points

nk=51;
nslab=101;
efermi=0.0; %% set Fermi Level
mu=0.0
delta=0.03
del_index=find(gs.wpos(:,2)>3);
% del_index=[]
kpoint1=[0.5,0.0];
kpoint1=kpoint1*gs.b2;
numEigs=100;
% [Energy1,~]=MTB.ham.get_slab_bands_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1,gs.a2,gs.b2,mu,delta);
% Energy2=MTB.ham.get_slab_q_sparse_BdG_del(gs.ham,gs.hopr2,del_index,nslab,nbands,numEigs,nrpts,kpoint1,gs.a2,gs.b2,mu,delta);
%%
% [Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
[Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG_del(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2,mu,delta);
% [Energy,kpath,kindex]=MTB.ham.get_slab_bands_BdG(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2,mu,delta);
MTB.plot.plot_bands(Energy,size(Energy,1),efermi,kpath,labels,kindex,"SrSnO-100slab")

%%
%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.25,0,0] %%
% shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end
%
MillerIndices=[0,1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
%
nslab=201;
delta=0.03;% 0.03~0.05
% del_index=find(gs.wpos(:,2)<0.1);
del_index=find(gs.wpos(:,2)>3);
%%
% del_index=[]
numEigs=50;
mus=linspace(-1,1,301);%mu for E_f-mu to E_f+mu of 100 points
kpoint=[0.0,0.0]*gs.b2;% Gamma Point
tic;
% Energy_G=MTB.ham.get_slab_mu_E_sparse_BdG(gs.ham,gs.hopr2,nslab,nbands,numEigs,nrpts,kpoint,gs.a2,mus,delta);
Energy_G=MTB.ham.get_slab_mu_E_sparse_BdG_del(gs.ham,gs.hopr2,del_index,nslab,nbands,numEigs,nrpts,kpoint,gs.a2,gs.b2,mus,delta);
%save("Energy-mus-X-100-300s.mat","Energy","mus");
toc

tic;
kpoint=[0.5,0.0]*gs.b2;% Gamma Point
% Energy_Y=MTB.ham.get_slab_mu_E_sparse_BdG(gs.ham,gs.hopr2,nslab,nbands,numEigs,nrpts,kpoint,gs.a2,mus,delta);
Energy_Y=MTB.ham.get_slab_mu_E_sparse_BdG_del(gs.ham,gs.hopr2,del_index,nslab,nbands,numEigs,nrpts,kpoint,gs.a2,gs.b2,mus,delta);
%save("Energy-mus-X-100-300s.mat","Energy","mus");
toc

% Energy=get_slab_mu_E_sparse_BdG_del(hamiltonian,hopping_r,del_index,nslab,nbands,numEigs,nrpts,kpoint,a,b,mus,delta)
%%
% load("data/NbSe2/E-u/Energy-mus-G-100-300s.mat")
% Energy_G=Energy;

% load("data/NbSe2/E-u/Energy-mus-X-100-300s.mat")
% Energy_Y=Energy;
figure()

% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
% mus=linspace(-0.5,0.3,501)
% mus=linspace(-0.5,0.3,201)
x=repmat(mus,50,1);
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
xlim([-0,0.6])
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi Surface-1             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate plane bands
%%
clc;
clear;
g=read_fplo("NbSe2");
knum=501;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Write the plane eigenvalue               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Occ=2;
% filename='NbSe2_1sEnk-501x501.dat';
% writeEnk(Enk,Kx,Ky,Occ,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plot the fermi surface by contour3          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设你有以下数据：
% E: nkx × nky × nbands 的能量数据
% kx_list, ky_list: 分别是 nkx × 1 和 nky × 1 的向量
% load('E_data.mat'); % 或你已经在工作区中

band_index = 1;  % 要绘制的能带索引
% Ef = 0;          % 费米能级（假设为0）

% 构造网格
% [kx, ky] = meshgrid(kx_list, ky_list);      % 注意 meshgrid 的顺序
kx=Kx;
ky=Ky;
% 获取目标能带的能量值（转置以匹配 meshgrid）
Ez = squeeze(Enk(:,:,band_index))';          

% 绘图
figure;
hold on;

% 1. 三维能带
% s = surf(kx, ky, Ez, 'EdgeColor', 'none');
% colormap turbo
shading interp
alpha(0.9);
Ez1 = squeeze(Enk(:,:,1))';
% surf(kx, ky, Ez1, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
Ez2 = squeeze(Enk(:,:,2))';
% surf(kx, ky, Ez2, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
lighting phong    % 或 gouraud，phong 更平滑
% shading interp    % 表面平滑
% for band = 1:2
%     Ez = squeeze(Enk(:,:,band))';
%     surf(kx, ky, Ez, 'EdgeColor', 'none', 'FaceAlpha', 0.6,'FaceColor', [1,0,0]);
% end


% Ef=0.0;
% 2. 费米面等高线（Ef）
% contour3(kx, ky, Ez, [Ef Ef], 'k', 'LineWidth', 2);  % Fermi contour

Ef_list = -0.2:0.2:0.6;  % 多个等高值
% Ef_list =[-0.2 -0.1]
Ef_list=[-0.089 -0.089]
% 2. 费米面等高线（Ef）
contour3(kx, ky, Ez1, Ef_list, 'LineColor', 'blue', 'LineWidth', 3);
contour3(kx, ky, Ez2, Ef_list, 'LineColor', 'red' ,'LineWidth', 1.5);

% 3. 视图与标签
view(3);
xlabel('k_x'); ylabel('k_y'); zlabel('E(k)');
title(['Bands and Fermi Surface ']);
colorbar;
% axis tight;
axis equal;
box on
view(0,90)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi by Green function     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the fermi surface by green function
%% G(k,mu)=(mu+i\delta-H)^{-1}
%% G(k,mu)=(mu+i\delta-E)^{-1}
%% A(k,mu)=-1\pi*[TrImg(G)]=sum
%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[1,0,0;0,1,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end

knum=501;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(gs,Kx,Ky,Kz);
mu=0.5;
eta=10^-2;
A = compute_spectral_function_kplane(Enk, mu, eta);
%% plot by pcolor
 figure()
 imagesc(A'); axis equal; axis off; colormap hot; colorbar;
%% plot by imagesc
 figure()
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
pcolor(Kx, Ky, A');  % A' 是因为 MATLAB 的列主序
shading interp;      % 平滑显示颜色块
colormap hot;        % 色图
colorbar;
axis equal tight;
title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu)]);
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('$k_y$', 'Interpreter', 'latex');
%% plot by surf
 figure()
surf(Kx, Ky, A', 'EdgeColor', 'none');
% colormap turbo;
colorbar;
view([30, 45]);       % 视角角度
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('$k_y$', 'Interpreter', 'latex');
zlabel('$A(k, \mu)$', 'Interpreter', 'latex');
title(['Spectral Function A(k, \mu = ', num2str(mu), ')']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculate the Fermi by Green function     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
clc;
clear;
g=read_fplo("NbSe2");

% g=roate_geometry([1,0,0],g);
g.atoms=[0,0,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx
%
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");

gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[]
for i=1:length(gs.orbnum_list)
    wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
    gs.wpos=[gs.wpos;wpos];
end
%%
knum=2001;
% % kxline=linspace(-0.5,0.5,knum);
% % kyline=linspace(-0.5,0.5,knum);
% % theta=pi/6;
% % R=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
% % kxline=linspace(-0.5,0.5,knum)*sqrt(3);
% % kyline=linspace(0,0,knum);
% kxline=linspace(0,0,knum);
% kyline=linspace(-0.5,0.5,knum);
% 
% kxline=linspace(-0.5,0.5,knum);
% kyline=linspace(0.5,0.5,knum);

kxline=linspace(0.5,0.5,knum);
kyline=linspace(-0.5,0.5,knum);

kline=[kxline',kyline',zeros(knum,1)];

kline=kline*gs.b;

% kline=(R*kline')'
% kline=kline*g.b;

mus=linspace(-0.4,0.2,1001);
% Enk=MTB.ham.get_kline_bands(gs,kline);
Enk=MTB.ham.get_kline_bands(g,kline);
eta=10^-2;
%
A = compute_spectral_function_kline_mu(Enk, mus, eta);

%plot by pcolor
 figure()
 imagesc(A); axis equal; axis off; colormap hot; colorbar;
 
%% plot by imagesc
 figure()
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
 % [Kx,Ky]=meshgrid(mus,kxline)
 [Kx,Ky]=meshgrid(mus,kyline)
pcolor(Kx, Ky,A);  % A' 是因为 MATLAB 的列主序
shading interp;      % 平滑显示颜色块
colormap hot;        % 色图
colorbar;
% axis equal tight;
% title(['Spectral Function A(k_x, k_y) at \mu = ', num2str(mu)]);
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('$k_y$', 'Interpreter', 'latex');
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
        g = MTB.read_poscar(g,"data/NbSe2/fplo/POSCAR");
        pos=textread("data/NbSe2/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/NbSe2/fplo/mydata-p1");
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

function C = construct_rotation_operator_periodic(wpos, a, b, R, k)
% wpos: Norb×3，轨道位置（笛卡尔坐标）
% a: 3×3 实空间晶格基矢
% R: 3×3 旋转矩阵
% k: 1×3 动量点
Norb = size(wpos, 1);
C = zeros(Norb);

% 平移搜索范围：-1, 0, 1
range = -2:2;
[X, Y, Z] = ndgrid(range, range, 0);  % 平移方向：xy方向±1
shifts = [X(:), Y(:), Z(:)];  % 共 3×3×1 = 9 个平移
Tlist = shifts * a;  % 将分数坐标转为笛卡尔坐标

wpos_rot = (R * wpos')';  % 所有轨道旋转后的位置
k=k*b;
for i = 1:Norb
    ri_rot = wpos_rot(i, 1:3);  % 旋转后的轨道位置

    min_dist = inf;
    j_best = -1;
    T_best = [0 0 0];

    for j = 1:Norb
        for t = 1:size(Tlist,1)
            T = Tlist(t,:);
            diff = wpos(j,:) - (ri_rot + T);
            d = norm(diff);
            if d < min_dist
                min_dist = d;
                j_best = j;
                T_best = T;
            end
        end
    end

    if min_dist > 1e-3
        error('No matching site found for orbital %d under rotation.', i);
    end


    % R=round((ri_rot-wpos(j_best,:))*inv(a));
    % phase = exp(1j * dot(k, R));
    % C(j_best, i) = phase;
    phase = exp(1j * dot(k, wpos(j_best,:) - ri_rot));
    C(j_best, i) = phase;

end
end



function R = rotation_matrix(theta)
% 绕 z 轴旋转 theta 弧度
R = eye(3);
R(1:2,1:2) = [cos(theta), -sin(theta); sin(theta), cos(theta)];
end

%% Enk(nk,nk,nband) Kx(nk,nk) Ky(nk,nk)
function writeEnk(Enk,Kx,Ky,Occ,filename)
    file=fopen(filename, 'w');
    fprintf(file, '%s\n', '# kx     ky     kz     Ev5     Ev4     Ev3     Ev2     Ev1     Ec1     Ec2     Ec3     Ec4     Ec5');
    for i=1:size(Enk,1)
        for j=1:size(Enk,2)
            fprintf(file, [repmat('%12.6f',1,9) '\n'],...
                Kx(i,j), Ky(i,j), 0.0, ...
                Enk(i,j,Occ-1), Enk(i,j,Occ), ...
                Enk(i,j,Occ+1), Enk(i,j,Occ+2), Enk(i,j,Occ+3), Enk(i,j,Occ+4));
        end
        fprintf(file,'\n');
    end
    fclose(file);
end



function A = compute_spectral_function_kplane(Enk, mu, eta)
% 输入参数
% Enk: [nkx, nky, nbands] 的能带数组
% mu: 能量点
% eta: Lorentzian 展宽参数，例如 eta = 0.01

    [nkx, nky, nbands] = size(Enk);
    A = zeros(nkx, nky);  % 初始化谱函数

    for ix = 1:nkx
        for iy = 1:nky
            for n = 1:nbands
                E = Enk(ix, iy, n);
                A(ix, iy) = A(ix, iy) + eta / ( (mu - E)^2 + eta^2 );
            end
        end
    end

    A = A / pi;  % 加上 prefactor 1/pi
end

function A = compute_spectral_function_kline_mu(Enk, mus, eta)
% 输入参数
% Enk: [nkx, nky, nbands] 的能带数组
% mu: 能量点
% eta: Lorentzian 展宽参数，例如 eta = 0.01

    [nkx, nbands] = size(Enk);
    nmu=length(mus);
    A = zeros(nkx, nmu);  % 初始化谱函数
    
    for ix = 1:nkx
        for imu = 1:nmu
            for n = 1:nbands
                E = Enk(ix, n);
                mu=mus(imu);
                A(ix, imu) = A(ix, imu) + eta / ( (mu - E)^2 + eta^2 );
            end
        end
    end

    A = A / pi;  % 加上 prefactor 1/pi
end


