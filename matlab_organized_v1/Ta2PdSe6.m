clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("NbSe2");
g.atoms=g.wpos*inv(g.a);
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
kpoint=[0.0,0.0,0.0];
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[1.0,-0.0,-0.0];
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
c=hk1-hk2;
max(c,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'L0','\Gamma','Y','L1','V'}; % labels for k
hkpoints={[-0.256553,-0.256553,0.0],...
          [0.0,0.0,0.0],...
          [0.5,-0.5,0.0],...
          [0.256553,-0.743447,0.0],...
          [0,-0.5,0],...
          [-0.256553,-0.256553,0.0]
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk");
hold on;
plot(kpath,Energy(57,:),'Color','blue','LineWidth',2);
plot(kpath,Energy(56,:),"Color",'red','LineWidth',2);

% filename='./data/SrSnO/data_all/SrSnO_bulk_band.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
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
g=read_fplo("Ta2PdSe6");
knum=101;
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

band_index = 56;  % 要绘制的能带索引
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
Ez1 = squeeze(Enk(:,:,56))';
surf(kx, ky, Ez1, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
Ez2 = squeeze(Enk(:,:,57))';
surf(kx, ky, Ez2, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("Ta2PdSe6");

MillerIndices=[1,-1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=51;
%%
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'};
hkpoints={[-0.5,0.0],...
          [0.0,0.0],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("Ta2PdSe6");
MillerIndices=[1,-1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
%
[nbands,~,nrpts]=size(gs.ham);
labels={'\Gamma','X','\Gamma'};
hkpoints={[0.0,0.0],...
          [0.5,0.0],...
          [0.0,0.0]};% hkpoints-high symmetry k points
%%
Np=1;
omegamin=-0.2;
omegamax=0.2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Check Time Reversal Symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.2,0.3,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.3,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(3),i*s2);

h1=inv(T)*conj(hk1)*T-hk2;
max(h1,[],'all')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check C3z Rotation Symmetry for        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CH(k)C^{}
theta_deg=120;
sigma_z=[1,0;0,1];
theta = theta_deg * pi / 180;
R = rotation_matrix(theta);
kpoint1=[0,0,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[exp(-1j*2*pi/3*-3/2),0,0,0,0,0;...
   0,exp(-1j*2*pi/3*3/2),0,0,0,0;...
   0,0,exp(-1j*2*pi/3*-5/2),0,0,0;...
   0,0,0,exp(-1j*2*pi/3*5/2),0,0;...
   0,0,0,0,exp(-1j*2*pi/3*-1/2),0;...
   0,0,0,0,0,exp(-1j*2*pi/3*1/2)]
%%
kpoint=[0.333333,0.333333,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C*hk1*C'-hk2;
%%
cn=eig(Psik(:,1:2)'*C*Psik(:,1:2))
Jz=angle(cn)/theta%-angle(exp(1j*pi*sum(kpoint)))/theta

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Check C2(2x+y)                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_deg=180;
sigma_z=[1,0;0,1];
theta = theta_deg * pi / 180;
R=[1,0,0;0,-1,0;0,0,-1];
kpoint1=[0,0,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[0,  1, 0, 0, 0, 0;...
   1, 0, 0, 0, 0, 0;...
   0, 0, 0,-1, 0, 0;...
   0, 0, -1, 0, 0, 0;...
   0, 0, 0, 0, 0,-1;...
   0, 0, 0, 0, -1, 0;];
%%
kpoint=[0.23,0.45,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C*hk1*C'-hk2

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Check C2(x+2y)                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_deg=180;
sigma_z=[1,0;0,1];
theta = theta_deg * pi / 180;
R=[1,0,0;0,-1,0;0,0,-1]
kpoint1=[0,0,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[0,  1, 0, 0, 0, 0;...
   1, 0, 0, 0, 0, 0;...
   0, 0, 0,-1, 0, 0;...
   0, 0, -1, 0, 0, 0;...
   0, 0, 0, 0, 0,-1;...
   0, 0, 0, 0, -1, 0;] 
%%
kpoint=[0.23,0.45,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C*hk1*C'-hk2

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Mirror y M_y               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("NbSe2");
[nbands,~,nrpts]=size(g.ham);
theta_deg=180;
sigma_z=[1,0;0,1];
theta = theta_deg * pi / 180;
R=[1,0,0;...
    0,-1,0;...
    0,0,1]
kpoint1=[0,0,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[0,  1, 0, 0, 0, 0;...
   -1, 0, 0, 0, 0, 0;...
   0, 0, 0, 1, 0, 0;...
   0, 0, -1, 0, 0, 0;...
   0, 0, 0, 0, 0,1;...
   0, 0, 0, 0, -1, 0;] 
%
kpoint=[0.23,0.45,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C*hk1*C'-hk2

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Mirror y M_xy               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=[-0.5,    -sin(2*pi/3),0.0;...
    -sin(2*pi/3), 0.5,   0.0;...
    0,0,1]
% kpoint1=[0.2,0.25,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[0, -1, 0, 0, 0, 0;...
   1, 0, 0, 0, 0, 0;...
   0, 0, 0,exp(1j*1*pi/3), 0, 0;...
   0, 0, exp(1j*2*pi/3), 0, 0, 0;...
   0, 0, 0, 0, 0,exp(1j*5*pi/3);...
   0, 0, 0, 0, exp(1j*4*pi/3), 0;];
% C=[0,-1, 0, 0, 0, 0;...
%    1, 0, 0, 0, 0, 0;...
%    0, 0, 0,-1, 0, 0;...
%    0, 0, 1, 0, 0, 0;...
%    0, 0, 0, 0, 0,-1;...
%    0, 0, 0, 0, 1, 0;] 
%
kpoint=[0.5,0.2,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
% kpoint2=[2/3,1/3,0.0];
% kpoint2=-[1/3,1/3,0]
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C'*hk1*C-hk2
% cn=eig(Psik(:,1:2)'*C*Psik(:,1:2))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Mirror y M_x               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=[-0.5,    sin(2*pi/3),0.0;...
    sin(2*pi/3), 0.5,   0.0;...
    0,0,1]
kpoint1=[0.2,0.25,0.0];
% C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
C=[0, 1, 0, 0, 0, 0;...
   -1, 0, 0, 0, 0, 0;...
   0, 0, 0,exp(1j*2*pi/3), 0, 0;...
   0, 0, exp(1j*1*pi/3), 0, 0, 0;...
   0, 0, 0, 0, 0,exp(1j*4*pi/3);...
   0, 0, 0, 0, exp(1j*5*pi/3), 0;]
%
kpoint=[0.5,0.0,0.0];
[~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
% [~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint*g.b)')'*inv(g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
h=C'*hk1*C-hk2

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Check Mirror M_y  for slab               %%%%
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
del_index=find(gs.wpos(:,2)>3);
%%
kpoint1=[0.1,0.0];
kpoint1b=kpoint1*gs.b2;
hk1=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1,gs.a2);
tic;
[Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q_odd(g.ham,g.hopr2,del_index,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
toc;

% hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2);
[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=kpoint1;
kpoint2b=kpoint2*gs.b2;
hk2=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2,gs.a2);
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
C=kron(P,C);
h=C'*hk1*C-hk2
max(h,[],'all')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Check Mirror M_y  for slab BdG           %%%%
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
[~,~,hk1]=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint1b,gs.a2,gs.b2,mu,delta)

[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=kpoint1;
kpoint2b=kpoint2*gs.b2;
% hk2=MTB.ham.get_slab_hk_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2,gs.a2);
[~,Psik2,hk2]=MTB.ham.get_slab_hk_BdG_at_q_odd_v2(gs.ham,gs.hopr2,del_index,nslab,nbands,nrpts,kpoint2b,gs.a2,gs.b2,mu,delta)
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

% max(hh,[],'all')
% m=eig(Psik1(:,1:size(P)/2)'*P*Psik1(:,1:size(P)/2));
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("NbSe2");
%%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=51;
%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs=read_fplo("NbSe2");
MillerIndices=[1,0,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

[nbands,~,nrpts]=size(gs.ham);
labels={'\Gamma','X','\Gamma'};
hkpoints={[0.0,0.0],...
          [0.5,0.0],...
          [0.0,0.0]};% hkpoints-high symmetry k points

Np=3;
omegamin=1;
omegamax=2.3;
omeganum=500;
omegas=linspace(omegamin,omegamax,omeganum);
nk=201;
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
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("Ta2PdSe6");

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
nslab=300;
delta=0.03;% 0.03~0.05
% del_index=find(gs.wpos(:,2)<0.1);
del_index=find(gs.wpos(:,2)>3);
%%
% del_index=[]
numEigs=50;
mus=linspace(0,0.5,401);%mu for E_f-mu to E_f+mu of 100 points
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
g=read_fplo("Ta2PdSe6");

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
knum=1001;
% % kxline=linspace(-0.5,0.5,knum);
% % kyline=linspace(-0.5,0.5,knum);
% % theta=pi/6;
% % R=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
% % kxline=linspace(-0.5,0.5,knum)*sqrt(3);
% % kyline=linspace(0,0,knum);
% kxline=linspace(0,0,knum);
% kyline=linspace(-0.5,0.5,knum);
% 
kxline=linspace(-0.5,0.5,knum);
kyline=linspace(0,0,knum);

kline=[kxline',kyline',zeros(knum,1)];

kline=kline*gs.b;

% kline=(R*kline')'
% kline=kline*g.b;

mus=linspace(-0.4,0.4,501);
Enk=MTB.ham.get_kline_bands(gs,kline);
% Enk=MTB.ham.get_kline_bands(g,kline);
eta=10^-2;
%
A = compute_spectral_function_kline_mu(Enk, mus, eta);

%plot by pcolor
 figure()
 imagesc(A); axis equal; axis off; colormap hot; colorbar;
 
%% plot by imagesc
 figure()
 % imagesc(A'); axis equal; axis off; colormap hot; colorbar;
 [Kx,Ky]=meshgrid(mus,kxline)
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
        g = MTB.read_poscar(g,"data/Ta2PdSe6/fplo/POSCAR");
        pos=textread("data/Ta2PdSe6/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/Ta2PdSe6/fplo/mydata-p1");
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


