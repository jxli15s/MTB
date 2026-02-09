%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       Read PtP2 Hamiltonian                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%p=parpool(8)
g = MTB.geometry("PtP2");
g = MTB.read_poscar(g,"data/PtP2/dp_orbital/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/PtP2/dp_orbital/wannier90_hr_p1.dat','data/PtP2/dp_orbital/wannier90_hr_p2.dat');
g.wpos=g.atoms;
%%
%% ===Calculate bulk bands======
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points
efermi=-3.5747;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kindex,"FeSe",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(30,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(32,:)-efermi,"Color",'red','LineWidth',2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs = MTB.geometry("FeSC");
gs = MTB.read_poscar(gs,"data/PtP2/sdsp_orbital/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/PtP2/sdsp_orbital/wannier90_hr_p1.dat','data/PtP2/sdsp_orbital/wannier90_hr_p2.dat');
gs.wpos=gs.atoms;%*gs.a;

MillerIndices=[1,0,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=50;
%%
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.0],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=51;
efermi=-3.5747;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("PtP2");
gs = MTB.read_poscar(gs,"data/PtP2/sdsp_orbital/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/PtP2/sdsp_orbital/wannier90_hr_p1.dat','data/PtP2/sdsp_orbital/wannier90_hr_p2.dat');
gs.wpos=gs.atoms%*gs.a;
MillerIndices=[1,0,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

[nbands,~,nrpts]=size(gs.ham);
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
labels={'X','\Gamma','X'};
Np=1;
efermi=-3.5747;
omegamin=efermi-1;
omegamax=efermi+1;
omeganum=300;
omegas=linspace(omegamin,omegamax,omeganum);
nk=201;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);
%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('magma')); %magma plasma inferno cividis inferno hot heat

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
caxis([1, 50])
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Calculate Corner   hexagonal               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("PtP2");
gs = MTB.read_poscar(gs,"data/PtP2/dp_orbital/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/PtP2/dp_orbital/wannier90_hr_p1.dat','data/PtP2/dp_orbital/wannier90_hr_p2.dat');
efermi=-3.5747;
% gs.atoms=kron(gs.atoms,ones(3,1));
% gs.atoms(:,3)=0;
gs.wpos=gs.atoms;%*gs.a;
%===Construct xy opened model======
%
n1=21;
n2=21;
tic;
gs = MTB.ham.get_xy_open_wannier(gs,n1,n2);
toc;
%
gs.ham=gs.ham{1};


%
% 1) 计算并平移中心到 (0,0,0)
center = mean(gs.wpos, 1);       % 1×3 向量 [mean_x, mean_y, mean_z]
gs.wpos = gs.wpos - [0.5,0.5,0.0]*gs.a;      % 所有坐标减去均值

figure()
hold on;
plot(gs.wpos(:,1),gs.wpos(:,2),'ro')

% plot(gs.wpos(:,1),gs.wpos(:,2),'ro')

% 2) 定义一个规则六边形（边长 a，可根据需求调整）
% 1)原始六边形顶点
% a = 69.3;
% theta = (0:5)' * 2*pi/6;
% hexVx = a * cos(theta);
% hexVy = a * sin(theta);
% 矩形
a = 12*sqrt(2)*2+1.8;
theta1 = (0:1)' * 2*pi/2;
theta2 = theta1 + pi/2;
theta= [theta1(:);theta2(:)];
theta=sort(theta);
hexVx = a * cos(theta);
hexVy = a * sin(theta);
% 三角形
% theta = (0:2)' * 2*pi/3;
% hexVx = a * cos(theta);
% hexVy = a * sin(theta);
% --- 在这里做 30° 旋转 ---
phi = 45;   % 角度制
R = [ cosd(phi), -sind(phi);
      sind(phi),  cosd(phi) ];

% 将所有顶点组成 2×6 矩阵，左乘 R
V = R * [hexVx'; hexVy'];  
hexVx = V(1, :)';   % 旋转后的 x
hexVy = V(2, :)';   % 旋转后的 y

% 3) 找出落在该六边形内的点 idxHex
xy = gs.wpos(:,1:2);             % 取 xy 平面坐标
inHex = inpolygon(xy(:,1), xy(:,2), hexVx, hexVy);
idxHex = find(inHex);

% %% 4) 可选：按 z 坐标近似等于 z0 进行筛选
% z0  = 9.3108;
% tol = 1e-6;
% idxZ = find( abs(gs.wpos(:,3) - z0) < tol );

% 5) 组合筛选条件（这里取交集，你也可以只用 idxHex 或 idxZ）
% idx = intersect(idxHex, idxZ);
idx=idxHex;

% 6) 截取新的 wpos 和 ham
gs.wpos = gs.wpos(idx, :);
gs.ham = gs.ham(idx, idx);
%
plot(gs.wpos(:,1),gs.wpos(:,2),'bo',LineWidth=3)


plot(hexVx,hexVy,'k-')
axis('equal')
%%
% [V, D] = eig(full(1/2*gs.ham{1}+1/2*gs.ham{1}'));
hk1=gs.ham;
[V, D] = eig(full(1/2*hk1+1/2*hk1'));
[E,ind]=sort(diag(D));
Psik=V(:,ind);
%%
figure()
plot(E-efermi,'ro')
%%
figure()
x=gs.wpos(1:end,1);
y=gs.wpos(1:end,2);
p=abs(Psik).^2;
scatter3(x, y, sum(p(:, 2407:2410), 2), 100, sum(p(:,2407:2410), 2), 'filled');
% scatter3(x, y, sum(p(:, 582:584), 2), 100, sum(p(:,582:584), 2), 'filled');
% view(2);
%%
for i=2200:2500
    fig = figure('Visible', 'off');
    % p(1:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i);
    % scatter(x(:), y(:), 100, abs(p(:,i)).^2, 'filled');
    scatter3(x(:), y(:),p(:,i), 100, p(:,i), 'filled');
    caxis([0, max(p(:,i))]);
    view(2)
    wxname="data/PtP2/dp_orbital/states/"+int2str(i)+'.png';
    saveas(fig, wxname)
end