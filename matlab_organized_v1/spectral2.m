%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
% g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);

% Step 4: Construct Supercell
n1 = 1; n2 = 1; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

n1 = 15; n2 = 1; % Supercell dimensions
gs = constructSupercell(g, n1, n2);
g.wpos(:,3)=0; gs.wpos(:,3)=0;
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy_ori,kpath_ori,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy_ori,nbands,efermi,kpath_ori,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
Vamp=0.0
gs.ham=gs.iniham+0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[7.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [7.5,0.0,0.0]};% hkpoints-high symmetry k points
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.3,0.3])
%%
E_grid = linspace(-0.5,0.5, 100);  % 能量网格
A_kE = zeros(length(kpath), length(E_grid));  % 初始化谱函数
sigma=0.05;
[~, idx] = min(abs(E_grid)); % 找到最小绝对值的索引
v1=linspace(0.0,0.0,idx+16);
v2=linspace(0.0,1,length(E_grid)-idx-15);
Vs=[v1,v2(2:end)];
%%
for iv=1:length(E_grid)
    Vamp=Vs(iv);
    gs.ham=gs.iniham+0;
    gs=moire_potential(g,gs,Vamp);
    Electric_field_in_evpA=0.0;
    [nbands,~,nrpts]=size(gs.ham);
    labels={'X','\Gamma','X'}; % labels for k
    hkpoints={[7.5,0.0,0.0],...
        [0.0,0.0,0.0],...
        [7.5,0.0,0.0]};% hkpoints-high symmetry k points
    % labels={'X','\Gamma','Y'}; % labels for k
    % hkpoints={[7.5,0.0,0.0],...
    %     [0.0,0.0,0.0],...
    %     [0.0,0.5,0.0]};% hkpoints-high symmetry k points
    tic;
    [Wk,psik,uksuper]=MTB.ham.get_bulk_unfolding_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,g,gs);
    [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
    toc;
    fprintf("Energy Loop %d  done! \n",iv)
    %%
    % Wk=abs(Wk)
    % if iv<idx+18
    %     Wk(61,:)=Wk(60,:)/8;
    %     Wk(63,:)=Wk(58,:)/8;
    %     Wk(65,:)=Wk(56,:)/8;
    %     Wk(67,:)=Wk(54,:)/8;
    % end
    for i=1:length(kpath)
        for j=1:nbands
        A_kE(i, iv) = A_kE(i, iv) + ...
            Wk(j,i) * exp(-((E_grid(iv) - Energy(j,i)) .^ 2) / sigma^2) / (sqrt(pi) * sigma);
        end
    end
end
%% =========== 2. 绘制 Band Unfolding ===========
% save('Spectral_weight.mat','kpath','Energy','Energy_ori',"E_grid","Wk",'A_kE')
load('Spectral_weight.mat')
figure;
% A_kE=abs(A_kE);
imagesc(kpath, E_grid-0.022, abs(A_kE'));
set(gca, 'YDir', 'normal'); % 反转 y 轴，使低能级在下方
% colormap('inferno'); % 使用和 Python 类似的 colormap
% colormap('hot')
% colormap(slanCM('RdBu'))
colormap(slanCM('inferno'))
% colormap(flipud(colormap));
shading interp

hold on;

% 画超胞能带
% for i = 1:nbands
%     plot(kpath,Energy(i,:)-0.022, 'c', 'LineWidth', 0.5, 'LineStyle', '--', 'Color', [0, 1, 1, 0.9]);
% end

for i = 1:size(g.ham,1)
    plot(kpath_ori,Energy_ori(i,:)-0.022, 'c', 'LineWidth', 1, 'LineStyle', '--', 'Color', [0.2, 1, 1, 0.9]);
end

for i=1:length(kindex)-2
     plot([kindex(i+1) kindex(i+1)],[-0.3 0.301],'--w','LineWidth',1)
end

plot(kpath,zeros(1,length(kpath)),'--w','LineWidth',1)
% xlim(kpath(1),kpath(length(kpath)/2),kpath(end))
xticks([0,kpath(101),kpath(end)])
xticklabels(labels)
xlim([kpath(101)/2,kpath(101)+kpath(101)/2])
ylim([-0.3,0.301])
% xlabel('$k$ (BZ)', 'Interpreter', 'latex');
ylabel('Energy-E_f (eV)');
% title('1D Band Unfolding for Double-Orbital Unit Cell');
colorbar;
caxis([0,300])

set(gca, 'FontSize', 14);
print('Spectral_Weight_off.png','-dpng','-r600')


%%
n1=3;
n2=3;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
% labels={'R','Y','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
labels={'M','\Gamma','K'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;

stepmax=10000;
minstepmax=11;
step=1;
knum=31;
u=(1*n1*n2)/nbands;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
U=20;

xinitial=zeros(nbands,nbands);
nsite=nbands/2;
for i = 1:nsite
    % nup=randn(1); 
    % ndn=1-nup;
    nup=0.5;
    ndn=0.5;
    sx=randn(1);
    sy=1-sx;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U;
end
% save("tmp.dat","tmp")
 % load("tmp.mat")
 % xinitial=tmp;

[xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);


gs.onsite_modify(xinitial);
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(1*n1*n2,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(1*n1*n2+1,:),"Color",'red','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);


% knum=81;
% band1=1;
% band2=1*n1*n2;
% [wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
spin_plot(gs,si)
%%
lattice_plot(gs)
%%
T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
%%
g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)]
% g.wpos=g.wpos
%%
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.2,0.2,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.2,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(4),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")



%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

gs = MTB.geometry("TaIrTe4_s");
%%
T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

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
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)]
% g.wpos=g.wpos
%%
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);

Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);

[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;
bandindex=1:120;


[Energy,Omega_k,kpath,kindex]=MTB.ham.get_bulk_bands_bcd(gs,hkpoints,nk,bandindex);

%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
% hold on;
% colormap(slanCM('RdBu'))
% scatter(kpath,Energy(62,:),[],(Omega_k(:,61)+Omega_k(:,62)),"filled");
% scatter(kpath,Energy(63,:),[],(Omega_k(:,63)+Omega_k(:,64)),"filled");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.2,0.2,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.2,-0.2,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
s2=[0  -1i
    1i  0];
T=kron(eye(60),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")

%%
%%
%%
n1=15;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;
%%
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(4*n1+1,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+2,:),"Color",'red','LineWidth',2);
%%

%%
clc;
clear;

g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');


T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)];
% g.wpos=g.wpos

n1=15;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
% labels={'R','Y','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;

stepmax=400;
minstepmax=11;
step=1;
knum=21;
u=(4*n1+2)/nbands;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
U=1;

xinitial=zeros(nbands,nbands);
nsite=nbands/2;
for i = 1:nsite
    nup=3; 
    ndn=3;
    sx=0.0;
    sy=0.0;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy)/8;-(sx+1j*sy)/8,nup].*U;
end
%save("tmp.dat","tmp")
 load("tmp.mat")
 xinitial=tmp;

[xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);

%%
gs.onsite_modify(xinitial*10);
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(4*n1+1,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+2,:),"Color",'red','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);


knum=81;
band1=1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
check_hermitian_elements(g)
%%
function check_hermitian_elements(obj)
    % H: Hamiltonian 矩阵，大小为 (nbands, nbands, nhopr)
    % g.hopr: hopping 位置矩阵，大小为 (nhopr, 3)

    [nbands, ~, nhopr] = size(obj.ham);
    
    % 遍历所有 (i, j)
    for i = 1:nbands
        for j = 1:nbands
            % 遍历所有 hopping 位置 R
            for k = 1:nhopr
                R = obj.hopr(k, :);  % 当前 hopping 位置
                Hij = obj.ham(i, j, k);  % H_ij(R)
                
                % 找到对应的 -R 位置
                idx_negR = find(all(obj.hopr == -R, 2));
                
                if isempty(idx_negR)
                    warning('No matching -R found for R = (%d, %d, %d)', R(1), R(2), R(3));
                    continue;
                end
                
                Hji_conj = conj(obj.ham(j, i, idx_negR)); % H^*_ji(-R)

                % 检查是否满足厄米性
                if abs(Hij - Hji_conj) > 1e-10
                    fprintf('Non-Hermitian element found at (i=%d, j=%d) for R = (%d, %d, %d)\n', ...
                            i, j, R(1), R(2), R(3));
                    fprintf('H(%d,%d, [%d %d %d]) = %f + %fi\n', i, j, R(1), R(2), R(3), real(Hij), imag(Hij));
                    fprintf('H*(%d,%d, [%d %d %d]) = %f + %fi\n', j, i, -R(1), -R(2), -R(3), real(Hji_conj), -imag(Hji_conj));
                    fprintf('Difference: %f + %fi\n\n', real(Hij - Hji_conj), imag(Hij - Hji_conj));
                    break;
                end
            end
        end
    end
    
    disp('Hamiltonian elements checked.');
end


%%
function T = getBasisTransformMatrix()
    % getBasisTransformMatrix Returns the basis transformation matrix
    %
    % Outputs:
    %   T - Matrix, the basis transformation matrix for Hamiltonian
    %
    % The transformation matrix `T` is used to change the basis of the
    % Hamiltonian matrices.
    % % 
    % % T = [1,0,0,0,0,0,0,0; ...
    % %      0,0,0,0,0,1,0,0; ...
    % %      0,0,1,0,0,0,0,0; ...
    % %      0,0,0,0,0,0,0,1; ...
    % %      0,1,0,0,0,0,0,0; ...
    % %      0,0,0,0,1,0,0,0; ...
    % %      0,0,0,1,0,0,0,0; ...
    % %      0,0,0,0,0,0,1,0];

    T = [1,0,0,0,0,0,0,0; ...
         0,0,0,0,1,0,0,0; ...
         0,1,0,0,0,0,0,0; ...
         0,0,0,0,0,1,0,0; ...
         0,0,1,0,0,0,0,0; ...
         0,0,0,0,0,0,1,0; ...
         0,0,0,1,0,0,0,0; ...
         0,0,0,0,0,0,0,1];
end

function ham = transformBasis(ham, T)
    % transformBasis Applies a basis transformation to the Hamiltonian
    %
    % Inputs:
    %   ham - 3D array, the original Hamiltonian in the initial basis
    %   T   - Matrix, the basis transformation matrix
    %
    % Outputs:
    %   ham - 3D array, the Hamiltonian in the transformed basis

    % Apply the transformation for each Hamiltonian slice
    for i = 1:size(ham, 3)
        ham(:, :, i) = T * ham(:, :, i) * inv(T);
    end
end

function wpos = setWannierPosition(g)
    % setWannierPosition Sets the Wannier positions for the geometry
    %
    % Inputs:
    %   g    - Struct, contains the geometry and atomic positions
    %
    % Outputs:
    %   wpos - Array, the Wannier positions for each orbital
    %
    % This function computes the Wannier positions based on the atomic
    % positions and lattice vectors.

    % Compute Wannier positions based on atomic positions
    wpos = g.atoms * g.a;
    wpos = [wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :); ...
            wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :)];
end

function gs = constructSupercell(g, n1, n2)
    % constructSupercell Constructs the supercell Hamiltonian
    %
    % Inputs:
    %   g  - Struct, contains the geometry and Hamiltonian of the material
    %   n1 - Integer, number of cells along the first lattice vector
    %   n2 - Integer, number of cells along the second lattice vector
    %
    % Outputs:
    %   gs - Struct, the supercell Hamiltonian and related properties

    % Generate supercell Hamiltonian
    gs = MTB.ham.get_supercell_wannier(g, n1, n2);
    % Store the initial Hamiltonian for reference
    gs.iniham = gs.ham+0;
end

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

function [xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
        end


        tap=sum(abs(xnew-xinitial),"all");

    
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end
        xinitial=xnew*0.8+0.2*xinitial;
        step=step+1;
    end
end

function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
    sitenum=size(Unk,2)/2;
    paulix=[0,1;1,0];pauliy=[0,-1i;1i,0];pauliz=[1,0;0,-1];
    parfor i=1:size(kindex,1)
        nki=abs(Unk(:,bandindex(i),kindex(i))).^2+nki;
        psik=reshape(Unk(:,bandindex(i),kindex(i)),[2,sitenum]);
        sxki=psik'*paulix*psik./2+sxki;
        syki=psik'*pauliy*psik./2+syki;
        szki=psik'*pauliz*psik./2+szki;
    end
    knum=size(Unk,3);
    ni=nki/knum;
    ni=reshape(ni,[2,size(ni,1)/2]);%first row spin up, second row spin down
    si=real([diag(sxki).';diag(syki).';diag(szki).']./knum);
end

function [Etot,kindex,bandindex]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
end

%%

function obj=add_elec(obj,Electric_field_in_evpA)
    % obj.wpos(:,3)=round(obj.wpos(:,3));
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        % obj.wpos(i,3)*Electric_field_in_evpA
    end 
end



function []=spin_plot(obj,si)
figure('Color','white')
pos_cart=obj.atoms;
pos_cart=pos_cart*obj.a;
h=zeros(1,3);
% plot(obj.wpos(:,1),obj.wpos(:,2),'o')
lattice_plot(obj)
hold on
si=si.';
si=si./sqrt(2);
h(2)=quiver(pos_cart(:,1),pos_cart(:,2),real(si(:,1)),real(si(:,2)),0,'Color','#007EC9','LineWidth',2,'DisplayName','xy');
hold on
% h(3)=quiver(pos_cart(:,1),pos_cart(:,2),zeros(size(si(:,1))),real(si(:,3)),0,'Color','#FF7E81','LineWidth',2,'DisplayName','z');
% legend(h(2:3))
set(gca,'Fontsize',18,'FontName','Times New Roman','linewidth',0.8)

end

function lattice_plot(obj)
lattice=obj.a;
sublattice=obj.atoms;
figure('Color','white')
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} ) ;
drawArrow([0,lattice(1,1)],[0,lattice(1,2)],'linewidth',2,'color','k')
hold on
drawArrow([0,lattice(2,1)],[0,lattice(2,2)],'linewidth',2,'color','k')
hold on
sublattice(:,1:2)=sublattice(:,1:2)*(lattice(1:2,1:2));


for i=1:size(sublattice,1)
    if sublattice(i,end)==1

        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#4DBEEE','MarkerSize',40,'MarkerFaceColor','k')
    elseif sublattice(i,end)==2
        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#4DBEEE','MarkerSize',40,'MarkerFaceColor','k')
    else
        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#0072BD','MarkerSize',40,'MarkerFaceColor','k')
    end
    text(sublattice(i,1)+0.05,sublattice(i,2)+0.05,num2str(i))
end
axis equal
axis off
end


function u_k_super=construct_supercell_wavefunction(kpoints,Psik, g,gs)
for ik=1:length(kpoints)
    S=gs.a*inv(g.a);
    S2u=[]
    num_wpos_u=size(g.wpos,1)
    num_wpos_s=size(gs.wpos,1)
    for i =1:num_wpos_s
        z=mod(i,num_wpos_u);
        if z==0
            z=num_wpos_u
        end
        S2u=[S2u,z]
    end
end
end

function g = initializeGeometry(materialName, dataPath)
    % initializeGeometry Initializes the material geometry and Hamiltonian
    %
    % Inputs:
    %   materialName - String, the name of the material (e.g., "TaIrTe4")
    %   dataPath     - String, the path to the directory containing input files
    %
    % Outputs:
    %   g            - Struct, contains the geometry and Hamiltonian information
    %
    % This function reads the material geometry from a POSCAR file and
    % constructs the Hamiltonian using Wannier90 HR files.

    % Create geometry object
    g = MTB.geometry(materialName);
    % Read geometry from POSCAR file
    g = MTB.read_poscar(g, fullfile(dataPath, "POSCAR"));
    % Read Wannier Hamiltonian data
    [g.ham, g.hopr] = MTB.wannier.read_hr(...
        fullfile(dataPath, "wannier90_hr_p1.dat"), ...
        fullfile(dataPath, "wannier90_hr_p2.dat"));
end