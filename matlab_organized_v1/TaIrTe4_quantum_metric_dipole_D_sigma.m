%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% load('xinitial_good.mat')
% Step 1: Initialize Geometry and Hamiltonian
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);

% Step 4: Construct Supercell
n1 = 15; n2 = 1; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

gs.iniham = gs.ham + 0; % Store the initial Ham
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

% Step 5: Neighbor Search
[result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs);
%
% Step 6: Initialize fixed Hartree-Fock States parameters
nec1 = (4 * n1 + 4) / n1 / 2 / 2 / 2;
nec2 = (4 * n1) / n1 / 2 / 2 / 2;
% xinitial_0 = diag(kron(ones(1, n1), [nec1 - 0.4, nec1 - 0.4, nec2 + 0.4, nec2 + 0.4, ...
%                                     nec1 + 0.4, nec1 + 0.4, nec2 - 0.4, nec2 - 0.4]));
xinitial_0 = diag(kron(ones(1, n1*n2), rand(1,8)));

[nbands,~,nrpts]=size(gs.ham);
% labels={'Y','\Gamma','X','R'}; % labels for k
% hkpoints={[0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
%
labels={'R','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
% labels={'Y','\Gamma','Y','R','X','\Gamma','-X'}; % labels for k
% hkpoints={[0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,-0.5,0.0],...
%           [0.5,-0.5,0.0],...
%           [0.5,-0.0,0.0],...
%           [0,0,0],...
%           [-0.5,0.0,0.0]};% hkpoints-high symmetry k points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Read and Write the wannier90_hr.dat for TI    %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=load("data/tit_hf/15x1/10nmd/U-epsilon-f/"+int2str(12)+"-15x1"+"/result_s1_eps"+int2str(13)+".00_r"+int2str(11)+".mat");
U0=result.U0;
U=result.V;
V=result.V;
xinitial=result.xinitial;
pairsU=result.pairsU;
pairsV=result.pairsV;
efermi=result.efermi;
modifyHam(gs, xinitial, V, V, pairsU, pairsV)
% Add Zeeman and Electric Field
%
Electric_field_in_evpA=0.1;
gs=add_elec(gs,Electric_field_in_evpA);

s3=[1  0
    0  -1];
Zeeman=kron(eye(60),0.001*s3);
gs.add_zeeman(Zeeman)
% plot the band along HSL
%%
Electric_field_in_evpA=0.00*0.529177; nk=251;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1+1,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
ylim([-0.1,0.1])
%%
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
%%
[nbands,~,nrpts]=size(gs.ham);
hkpoints={[-0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
labels={'-X','\Gamma','X'};
Np=2;
omegamin=efermi-0.1;
omegamax=efermi+0.1;
omeganum=50;
omegas=linspace(omegamin,omegamax,omeganum);
nk=51;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);
%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('inferno')); %magma plasma inferno cividis inferno hot heat
% caxis([0, 32])
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
% caxis([0, 32])
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
% save("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K_band.mat",'Energy','kpath','kindex','nbands','efermi','Dk_core','gk','TDos_new','efermi')
%%
% Calculate the Ham on K-mesh plane
knum  = 200;
kxline = [-0.5,0.5];
kyline = [-0.5,0.5];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]  = MTB.ham.get_bulk_plane_bands(g, Kx,Ky,Kz);
[Hamk,Unk,Enk]  = MTB.ham.get_bulk_plane_bands_with_Ham(gs, Kx,Ky,Kz);
efermi=calculate_ef(Enk(:),0.5);
Enk=Enk-efermi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the DOS and TDOS for the carrier density  %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plottap=2;
nk=knum;
Enum=3000;
Emin=-0.1;
Emax=0.1;
eps=(Emax-Emin)/Enum*10;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
Dos=Dos*(Emax-Emin)/Enum;
TDos_new=TDos/norm(cross(gs.a(1,:),gs.a(2,:)))*10^16;

% plot the Dos and TDos
% figure('Color','white')
% hold on;
% plot(Eaxis,Dos*10^15,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
% xlabel('E(eV)')
% ylabel('n (cm^{-2})')
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
% legend('Density')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       calculate the QM and QMD on this plane          %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. dk_vecs: 对应两个方向的 k 步长, 比如沿 b1, b2
% ===== 定义 dk_list (2D 情况) =====
b1 = gs.b(1,:);    % 你的 reciprocal vectors，如果是别的字段，就替换
b2 = gs.b(2,:);
dk_vecs = [b1/knum;   % 对应 "kx 方向"的单位向量
           b2/knum];  % 对应 "ky 方向"的单位向量
% 3. 其它参数
% band_list =55:66;               % 比如 valence top + conduction bottom…
band_list = 59:62;
% band_list =1:8;  
Emin=-0.1;
Emax=0.1;
NEF=3000;
Ef_list   = linspace(Emin,Emax,NEF);
AreaBZ = norm(gs.b(1,:))*norm(gs.b(2,:));
eta       = 1e-4;
weights   = [];               % 让函数内部设置均匀 AreaBZ/(Nkx*Nky)
deltaE_reg = 1e-5;

TK=30;
sigma_abc = MTB.ham.get_sigma_quantum_metric_dipole( ...
    Hamk, Unk, Enk, band_list, dk_vecs, Ef_list, ...
    AreaBZ, eta, weights, deltaE_reg,TK);

[Dk_core,gk] = MTB.ham.get_Dk_qmd_plainD_core(Hamk, Unk, Enk, band_list, dk_vecs, deltaE_reg);

%
file_name="QMD_B_1meV_30K.mat";
save(file_name,"sigma_abc","Dk_core","gk","gs","Kx","Ky","Dos","TDos","TDos_new","TK","Ef_list","band_list","Enk","efermi","-v7.3");
%%
siamg_abc_all=zeros(2,2,2,3000,10);
for i =1:9
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_0"+int2str(i)+"meV_30K.mat";
load(filename)
sigma_abc_all(:,:,:,:,i)=sigma_abc;
end
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat";
load(filename)
sigma_abc_all(:,:,:,:,10)=sigma_abc;
%%
filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_sigma_TDos_all_30K.mat";
load(filename)
Emin=-0.1;
Emax=0.1;
NEF=3000;
Ef_list   = linspace(Emin,Emax,NEF);
%%
figure()
abc=[1,2,2];
a=abc(1);b=abc(2);c=abc(3);
peak_1=[]
% figure()
hold on
N = 201;
cmap = [linspace(0,1,N)', zeros(N,1), linspace(1,0,N)'];  % [R,G,B]：红->蓝
for i=100:201
    y=1000*squeeze(sigma_abc_all(a,b,c,:,i))*pi;
    peak_1=[peak_1,y(1755)];
    % plot(Ef_list,y,'LineWidth',1.5,'Color',cmap(i,:))
    dos=squeeze(TDos_new_all(1,:,i));
    gapMask=(abs(dos)<7e7);
    dos(gapMask)=NaN;
    plot(dos,y,'LineWidth',1.5,'Color',cmap(i,:))
    plot(dos,-y,'LineWidth',1.5,'Color',cmap(201-i+1,:))
end

T_list=1:N
legend(arrayfun(@(T)sprintf('%g meV',T/10), T_list, 'UniformOutput', false), ...
       'Location','best');

box on;
%%
figure(); hold on
N = 10;
cm = slanCM('RdBu');
idx = round(linspace(1,size(cm,1),N));
cmap = cm(idx,:);
cmap = flipud(cmap); % 如需反向

colororder(gca, cmap);   % R2019b+，老版本可用 set(gca,'ColorOrder',cmap)

for i = 1:10
    plot(Ef_list, 1000*squeeze(sigma_abc_all(a,b,c,:,i)), 'LineWidth', 1.5);
end
% T_list=1:N
% legend(arrayfun(@(T)sprintf('%g meV',T/10), T_list, 'UniformOutput', false), ...
%        'Location','best');

%%
figure()
plot(1:N,peak_1,'-o')
%%
z=1000*squeeze(-sigma_abc_all(2,1,1,:,:));
kx=squeeze(TDos_new_all);
Blist=-100:100;
ky=kron(ones(3000,1),Blist);

figure()
surf(kx,ky,z,'EdgeColor','none')
% pcolor(z')
colorbar
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
clim([-0.4,0.4])
xlim([-10e12,10e12])
shading interp
view(2)


%%
load("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat")
% filename="/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_09meV_30K.mat";
% sigma_abc=squeeze(sigma_abc_all(:,:,:,:,1))
% TDos_new=squeeze(TDos_new_all(1,:,1))
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            % plot(Ef_list,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            plot(TDos_new,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            % filename="sigma_abc_"+int2str(x)+int2str(y)+int2str(z)+".dat";
            % outlist=[Ef_list.',reshape(sigma_abc(x,y,z,:),300,1)];
            % writeoutput(filename,outlist)
        end
    end
end
% plot(Ef_list,1000*sigma_yxx,'ko')
ylabel('mA/V^2')
xlabel('Energy(eV)')

%%
theta = linspace(0,2*pi,721);   % 0..360 deg
E0 = 1;                         % 你也可以取 rms 或 1，反正只差整体系数

s_yxx = 0.0;    % sigma_yxx
s_xyy = -1;    % sigma_xyy

Jx = s_xyy * sin(theta).^2.*sin(theta);
Jy = s_yxx * cos(theta).^2.*cos(theta);

% Jpar = Jx.*cos(theta) + Jy.*sin(theta);   % longitudinal along E

figure; hold on
plot(theta*180/pi, Jx+Jy, 'LineWidth', 1.5)
% plot(theta*180/pi, Jy, 'LineWidth', 1.5)
% plot(theta*180/pi, Jpar,'LineWidth', 1.5)
yline(0,'k--','LineWidth',1.2)
xlabel('\theta (deg)'); xlim([0 360])
legend('J_x','J_y','J_{||}','Location','best')
box on;

figure()
polarplot(theta,Jx+Jy)
%%
abc=[2,1,1];
a=abc(1);b=abc(2);c=abc(3);
knum=500;
band1=zeros(knum/2,4);
banddk1=zeros(knum/2,4);
bandgxx1=zeros(knum/2,4);
bandgyx1=zeros(knum/2,4);

idx=0;
for i=knum:-1:knum/2+1
    idx=idx+1;
    for j=59:62
        band1(idx,j-58)=Enk(i,i,j);
        banddk1(idx,j-58)=Dk_core(i,i,a,b,c,j-54);
        bandgxx1(idx,j-58)=gk(i,i,1,1,j-54);
        bandgxx1(idx,j-58)=gk(i,i,2,1,j-54);
    end
end

band2=zeros(knum/2,4);
banddk2=zeros(knum/2,4);
bandgxx2=zeros(knum/2,4);
bandgyx2=zeros(knum/2,4);
idx=0
for i=knum/2+1:knum
    idx=idx+1;
    for j=59:62
        band2(idx,j-58)=Enk(knum/2+1,i,j);
        banddk2(idx,j-58)=Dk_core(knum/2+1,i,a,b,c,j-54);
        bandgxx2(idx,j-58)=gk(knum/2+1,i,1,1,j-54);
        bandgyx2(idx,j-58)=gk(knum/2+1,i,2,1,j-54);
    end
end

band3=zeros(knum/2,4);
banddk3=zeros(knum/2,4);
bandgxx3=zeros(knum/2,4);
bandgyx3=zeros(knum/2,4);
idx=0
for i=knum:-1:knum/2+1
    idx=idx+1
    for j=59:62
        band3(idx,j-58)=Enk(knum/2*3+1-i,i,j);
        banddk3(idx,j-58)=Dk_core(knum/2*3+1-i,i,a,b,c,j-54);
        bandgxx3(idx,j-58)=gk(knum/2*3+1-i,i,1,1,j-54);
        bandgyx3(idx,j-58)=gk(knum/2*3+1-i,i,2,1,j-54);
    end
end

band=[band1;band2;band3];
banddk=[banddk1;banddk2;banddk3];
bandgxx=[bandgxx1;bandgxx2;bandgxx3];
bandgyx=[bandgyx1;bandgyx2;bandgyx3];
%%
figure()
kdist=kpath(1:3*knum/2);
msize=10;
Dmax=max(abs(banddk),[],'all');
hold on;
for i=1:4
% scatter(kdist,band(:,i),msize,banddk(:,i)/Dmax,'filled')
scatter(kdist,band(:,i),msize,bandgxx(:,i),'filled')
% plot(squeeze(Enk(100,:,i+54)))
end

kk=kindex;
linesize=1;
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',2)
for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[-0.1 0.1],'--k','LineWidth',linesize)
end

grid off
box on
colorbar
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp
% clim([-0.5,0.5])

% colorbar
% colormap(sky)
% clim([0,5000])
%%
 Ef=-0.0;
 band_list=55:66;
 Enk_sel = Enk(:,:,band_list);
 Nb_sel=length(band_list);
 eta=1e-3;
 delta_n = (1/pi)*eta ./ ((Enk_sel - Ef).^2 + eta^2); % Nkx x Nky x Nb_sel
 [Nkx,Nky]=size(Kx);
Dmap = zeros(Nkx,Nky);
for in=1:Nb_sel
    Dmap = Dmap + squeeze(Dk_core(:,:,1,2,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
end

% Dmap=squeeze(Dk_core(:,:,1,2,1,6)) .* delta_n(:,:,6) / (2*pi)^2/10^19;

D_yxx = sum(Dmap,'all')

figure()
surf(Kx,Ky,Dmap,'EdgeColor','none')
colorbar
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp
% clim([-1,1])
% zlim([-1,1])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                              初始几何和Ham                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


function [result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs)
    % findNeighbors Finds neighboring atoms for interactions
    %
    % Inputs:
    %   gs - Struct, contains the supercell Hamiltonian and Wannier positions
    %
    % Outputs:
    %   result_matrices - Cell array, contains neighbor data for various distances
    %   pairsU          - Cell array, contains pairs for on-site and off-site interactions (U)
    %   pairsV          - Cell array, contains pairs for on-site and off-site interactions (V)

    result_matrices = cell(1, 7);
    % Loop through different neighbor distances
    for i = 1:7
        result_matrices{i} = find_neighbor_data(gs.wpos(1:2:end, :), gs.a, i, 2);
    end

    % Filter and group pairs
    pairs_onsite = result_matrices{1}(result_matrices{1}(:, 1) == result_matrices{1}(:, 2), :);
    pairs_onsite_nn = result_matrices{1}(result_matrices{1}(:, 1) ~= result_matrices{1}(:, 2), :);

    pairs_onsite_1 = pairs_onsite(mod(pairs_onsite(:,1),4)==1,:);
    pairs_onsite_2 = pairs_onsite(mod(pairs_onsite(:,1),4)==2,:);
    pairs_onsite_3 = pairs_onsite(mod(pairs_onsite(:,1),4)==3,:);
    pairs_onsite_4 = pairs_onsite(mod(pairs_onsite(:,1),4)==0,:);
    pairs_onsite_nn_13 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==1,:);
    pairs_onsite_nn_24 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==0,:);


    pairs_offsite_nn = result_matrices{2};
    pairs_offsite_nn_13=result_matrices{2}(mod(result_matrices{2}(:, 1),2)~=0 & mod(result_matrices{2}(:, 2),2)~=0, :);
    pairs_offsite_nn_24=result_matrices{2}(mod(result_matrices{2}(:, 1),2)==0 & mod(result_matrices{2}(:, 2),2)==0, :);


    pairs_offsite_nnn = result_matrices{3};
    pairs_offsite_nnn_12=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==2,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_14=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_34=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==3,:)];
    pairs_offsite_nnn_32=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==3,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==2,:)];

    pairs_offsite_nnnn = result_matrices{4};
    pairs_offsite_nnnn_12=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==2,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_14=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_34=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==3,:)];
    pairs_offsite_nnnn_32=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==3,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==2,:)];

    pairs_offsite_nnnnn = result_matrices{5};
    pairs_offsite_nnnnn_13=result_matrices{5}(mod(result_matrices{5}(:, 1),2)~=0 & mod(result_matrices{5}(:, 2),2)~=0, :);
    pairs_offsite_nnnnn_24=result_matrices{5}(mod(result_matrices{5}(:, 1),2)==0 & mod(result_matrices{5}(:, 2),2)==0, :);


    % Construct pairs for U and V
    pairsU0 = {pairs_onsite_1, pairs_onsite_2, pairs_onsite_3, pairs_onsite_4};
    pairsU = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34,...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
    pairsV = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
end

function [xinitial, ni, si, efermi] = runHartreeFock(gs, xinitial, U, V, knum, stepmax, critial)
    % runHartreeFock Runs the Hartree-Fock self-consistent calculation
    %
    % Inputs:
    %   gs       - Struct, the supercell Hamiltonian and related properties
    %   xinitial - Cell array, initial guesses for Hartree-Fock states
    %   U        - Array, on-site Coulomb interaction values
    %   V        - Array, off-site Coulomb interaction values
    %   knum     - Integer, number of k-points along each direction
    %   stepmax  - Integer, maximum number of Hartree-Fock steps
    %   critial  - Float, convergence criterion for self-consistency
    %
    % Outputs:
    %   xinitial - Cell array, converged Hartree-Fock states
    %   ni       - Array, converged occupation numbers
    %   si       - Array, spin polarization values
    %   efermi   - Float, Fermi energy of the system

    % Generate 2D k-mesh
    [Kx, Ky, Kz] = gs.get_Bulk2Dkmesh([0, 1], [0, 1], knum);
    kpoints = [Kx(:), Ky(:), Kz(:)];

    % Set electric field to zero (default)
    Electric_field_in_evpA = 0;

    % Run the Hartree-Fock solver
    [xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, ...
        Electric_field_in_evpA, xinitial, stepmax, critial, U, V, 1e-10, 4);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Function to find the n-th NN neighbor pairing          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result_matrix = find_neighbor_data(coords_cartesian, lattice_vectors, n, dimensionality)
    % 动态处理二维或三维晶格，计算第 n 近邻的原子对信息
    %
    % 输入:
    % coords_cartesian: 原胞内的笛卡尔坐标 (n_atoms x 3)
    % lattice_vectors: 晶格基矢量 (3 x 3)
    % n: 第 n 近邻
    % dimensionality: 2 表示二维晶格，3 表示三维晶格
    %
    % 输出:
    % result_matrix: (num_pairs x 6) 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 1. 将笛卡尔坐标转换为分数坐标
    inv_lattice = inv(lattice_vectors);
    coords_fractional = coords_cartesian * inv_lattice;

    % 2. 构造超胞
    if dimensionality == 2
        % 仅扩展 x 和 y
        [super_coords, super_indices] = construct_supercell_2d(coords_fractional, lattice_vectors, 1);
    elseif dimensionality == 3
        % 扩展 x, y, z
        [super_coords, super_indices] = construct_supercell(coords_fractional, lattice_vectors, 1);
    else
        error('Dimensionality must be 2 or 3');
    end

    % 3. 在 Non-PBC 条件下计算超胞的距离矩阵
    if dimensionality == 2
        % 仅考虑 x-y 平面距离
        dist_matrix_non_pbc = compute_distance_matrix_2d(super_coords, lattice_vectors, false);
    elseif dimensionality == 3
        % 考虑完整三维距离
        dist_matrix_non_pbc = compute_distance_matrix(super_coords, lattice_vectors, false);
    end

    % 4. 提取 unique 距离（非零，考虑浮点误差）
    raw_distances = triu(dist_matrix_non_pbc);
    % raw_distances = raw_distances(raw_distances > 0); % 去掉 0 距离
    rounded_distances = round(raw_distances, 2); % 保留两位小数，分组距离
    unique_distances = unique(rounded_distances, 'sorted');

    % 5. 找到第 n 近邻的距离
    if n > length(unique_distances)
        error('第 %d 近邻超出最大可能的距离范围', n);
    end
    nth_distance = unique_distances(n);

    % 6. 找到满足第 n 近邻距离的原子对
    [pair_i, pair_j] = find(abs(dist_matrix_non_pbc - nth_distance) < 1e-2); % 容忍浮点误差

    % 7. 直接生成结果矩阵
    result_matrix=find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors);
end

function unique_pairs = find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors)
    % 找到所有第 n 近邻的原子对并去除重复
    %
    % 输入:
    % pair_i, pair_j: 满足第 n 近邻条件的原子对索引
    % super_coords: 超胞中原子的分数坐标
    % super_indices: 超胞中原子的索引和周期性偏移
    % nth_distance: 第 n 近邻的距离
    % lattice_vectors: 晶格基矢量
    %
    % 输出:
    % unique_pairs: 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 初始化结果存储
    num_pairs = length(pair_i);
    all_pairs = zeros(num_pairs, 6); % [i, j, frax, fray, fraz, dis]

    % 遍历所有原子对
    for k = 1:num_pairs
        i = pair_i(k); % 原子 i
        j = pair_j(k); % 原子 j

        % 原胞内的原子编号
        atom_i = super_indices(i, 1);
        atom_j = super_indices(j, 1);

        % 计算周期性偏移（分数坐标差）
        delta_r = super_coords(i, :) - super_coords(j, :);

        % 将分数坐标差转换为笛卡尔坐标，用于计算实际距离
        delta_cartesian = delta_r * lattice_vectors;
        distance = sqrt(sum(delta_cartesian.^2));

        % 检查是否接近第 n 近邻距离
        if abs(distance - nth_distance) < 1e-1
            % 记录原子对信息
            all_pairs(k, :) = [atom_i, atom_j, super_indices(j,2:end)-super_indices(i,2:end), nth_distance];
        end
    end

    % 去除重复的原子对（如 i->j 和 j->i）
    % 按原子对的编号排序，并提取唯一值
    [~, unique_rows] = unique(all_pairs, 'rows');
    unique_pairs = all_pairs(unique_rows, :); % 提取唯一行
end

function dist_matrix = compute_distance_matrix(coords, lattice_vectors, pbc)
    % 计算三维距离矩阵，可选是否使用周期性边界条件 (PBC)
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            delta_cartesian = delta_r * lattice_vectors; % 转为笛卡尔坐标
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function dist_matrix = compute_distance_matrix_2d(coords, lattice_vectors, pbc)
    % 计算二维距离矩阵（仅考虑 x 和 y 方向），可选 PBC
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            % 仅保留 x 和 y 方向的笛卡尔坐标
            delta_cartesian = delta_r * lattice_vectors; 
            delta_cartesian = delta_cartesian(:, 1:2); % x 和 y 方向
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function [super_coords, super_indices] = construct_supercell(coords, lattice_vectors, scale)
    % 构造三维超胞
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            for k = -scale:scale
                offset = [i, j, k];
                super_coords = [super_coords; coords + offset];
                for a = 1:n_atoms
                    super_indices = [super_indices; a, i, j, k];
                end
            end
        end
    end
end

function [super_coords, super_indices] = construct_supercell_2d(coords, lattice_vectors, scale)
    % 构造二维超胞（仅扩展 x 和 y 方向）
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            offset = [i, j, 0]; % z 方向固定为 0
            super_coords = [super_coords; coords + offset];
            for a = 1:n_atoms
                super_indices = [super_indices; a, i, j, 0];
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   Function to run hartree                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xinitial_U, xinitial_V] = initializeStates(U, V, pairsU, pairsV, nbands,xinitial_0)

    % initializeStates initializes the xinitial_U and xinitial_V arrays
    % Inputs:
    %   U       - Cell array or vector for U
    %   V       - Cell array or vector for V
    %   pairsU  - Cell array of pair matrices corresponding to U
    %   pairsV  - Cell array of pair matrices corresponding to V
    %   nbands  - Number of bands
    % Outputs:
    %   xinitial_U - Cell array containing initialized 3D matrices for U
    %   xinitial_V - Cell array containing initialized 3D matrices for V

    % Initialize empty cell arrays
    xinitial_U = {};
    xinitial_V = {};
    
    % Populate xinitial_U
    if ~isempty(U)
        for i = 1:length(U)
            pair = pairsU{i}; % Get the (i-1)th pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_idx = 1:num_pair
                % 提取原子对索引和分数坐标
                atom_i = pair(pair_idx, 1); % 原子 i
                atom_j = pair(pair_idx, 2);
                x=U(i);
                nj=xinitial_0(2*atom_j-1,2*atom_j-1)+xinitial_0(2*atom_j,2*atom_j);
                correlation(pair_idx,2*atom_i-1, 2*atom_i-1) = x*nj;%
                correlation(pair_idx,2*atom_i, 2*atom_i) = x*nj;%
            end
            xinitial_U = [xinitial_U, correlation];
        end
    end
    
    % Populate xinitial_V
    if ~isempty(V)
        for i = 1:length(V)
            pair = pairsV{i}; % Get the ith pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_index = 1:num_pair
                % 提取原子对索引  
                atom_i = pair(pair_index, 1); % 原子 i
                atom_j = pair(pair_index, 2); % 原子 j
                row_up_i = 2 * atom_i - 1;
                row_dn_i = 2 * atom_i;
                col_up_j = 2 * atom_j - 1;
                col_dn_j = 2 * atom_j;
                x=V(i)*0;
                correlation(pair_index,row_up_i, col_up_j) = -x; % 上自旋部分
                correlation(pair_index,row_dn_i, col_dn_j) = -x; % 下自旋部分
                correlation(pair_index,row_dn_i, col_up_j) = 0;
                correlation(pair_index,row_up_i, col_dn_j) = 0;
                correlation(pair_index,col_up_j, row_up_i) = -x; % 上自旋部分
                correlation(pair_index,col_dn_j, row_dn_i) = -x; % 下自旋部分
                correlation(pair_index,col_dn_j, row_up_i) = 0;
                correlation(pair_index,col_up_j, row_dn_i) = 0;

            end
            xinitial_V = [xinitial_V, correlation];
        end
    end
end


function potentials = calculate_dual_gate_potentials(r, d, epsilon, n_max)
    % calculate_dual_gate_potentials: Compute U and V values for various distances.
    % 
    % Inputs:
    %   d       - Distance between the gates (in meters)
    %   epsilon - Dielectric constant
    %   n_max   - Maximum number density
    % 
    % Output:
    %   potentials - Struct containing distances and their corresponding potentials

    % Define atomic distances (in meters)
    distances = [r, 3.77, 4.4, 6.91, 7.54, 8.66, 10.17] * 10^(-10);

    % Initialize results
    potentials = struct();
    potentials.distances = distances;
    potentials.values = zeros(size(distances));

    % Calculate dual gate potential for each distance
    for i = 1:length(distances)
        potentials.values(i) = dual_gate_potential(distances(i), d, epsilon, n_max);
    end

    % Display results
    fprintf('Distance (m) \t Potential (V)\n');
    for i = 1:length(distances)
        fprintf('%.2e \t %.2f\n', distances(i), potentials.values(i));
    end
end


function V = dual_gate_potential(r, d, epsilon, n_max)
    % Calculate dual-gate Coulomb potential using image charge method
    % r: radial distance
    % d: distance to the gates
    % q: charge
    % epsilon: dielectric constant
    % n_max: number of image charges to consider
    e=1.602176634*10^(-19);
    epsilon_0=8.854187817*10^(-12);
    k_e=1/(4*pi*epsilon_0);
    V = 0; % Initialize potential
    for n = -n_max:n_max
        V = V + ((-1)^n) / sqrt(r^2 + (2*n*d)^2);
    end
    % V = e^2/(4*pi*epsilon_0)/e/epsilon * V; % Final potential meV*m
    V=e^2*k_e*1/epsilon/e * V; %to eV
end

function [xinitial,ni,si,efermi]=runhartreev8(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)

        [xinitial, metaData] = flattenNestedCell(xinitial);        % 展平操作
        objective = @(xinitial) one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV);
        xinitial = quasi_newton(objective, xinitial, 7, critial, stepmax);
        xinitial = restoreNestedCell(xinitial, metaData);  % 复原操作
        modifyHam(gs, xinitial, U, V, pairsU, pairsV)  %修改Ham
        nbands=size(gs.ham,1);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function xnew = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)
    % 1. 基本初始化
    gs.ham=gs.iniham+0; 
    nbands = size(gs.ham, 1);
    % 恢复操作
    xinitial = restoreNestedCell(xinitial, metaData);
    xinitial_0=xinitial{1};
    if ~isempty(U)
         xinitial_U=xinitial{2};
    end 
    if ~isempty(V)
        xinitial_V=xinitial{3};
    end 

    gs.onsite_modify(xinitial_0);

    if ~isempty(U)
        for i = 1:length(U)
            pairs = pairsU{i};
            xinitial_1 = xinitial_U{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_1(j,:,:)));
            end
        end
    end

    if ~isempty(V)
        for i = 1:length(V)
            pairs = pairsV{i};
            xinitial_2 = xinitial_V{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_2(j,:,:)));
            end
        end
    end


    
    % 2. 计算波函数和能量 (含电场影响)
    [Unk, Enk] = MTB.ham.get_bulk_plane_bands_add_electric(gs, Electric_field_in_evpA, Kx, Ky, Kz);
    Unk = reshape(Unk, [nbands, nbands, knum^2]);
    Enk = reshape(Enk, [knum^2, nbands]);
    C_k_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);
    C_k = calculate_C_k_u1u2(Unk, Enk, u1, u2);

    % 3. 计算 U 和 V 的关联函数
    onsite_correlation_U =  zeros(size(xinitial_0));
    offsite_correlation_U = {};
    offsite_correlation_V = {};
    total_correlation = {};
   
    % 3.1 计算 U[0] 的关联函数 (onsite)
    if ~isempty(U0)
        for i = 1:length(U0)
            onsite_U = calculate_onsite_matrix(C_k_avg,pairsU0{i}, U0(i));
            % onsite_U = onsite_U*0.8+xinitial_0*0.2;
            onsite_correlation_U = onsite_correlation_U + onsite_U;
        end
        % To keep the TR symmetry
        % % diag_terms = diag(onsite_correlation_U)+diag(C_k_avg);
        % % onsite_correlation_U(1:nbands+1:end)=diag_terms/2;
        % % onsite_correlation_U=(onsite_correlation_U+onsite_correlation_U')/2;
    end
    
    % 3.2 计算 U[1:] 的关联函数 (offsite)
    if ~isempty(U)
        for i = 1:length(U)           
            offsite_U = calculate_offsite_U(C_k, kpoints, pairsU{i}, gs.a, U(i));
            % offsite_U = offsite_U*0.8 + xinitial{2}{i}*0.2;
            offsite_correlation_U = [offsite_correlation_U,offsite_U];
        end
    end
    
    % 3.3 计算 V 的关联函数 (offsite_V)
    if ~isempty(V)
        for i = 1:length(V)
            offsite_V = calculate_offsite_V(C_k, kpoints, pairsV{i}, gs.a, V(i));
            % offsite_V = offsite_V*0.8 + xinitial{3}{i}*0.2;         
            offsite_correlation_V = [offsite_correlation_V, offsite_V];
        end
    end
    total_correlation = {onsite_correlation_U,offsite_correlation_U,offsite_correlation_V};

    [xnew, ~] = flattenNestedCell(total_correlation);% 展平成列向量
end

function modifyHam(gs, xinitial, U, V, pairsU, pairsV)
    % modifyStates modifies onsite and offsite states using gs object
    % Inputs:
    %   gs        - Object containing methods `onsite_modify` and `offsite_modify`
    %   xinitial  - Cell array containing initial state matrices
    %   U         - Cell array or vector for U
    %   V         - Cell array or vector for V
    %   pairsU    - Cell array of pair matrices corresponding to U
    %   pairsV    - Cell array of pair matrices corresponding to V

    % Handle onsite modification for xinitial{1}
    xinitial_0 = xinitial{1};
    gs.onsite_modify(xinitial_0);

    % Handle offsite modification for U
    if ~isempty(U)
        xinitial_U = xinitial{2}; % Extract second cell for U states
        for i = 1:length(U)
            pairs = pairsU{i};  % Get pairs for the (i-1)th entry
            xinitial_1 = xinitial_U{i}; % Get initial state for U
            for j = 1:size(pairs, 1)
                % Modify offsite states for U
                gs.offsite_modify([0,0,0], squeeze(xinitial_1(j, :, :)));
            end
        end
    end

    % Handle offsite modification for V
    if ~isempty(V)
        xinitial_V = xinitial{3}; % Extract third cell for V states
        for i = 1:length(V)
            pairs = pairsV{i};    % Get pairs for the ith entry
            xinitial_2 = xinitial_V{i}; % Get initial state for V
            for j = 1:size(pairs, 1)
                % Modify offsite states for V
                gs.offsite_modify(pairs(j, 3:5), squeeze(xinitial_2(j, :, :)));
            end
        end
    end
end

function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
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

function [Etot,kindex,bandindex,efermi]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    efermi = max(a,[],"all");
    % rule out the bands below the fermi level
    % % a=a(ceil(size(Enk(:),1)*0.5)+1:end);
    % % b=b(ceil(size(Enk(:),1)*0.5)+1:end);

    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end

function C_total_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)，每行是一个 k 点的能带
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    % u1: 填充比例下限 (0 <= u1 <= 1)
    % u2: 填充比例上限 (0 <= u2 <= 1)
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)

    % 获取维度信息
    [knum2, nbands] = size(Enk);  % Enk 的维度
    C_total = zeros(nbands, nbands); % 初始化总关联函数矩阵

    % 1. 计算 u1 和 u2 对应的费米能级
    E_low = calculate_ef(Enk, u1); % u1 对应的费米能级
    E_high = calculate_ef(Enk, u2); % u2 对应的费米能级

    %fprintf('能级范围: %.6f eV 到 %.6f eV\n', E_low, E_high);

    % 2. 遍历每个 k 点，累加所有 C^(m)
    for m = 1:knum2
        % 提取第 m 个 k 点的本征值和本征矢量
        E_k = Enk(m, :);          % 第 m 个 k 点的能带本征值，维度 (1, nbands)
        U_k = Unk(:, :, m);       % 第 m 个 k 点的本征矢量矩阵，维度 (nbands, nbands)

        % 筛选能级：处于 [E_low, E_high] 区间的态置为 1，其余置为 0
        W_k = diag(double(E_k >= E_low & E_k <= E_high)); % 维度 (nbands, nbands)

        % 计算关联函数矩阵 C^(m)
        C_m = conj(U_k) * W_k * transpose(U_k);

        % 累加当前 k 点的关联函数矩阵
        C_total = C_total + C_m;
    end

    % 3. 求平均：除以总的 k 点数目 knum^2
    C_total_avg = C_total / knum2;

    %fprintf('关联函数计算完成，结果已对 k 点求平均。\n');
end

function C_k = calculate_C_k_u1u2(U_k, E_k, u1, u2)
    % 计算关联函数矩阵 C_k，使用 u1 和 u2 提取占据能级范围
    %
    % 输入:
    % U_k: 本征矢量矩阵，维度 (nbands, nbands, nk)
    % E_k: 本征值矩阵，维度 (nk, nbands)
    % u1, u2: 占据比例范围
    %
    % 输出:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)

    % 获取维度
    [nbands, ~, nk] = size(U_k);
    C_k = zeros(nbands, nbands, nk);

    % 确定占据能级范围 [E_low, E_high]
    total_states = numel(E_k);
    E_sorted = sort(E_k(:)); % 将所有能级排序
    num_u1 = ceil(total_states * u1);
    num_u2 = ceil(total_states * u2);
    E_low = E_sorted(num_u1);
    E_high = E_sorted(num_u2);

    % 遍历每个 k 点
    for k = 1:nk
        % 提取当前 k 点的本征值和本征矢量
        E_k_point = E_k(k, :);
        U_k_point = U_k(:, :, k);

        % 构建占据函数 f(E_k) 范围
        f_k = double(E_k_point >= E_low & E_k_point <= E_high);

        % 计算关联函数矩阵 C_k
        C_k(:, :, k) = conj(U_k_point) * diag(f_k) * transpose(U_k_point);
    end
end

function onsite_U = calculate_onsite_matrix(correlation, pairs, U)

    % 计算 Onsite 矩阵，仅保留对角线上的 2x2 block 并进行变换
    %
    % 输入:
    % correlation: 关联函数矩阵，维度 (2*nbands, 2*nbands)
    % U: Hubbard U 参数
    %
    % 输出:
    % onsite_matrix: Onsite 矩阵，维度 (2*nbands, 2*nbands)

    % 获取轨道数目
    total_bands = size(correlation, 1);
    if mod(total_bands, 2) ~= 0
        error('关联函数矩阵的行数必须是偶数，每个轨道对应上下自旋');
    end
    nbands = total_bands / 2;

    % 定义变换矩阵 T
    T = [0, -1; 1, 0];

    % 初始化输出矩阵
    onsite_matrix = zeros(total_bands, total_bands);
    num_pairs = size(pairs,1);

    % 遍历每个对角线的 2x2 block
       % 遍历每个原子对
    for pair_idx = 1:num_pairs
        i = pairs(pair_idx,1);
        % 提取对角线上的 2x2 block
        block = correlation(2*i-1:2*i, 2*i-1:2*i);

        % 应用变换矩阵 T
        transformed_block = T * block * T';

        % 保存到输出矩阵的对应位置
        onsite_matrix(2*i-1:2*i, 2*i-1:2*i) = U * transformed_block;
        % To keep the TR symmetry
    end
    onsite_U=onsite_matrix;
end

function offsite_U = calculate_offsite_U(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_U = zeros(num_pairs, nbands, nbands);
    % offsite_U_total = zeros(nbands,nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义行/列索引
        index_i = 2 * atom_i - 1; 
        index_j = 2 * atom_j - 1; 


        % 遍历所有 k 点，计算相位修正的关联值
        ni_up_sum = 0;
        ni_dn_sum = 0;
        nj_up_sum = 0;
        nj_dn_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

            % 提取 C_k 中的元素
            % ni_up = C_k(index_i, index_i, k) * phase_factor; % 上自旋
            % ni_dn = C_k(index_i+1, index_i+1, k) * phase_factor; % 下自旋
            nj_up = C_k(index_j, index_j, k) * phase_factor; % 上自旋
            nj_dn = C_k(index_j+1, index_j+1, k) * phase_factor; % 下自旋
            % 累加相位修正的贡献
            % ni_up_sum = ni_up_sum + ni_up;
            % ni_dn_sum = ni_dn_sum + ni_dn;
            nj_up_sum = nj_up_sum + nj_up;
            nj_dn_sum = nj_dn_sum + nj_dn;
        end

        % 平均化，并乘以 U
        % ni_up_sum_avg = real(ni_up_sum/nk) * U;
        % ni_dn_sum_avg = real(ni_dn_sum/nk) * U;
        nj_up_sum_avg = real(nj_up_sum/nk) * U;
        nj_dn_sum_avg = real(nj_dn_sum/nk) * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(index_i, index_i) = nj_dn_sum_avg+nj_up_sum_avg; %
        correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg+nj_dn_sum_avg; %
        % correlation_temp(index_j, index_j) = ni_dn_sum_avg+ni_up_sum_avg; %
        % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg+ni_dn_sum_avg; %
        % To keep the PT symmetry
        % % % correlation_temp(index_i, index_i) = nj_dn_sum_avg/2+nj_up_sum_avg/2; %
        % % % correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg/2+nj_dn_sum_avg/2; %
        % % % correlation_temp(index_j, index_j) = ni_dn_sum_avg/2+ni_up_sum_avg/2; %
        % % % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg/2+ni_dn_sum_avg/2; %       
        % 将当前原子对的矩阵存储到输出列表中
        offsite_U(pair_idx, :, :) = correlation_temp;
        % offsite_U_total = offsite_U_total+correlation_temp;
    end
end

function offsite_V = calculate_offsite_V(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_V = zeros(num_pairs, nbands, nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义上下自旋的行/列索引
        col_dn_i = 2 * atom_i;       % 下自旋 (dn) 行
        col_up_i = 2 * atom_i - 1;   % 上自旋 (up) 行
        row_dn_j = 2 * atom_j;       % 下自旋 (dn) 列
        row_up_j = 2 * atom_j - 1;   % 上自旋 (up) 列

        % 遍历所有 k 点，计算相位修正的关联值
        dndn_sum = 0;
        upup_sum = 0;
        updn_sum = 0;
        dnup_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

            % 提取 C_k 中的元素
            dndn = C_k(row_dn_j, col_dn_i, k) * phase_factor; % 下自旋-下自旋
            upup = C_k(row_up_j, col_up_i, k) * phase_factor; % 上自旋-上自旋
            updn = C_k(row_up_j, col_dn_i, k) * phase_factor;
            dnup = C_k(row_dn_j, col_up_i, k) * phase_factor;

            % 累加相位修正的贡献
            dndn_sum = dndn_sum + dndn;
            upup_sum = upup_sum + upup;
            updn_sum = updn_sum + updn;
            dnup_sum = dnup_sum + dnup;
        end

        % 平均化，并乘以 U
        dndn_avg = dndn_sum / nk * U;
        upup_avg = upup_sum / nk * U;
        updn_avg = updn_sum / nk * U;
        dnup_avg = dnup_sum / nk * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(col_up_i, row_up_j) = -upup_avg; % 上自旋部分
        correlation_temp(col_dn_i, row_dn_j) = -dndn_avg; % 下自旋部分
        correlation_temp(col_dn_i, row_up_j) = -updn_avg; % 上自旋部分
        correlation_temp(col_up_i, row_dn_j) = -dnup_avg; % 下自旋部分

        % 将当前原子对的矩阵存储到输出列表中
        offsite_V(pair_idx, :, :) = correlation_temp;
    end
end

function projector = get_projector(gs, kpoints,Unk, whichbands)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % gs: geometry
    % kpoints: 面内mesh 维度(knum^2,3)
    % whichbands: 选择投影的band 维度(1,nbands)
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    %
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)
    [nbands,~,knum2]=size(Unk);
    projectors_k = zeros(nbands,nbands,knum2);
    % 构造 projectors_k
    for i = 1:knum2
        projectors_k(:, :, i) = diag(which_bands); % 按对角阵填充
    end

    % 应用投影算符
    vs_conj = conj(Unk);
    vs_transpose = permute(Unk, [1, 2, 1]); % 转置第一和第二维
    projectors_k = pagemtimes(pagemtimes(vs_conj, projectors_k), vs_transpose);


    rvectors=gs.hopr;
    expr=exp(1j*(kpoints*(rvectors*gs.a).'));
    % 获取维度信息
    [knum2, nbands] = size(Enk);  % Enk 的维度

    
    C_total = zeros(nbands, nbands); % 初始化总关联函数矩阵

    % 1. 计算 u1 和 u2 对应的费米能级
    E_low = calculate_ef(Enk, u1); % u1 对应的费米能级
    E_high = calculate_ef(Enk, u2); % u2 对应的费米能级

    %fprintf('能级范围: %.6f eV 到 %.6f eV\n', E_low, E_high);

    % 2. 遍历每个 k 点，累加所有 C^(m)
    for m = 1:knum2
        % 提取第 m 个 k 点的本征值和本征矢量
        E_k = Enk(m, :);          % 第 m 个 k 点的能带本征值，维度 (1, nbands)
        U_k = Unk(:, :, m);       % 第 m 个 k 点的本征矢量矩阵，维度 (nbands, nbands)

        % 筛选能级：处于 [E_low, E_high] 区间的态置为 1，其余置为 0
        W_k = diag(double(E_k >= E_low & E_k <= E_high)); % 维度 (nbands, nbands)

        % 计算关联函数矩阵 C^(m)
        C_m = conj(U_k) * W_k * transpose(U_k);

        % 累加当前 k 点的关联函数矩阵
        C_total = C_total + C_m;
    end

    % 3. 求平均：除以总的 k 点数目 knum^2
    C_total_avg = C_total / knum2;

    %fprintf('关联函数计算完成，结果已对 k 点求平均。\n');
end


function [band_gap, VBM, CBM, VBM_k, CBM_k] = compute_band_gap_1d(En, valence_band_index, conduction_band_index)
% compute_band_gap_1d: 计算指定价带和导带之间的能隙及对应的 k 点位置
%
% 输入参数:
%   En - 能量矩阵，维度为 (knum, nbands)，包含能量数据
%   valence_band_index - 价带索引 (整数)
%   conduction_band_index - 导带索引 (整数)
%
% 输出参数:
%   band_gap - 带隙 (CBM - VBM)
%   VBM - 价带最高点能量值
%   CBM - 导带最低点能量值
%   VBM_k - 价带最高点对应的 k 点位置 (整数)
%   CBM_k - 导带最低点对应的 k 点位置 (整数)

    % 获取输入矩阵的维度
    [knum, nbands] = size(En);
    assert(valence_band_index <= nbands && conduction_band_index <= nbands, ...
        '价带或导带索引超出能量矩阵的范围');

    % Step 1: 计算 VBM 和 CBM 及其索引
    [VBM, VBM_k] = max(En(:, valence_band_index)); % 价带最高点
    [CBM, CBM_k] = min(En(:, conduction_band_index)); % 导带最低点

    % Step 2: 计算带隙
    band_gap = CBM - VBM;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Function for scf convergency                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xn = quasi_newton(func, x0, q, tol, maxsteps, history)
    % 默认参数设置
    if nargin < 3, q = 10; end
    if nargin < 4, tol = 1e-7; end
    if nargin < 5, maxsteps = 100; end
    if nargin < 6, history = {}; end

    % 初始化变量
    f0 = func(x0).'; % 初始梯度，确保为列向量
    n = length(f0); % 确保 x0 和 func 输出长度一致

    alpha = 1e-3; % 初始步长
    U = zeros(n, q); % 初始化 U 矩阵
    V = zeros(n, q); % 初始化 V 矩阵

    xn = x0(:); % 确保 x0 是列向量
    fn = f0; % 初始化梯度为列向量

    disp('Starting quasi-Newton...');

    % === 预热阶段：初始化 U 和 V ===
    for i = 1:q
        fn = func(xn).'; % 确保梯度为列向量
        xn_1 = xn - alpha * fn; % 初步更新变量

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 更新 U 和 V 矩阵
        U(:, i) = fn - xn; % 梯度差
        V(:, i) = func(fn).' - fn; % 梯度映射差
        xn = xn_1; % 更新变量
    end

    % % % % % % disp('Warmed up.');

    % === 主迭代循环 ===
    for i = 1:maxsteps
        fn = func(xn).'; % 确保梯度为列向量

        % 检查维度一致性
        xn = xn(:); % 确保 xn 是列向量
        fn = fn(:); % 确保 fn 是列向量
        if size(U, 1) ~= length(xn)
            error('Dimension mismatch: U should have %d rows, but got %d.', length(xn), size(U, 1));
        end

        % 构造矩阵 C 和向量 b
        C = transpose(U) * U - transpose(U) * V; % 计算 C 矩阵
        b = transpose(U) * (xn - fn); % 计算 b 向量

        % 判断 C 是否病态，必要时正则化
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        % 求解更新方向 delta
        delta = V * (C \ b);

        % 更新变量
        xn_1 = fn - delta;

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 存储历史记录
        history{end+1} = xn_1;

        % 更新 U 和 V 矩阵
        fn_1 = func(xn_1).'; % 确保梯度为列向量
        U = [U(:, 2:end), fn_1 - xn_1]; % 滑动窗口更新 U
        V = [V(:, 2:end), func(fn_1).' - fn_1]; % 滑动窗口更新 V

        % 检查收敛条件
        dx = max(abs(xn_1 - xn)); % 计算变量更新幅度
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Converged at iteration %d\n', i);
            xn = xn_1.'; % 最终结果转置为行向量
            return;
        end

        xn = xn_1; % 更新当前变量
    end

    % 未收敛提示
    disp('Did not converge within the maximum steps.');
    xn = xn.'; % 如果未收敛，也转置为行向量
end

function [flattenedVector, metaData] = flattenNestedCell(C)
    flattenedVector = []; % 用于存储展平后的行向量
    metaData = {};        % 用于存储路径和尺寸信息

    % 递归处理单元数组
    function processElement(element, path)
        if iscell(element)
            % 如果是单元数组，递归处理每个子元素
            for i = 1:numel(element)
                processElement(element{i}, [path, i]);
            end
        else
            % 如果是普通数组，展平并记录信息
            flattenedVector = [flattenedVector, element(:)']; % 展平并拼接
            metaData{end+1} = struct('Path', {path}, 'Size', size(element)); % 记录路径和尺寸
        end
    end

    % 开始递归处理
    processElement(C, []);
end

function restoredCell = restoreNestedCell(flattenedVector, metaData)
    restoredCell = {};    % 初始化恢复后的单元数组
    currentIndex = 1;     % 当前展平向量的索引

    % 递归赋值函数
    function target = assignElement(target, path, sizeInfo)
        if isempty(path)
            % 如果路径为空，说明已经到达最终节点
            numElements = prod(sizeInfo); % 当前数组的总元素数
            reshapedArray = reshape(flattenedVector(currentIndex:currentIndex + numElements - 1), sizeInfo); % 恢复形状
            currentIndex = currentIndex + numElements; % 更新索引位置
            target = reshapedArray; % 直接返回恢复后的数组
        else
            % 处理嵌套单元数组
            idx = path(1);
            if numel(target) < idx || ~iscell(target{idx})
                target{idx} = {}; % 确保路径上的每一级是单元数组
            end
            target{idx} = assignElement(target{idx}, path(2:end), sizeInfo); % 递归处理下一层路径
        end
    end

    % 遍历元数据并恢复每个元素
    for i = 1:numel(metaData)
        restoredCell = assignElement(restoredCell, metaData{i}.Path, metaData{i}.Size);
    end
end

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
                if abs(Hij - Hji_conj) > 1e-5
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