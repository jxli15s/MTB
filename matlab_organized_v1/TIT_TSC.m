clear;
clear all;
%p=parpool(8)
% g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
% g.wpos=[];
% g.wpos=g.atoms*g.a;
% Step 1: Initialize Geometry and Hamiltonian
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);
%
[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
%
Electric_field_in_evpA=0.00;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Check Time Reversal Symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.2,0.4,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.4,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(4),1i*s2);

h1=inv(T)*conj(hk1)*T-hk2;
max(h1,[],'all')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Inversion Symmetry         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];
kpoint=[0.25,0.3,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.25,-0.3,0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
Pblk=blkdiag(d1,d2);
n1=1;
P=kron(fliplr(eye(n1)),Pblk);
%%
pos=eye(size(g.wpos,1));
kpoint=kpoint*g.b;
for i=1:size(g.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    % fprintf('%.6f\n',pos(i,i));
end
P=P*pos;
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
a=Psik'*P*Psik; %% P \psi=e \psi    e=\psi.T.conj P \psi
c=eig(a);
%%
a=eig(Psik(:,5:6)'*P*Psik(:,5:6));
sum(a)
%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(4,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
PnkZ2 = zeros(4,size(kpoints,1));
Z2 = zeros(size(kpoints,1),1);
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];
d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
P=blkdiag(d1,d2);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:)
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:4)'*P*Psik(:,1:4));
    Z4(i)=sum(Pnk(:,i));
    % PnkZ2=eig(Psik(:,1:2)'*P*Psik(:,1:2))
end
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate slab bands                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs = MTB.geometry("TaIrTe4");
gs = MTB.read_poscar(gs,"data/TaIrTe4_2d_tb/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
gs.wpos=gs.atoms*gs.a;

MillerIndices=[1,0,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=101;

[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.0],...
          [0.5000000000,0.0000000000],...
          [0.0,0.0]};% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Calculate surface states                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
gs = MTB.geometry("TaIrTe4");
gs = MTB.read_poscar(gs,"data/TaIrTe4_2d_tb/POSCAR");
[gs.ham,gs.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
gs.wpos=gs.atoms*gs.a;
MillerIndices=[1,0,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;

[nbands,~,nrpts]=size(gs.ham);
hkpoints={[0.0,0.0000000000],...
          [0.5000000000,0.0000000000],...
          [0.0,0.0000000000]};% hkpoints-high symmetry k points
labels={'\Gamma','X','\Gamma'};
Np=2;
omegamin=-1;
omegamax=1;
omeganum=500;
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
%%
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