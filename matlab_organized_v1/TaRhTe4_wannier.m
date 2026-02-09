clear;
clear all;
%parpool(8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Read the structure and hop from wannier  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = MTB.geometry("TaRhTe4");
g = MTB.read_poscar(g,"data/TaRhTe4/wannier/POSCAR-TaRhTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaRhTe4/wannier/wannier90_hr_p1.dat','data/TaRhTe4/wannier/wannier90_hr_p2.dat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate the bulk band structure      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
nk=51;
efermi=-1.6834;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-bulk",0)
hold on;
plot(kpath,Energy(60,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(61,:)-efermi,"Color",'red','LineWidth',2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   Calculate DOS                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("./data/TaRhTe4/data/TaIrTe4_bulk-400x400_wannier.mat")
nk=400;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
% [~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);


plottap=1;
Enum=10 1;
Emin=-0.2;
Emax=0.1;
eps=0.005;
efermi=-1.6834;

[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Plot the 3D bands for VHS               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("./data/TaRhTe4/data/bulk-plane-201-GammaCenter.mat")
% load("./data/TaRhTe4/data/TaIrTe4_bulk-400x400.mat")
load("./data/TaRhTe4/data/TaIrTe4_bulk-400x400_wannier.mat")
nk=400;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

figure('Color','white')
surf(Kx,Ky,Enk(:,:,76))
hold on;
surf(Kx,Ky,Enk(:,:,77))
colorbar;
% shading flat;
shading interp;
lighting gouraud;
 light('Position', [0 0 1], 'Style','infinite');
 material shiny;
colormap(slanCM('RdBu'))
hold on;

numContours = 1000;  % 
contourLevels = linspace(-0.3, 0.3, numContours);  %

contour3(Kx,Ky,Enk(:,:,77),contourLevels,'LineColor',[0.5 0.5 0.5]); 
contour3(Kx,Ky,Enk(:,:,76), contourLevels, 'LineColor',[0.5 0.5 0.5]); 
% contour3(Kx,Ky,Enk(:,:,76), 30, 'LineColor',[0.5 0.5 0.5]); 
zlim([-0.3 0.3]);
ylim([-0.2 0.2]);
caxis([-0.4 0.3])
axis 
view(-92,18)
% valence position 0.0 0.0506 0.0  -0.1407
% conduc 0.056 0.0177 0.0  0.01152
