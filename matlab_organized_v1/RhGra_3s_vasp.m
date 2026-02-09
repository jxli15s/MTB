% clear;
% clear all;
g = MTB.geometry("Rgra_3s");
g = MTB.read_poscar(g,"data/Graphene/3s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/3s/wannier90_hr_p1.dat','data/Graphene/3s/wannier90_hr_p2.dat');
%%
%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'\Gamma','K','M'}; % labels for k
hkpoints={[0.0,0.0,0.000000],...
          [1/3,2/3,0.000000],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=-2.8077 ;
%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)
%% Calculate Berry Curvature by LOOP method
knum=501;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
plottap=1;
bandindex=1:3;
[Omega_k,KX,KY] = MTB.ham.get_Berry_curvature(bandindex,Unk,Kx,Ky,plottap);