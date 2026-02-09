clear;
clear all;
% parpool('local',4)

g = MTB.geometry("SrSnO");
g = MTB.read_poscar(g,"data/SrSnO/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/SrSnO/wannier90_hr_p1.dat","data/SrSnO/wannier90_hr_p2.dat");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
%%

[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5]};% hkpoints-high symmetry k points
nk=51;
nslab=50;
efermi=3.6616; %% set Fermi Level
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"SrSnO-100slab")


