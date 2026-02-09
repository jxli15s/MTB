clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/Graphene/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');

% [g.ham,g.hopr] = MTB.read_hr('data/Graphene/Graphene_hr.dat');
% g.atoms=[0.5,0.5,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[1,1];
g.get_suborbidx
%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'K','\Gamma','M','K'}; % labels for k
hkpoints={[0.333333,0.333333,0.000000],...
          [0.000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=-1.2533;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)
%%
Umatrix=[2,1,0;0,2,0;0,0,1];
shift=[0.0,0,0] %%
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"brute");
%%
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','Y','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.0,0.000000],...
          [0.000000000,0.0000000000,0.000000],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]...
          [0.000000000,0.0000000000,0.000000]};% hkpoints-high symmetry k points
% labels={'K','\Gamma','M','K'}; % labels for k
% hkpoints={[0.333333,0.333333,0.000000],...
%           [0.000000000,0.0000000000,0.000000],...
%           [0.5,0.0,0.0],...
%           [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=-1.2533;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Haldane-bulk",0)
% plot_geometry_sub(gs,g)
%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/Graphene/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');

% [g.ham,g.hopr] = MTB.read_hr('data/Graphene/Graphene_hr.dat');
% g.atoms=[0.5,0.5,0]
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[1,1];
g.get_suborbidx