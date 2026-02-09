%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Construct the g.ham                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("Rgra_5s");
g = MTB.read_poscar(g,"data/Graphene/5s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/5s/wannier90_hr_p1.dat','data/Graphene/5s/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
g.Rcart = g.hopr * g.a;
%%
G = [0,0]+[1/3,2/3]*0.9; 
K=[1/3,2/3]; 
M=K+([0.5,0.5]-[1/3,2/3])*0.1;
nodes = [G;K;M];
labels = {'\Gamma','K','M','\Gamma'};
tbHFMF.plot_band_path(g, nodes, 80, labels);



%%
