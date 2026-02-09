clear;
clear all;
% parpool('local',4)

g = MTB.geometry("SrSnO");
g = MTB.read_poscar(g,"data/SrSnO/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/SrSnO/wannier90_hr_p1.dat","data/SrSnO/wannier90_hr_p2.dat");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%

nslab=300;
delta=0.03;% 0.03~0.05
numEigs=30;
efermi=3.6616;
mus=linspace(3.3616,3.9616,101);%mu for E_f-mu to E_f+mu of 100 points
kpoint=[0.0,0.0];% Gamma Point
tic;
Energy=MTB.ham.get_slab_mu_E_sparse_BdG(g.ham,g.hopr2,nslab,nbands,numEigs,nrpts,kpoint,g.a2,mus,delta);
save("Energy-mus.mat","Energy");
toc
%%

x=repmat(mus,30,1)
load("Energy-mus.mat")
% plot(x',Energy','*-','Color','red')
plot(x'-efermi,Energy','*-')
%%
load("data/SrSnO/1000s/SrSnO-Energy-mus-1000s-501mu-003.mat")
figure()
mus=linspace(3.3616,3.9616,501);
x=repmat(mus,30,1)
plot(x'-3.6616,Energy','o')