 clear;
 clear all;
 p=parpool('local',8)
%delete(p)
%% Read POSCAR and Wannier
 g = MTB.geometry("MnBiTe");
 g = MTB.read_poscar(g,"data/MnBiTe-xu/POSCAR-MnBiTe");
 [g.ham,g.hopr] = MTB.wannier.read_hr("data/MnBiTe-xu/wannier90_hr_p1.dat","data/MnBiTe-xu/wannier90_hr_p2.dat");
 [nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
 efermi=0.0258; %% set Fermi Level 0.0258

%% Calculate Berry Curvature by Quantum metric add electric fields

g.readwpos("data/MnBiTe-xu/wpos_MnBiTe");
g.wpos=g.wpos*27;

nk=5;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
step=4;

Enum=51;
Emin=-0.2;
Emax=0.2;
nsband=120:121
Electric_field_list=linspace(-0.003,0.003,step);
Enk_all=zeros(nk,nk,nbands,step);
Omega_k_all=zeros(nk,nk,size(nsband,2),step);

for i=1:step
    tic
    Electric_field_in_evpA=Electric_field_list(i);
    [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
    [Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk-efermi,Unk,nsband,0.03,0);
    Enk_all(:,:,:,i)=Enk-efermi;
    Omega_k_all(:,:,:,i)=Omega_k;
    toc
end
save("Omega_k_all.mat","Omega_k_all","-v7.3")
save("Enk_all.mat","Enk_all","-v7.3")

%% 
%% Calculate bulk bands
tic;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'K','\Gamma','K'}; % labels for k
hkpoints={[0.333333,0.333333,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [-0.333333,-0.333333,0.000000]};% hkpoints-high symmetry k points
nk=21;

[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
efermi=0.0258; %% set Fermi Level 0.0258
toc
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"MnBiTe-bulk",0)
toc
%%
Electric_field_in_evpA=0
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);