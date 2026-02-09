clear;
clear all;
% parpool('local',4)

g = MTB.geometry("SrSnO");
g = MTB.read_poscar(g,"data/SrSnO/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr("data/SrSnO/wannier90_hr_p1.dat","data/SrSnO/wannier90_hr_p2.dat");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/SrSnO/upupdwdw/wannier90_hr_p1.dat","data/SrSnO/upupdwdw/wannier90_hr_p2.dat");
% MillerIndices=[1,1,1];
% Umatrix = g.MillerIndicestoumatrix(MillerIndices);
% Urot = g.surfab;
%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.5,0.5],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=3.6616; %% set Fermi Level
%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")


%%
kpoint=[0.0,0.0,0.0]
[Energy,Psik]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
%%
P=eye(24);
P(19:end,:)=-P(19:end,:);
a=Psik'*P*Psik
%%
s2=[0  -1i
    1i  0];
T=kron(eye(12),i*s2)
a=Psik'*T*Psik;
%%
kpoint=[0.0,0.0,0.0]
[Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
h1=inv(T)*conj(hk)*T