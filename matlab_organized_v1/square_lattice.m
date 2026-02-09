clear;
clear all;
delta=2.8;
%get the Kai 3band Ham
g=get_Kai_3band(delta);

% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Calculate the BC in xy plane            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=51;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
% Calculate Berry Curvature by LOOP method
plottap=1;
bandindex=3;
[Omega_k,KX,KY] = MTB.ham.get_Berry_curvature(bandindex,Unk,Kx,Ky,plottap);
Chern=sum(Omega_k,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Get the Wilson Loop              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=51;
kx=linspace(0,1,knum);
band1=3;
band2=3;
[wx1,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx2,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Get slab bands for edges           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=get_Kai_3band(delta);
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=51;
%
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'};
hkpoints={[-0.5,0.0],...
          [0.0,0.0],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0);


%%
function g=get_Kai_3band(delta)

g = MTB.geometry("square");

%set lattice constant
g.a=[1,0,0;...
     0,1,0;...
    0,0,10];
%get rec lattce
g.b=inv(g.a)*2*pi;
%set atoms positions
g.atoms=[0,0,0];
%set wannier orbiral positions
g.wpos=[0,0,0;...
        0,0,0;...
        0,0,0];

tdd = 1;
tpd = 1;
tpp = 1;
% delta=2.8;
onsite_d=-4*tdd+2*tpp+delta-2*tpp*delta/(4*tpp+delta);
tpp_prime=tpp*delta/(4*tpp+delta);

% Initialize Hamiltonian matrix elements based on tight-binding model
g.ham=zeros(3,3,5);
g.hopr=zeros(5,3);

g.hopr(1,:)=[0,0,0];
g.hopr(2,:)=[1,0,0];
g.hopr(3,:)=[-1,0,0];
g.hopr(4,:)=[0,1,0];
g.hopr(5,:)=[0,-1,0];

g.ham(:,:,1)=[onsite_d,   0   ,      0;...
                     0,   0   ,      1j*delta;...
                     0,   -1j*delta  ,       0];
g.ham(:,:,2)=[-tdd, tpd, 0;...
              -tpd, tpp, 0;...
                 0,   0, -tpp_prime];
g.ham(:,:,3)=[-tdd, -tpd, 0;...
              tpd, tpp, 0;...
                 0,   0, -tpp_prime];
g.ham(:,:,4)=[-tdd, 0, tpd;...
              0, -tpp_prime, 0;...
              -tpd,   0, tpp];
g.ham(:,:,5)=[-tdd, 0, -tpd;...
              0, -tpp_prime, 0;...
              tpd,   0, tpp];

end
