clc;
clear;
g = MTB.geometry("Square-octagon");
a=3.47; b=3.47; c=10;
g.a=[a,0,0;0,b,0;0,0,c];
a1=g.a(1,:);
a2=g.a(2,:);
a3=g.a(3,:);
omega=dot(a1,cross(a2,a3));
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
b=[b1;b2;b3];
g.b = b;

g.atoms=[0.25, 0.00, 0.50;...
         0.00, 0.25, 0.50;...
        -0.25, 0.00, 0.50;...
         0.00,-0.25, 0.50];
g.wpos=g.atoms*g.a;
t1=-2;
t2=3;
nrps=5;
g.ham=zeros(size(g.wpos,1),size(g.wpos,1),nrps);
g.hopr=zeros(nrps,3);
%%
g.hopr(1,:)=[0,0,0];
g.ham(1,2,1)=t2;
g.ham(2,1,1)=t2;
g.ham(2,3,1)=t2;
g.ham(3,2,1)=t2;
g.ham(3,4,1)=t2;
g.ham(4,3,1)=t2;
g.ham(1,4,1)=t2;
g.ham(4,1,1)=t2;
g.hopr(2,:)=[1,0,0];
g.ham(1,3,2)=t1;
g.hopr(3,:)=[0,1,0];
g.ham(2,4,3)=t1;
g.hopr(4,:)=[-1,0,0];
g.ham(3,1,4)=t1;
g.hopr(5,:)=[0,-1,0];
g.ham(4,2,5)=t1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','M','X'}; % labels for k
hkpoints={[0.5,0.0,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.5,0.0],...
          [0.5,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=101;
efermi=0;
%[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")






