clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Hexagonal");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

% gs = MTB.geometry("TaIrTe4_s");
%%
g.a=[1,0,0;...
    1/2,sqrt(3)/2,0;...
    0,0,10]
g.b=inv(g.a)*2*pi
g.atoms=[0,0,0]
g.wpos=[0,0,0;...
        0,0,0]
g.ham=zeros(2,2,7);
g.hopr=zeros(7,3);
for i=2:7
    g.ham(:,:,i)=[1,0;0,1];
end
g.hopr(1,:)=[0,0,0];
g.hopr(2,:)=[1,0,0];
g.hopr(3,:)=[-1,0,0];
g.hopr(4,:)=[0,1,0];
g.hopr(5,:)=[0,-1,0];
g.hopr(6,:)=[-1,1,0];
g.hopr(7,:)=[1,-1,0];
%%
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'M','\Gamma','K'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

%%
n1=1;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'M','\Gamma','K'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;
%%
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

%%

n1=3;
n2=3;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
% labels={'R','Y','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
labels={'M','\Gamma','K'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.333333,0.333333,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;

stepmax=10000;
minstepmax=11;
step=1;
knum=31;
u=(1*n1*n2)/nbands;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
U=20;

xinitial=zeros(nbands,nbands);
nsite=nbands/2;
for i = 1:nsite
    % nup=randn(1); 
    % ndn=1-nup;
    nup=0.5;
    ndn=0.5;
    sx=randn(1);
    sy=1-sx;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U;
end
% save("tmp.dat","tmp")
 % load("tmp.mat")
 % xinitial=tmp;

[xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);


gs.onsite_modify(xinitial);
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(1*n1*n2,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(1*n1*n2+1,:),"Color",'red','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);


% knum=81;
% band1=1;
% band2=1*n1*n2;
% [wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
spin_plot(gs,si)
%%
lattice_plot(gs)
%%
T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
%%
g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)]
% g.wpos=g.wpos
%%
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.2,0.2,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.2,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(4),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")



%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

gs = MTB.geometry("TaIrTe4_s");
%%
T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;

g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)]
% g.wpos=g.wpos
%%
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);

Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);

[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;
bandindex=1:120;


[Energy,Omega_k,kpath,kindex]=MTB.ham.get_bulk_bands_bcd(gs,hkpoints,nk,bandindex);

%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
% hold on;
% colormap(slanCM('RdBu'))
% scatter(kpath,Energy(62,:),[],(Omega_k(:,61)+Omega_k(:,62)),"filled");
% scatter(kpath,Energy(63,:),[],(Omega_k(:,63)+Omega_k(:,64)),"filled");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint=[0.2,0.2,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.2,-0.2,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
s2=[0  -1i
    1i  0];
T=kron(eye(60),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")

%%
%%
%%
n1=15;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;
%%
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(4*n1+1,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+2,:),"Color",'red','LineWidth',2);
%%

%%
clc;
clear;

g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');


T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)];
% g.wpos=g.wpos

n1=15;
n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
% labels={'R','Y','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=-0.00;
nk=21;

stepmax=400;
minstepmax=11;
step=1;
knum=21;
u=(4*n1+2)/nbands;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
U=1;

xinitial=zeros(nbands,nbands);
nsite=nbands/2;
for i = 1:nsite
    nup=3; 
    ndn=3;
    sx=0.0;
    sy=0.0;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy)/8;-(sx+1j*sy)/8,nup].*U;
end
%save("tmp.dat","tmp")
 load("tmp.mat")
 xinitial=tmp;

[xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);

%%
gs.onsite_modify(xinitial*10);
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

hold on;
plot(kpath,Energy(4*n1+1,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+2,:),"Color",'red','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);


knum=81;
band1=1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);

%%
function [xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
        end


        tap=sum(abs(xnew-xinitial),"all");

    
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end
        xinitial=xnew*0.8+0.2*xinitial;
        step=step+1;
    end
end

function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
    sitenum=size(Unk,2)/2;
    paulix=[0,1;1,0];pauliy=[0,-1i;1i,0];pauliz=[1,0;0,-1];
    parfor i=1:size(kindex,1)
        nki=abs(Unk(:,bandindex(i),kindex(i))).^2+nki;
        psik=reshape(Unk(:,bandindex(i),kindex(i)),[2,sitenum]);
        sxki=psik'*paulix*psik./2+sxki;
        syki=psik'*pauliy*psik./2+syki;
        szki=psik'*pauliz*psik./2+szki;
    end
    knum=size(Unk,3);
    ni=nki/knum;
    ni=reshape(ni,[2,size(ni,1)/2]);%first row spin up, second row spin down
    si=real([diag(sxki).';diag(syki).';diag(szki).']./knum);
end

function [Etot,kindex,bandindex]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
end

%%

function obj=add_elec(obj,Electric_field_in_evpA)
    % obj.wpos(:,3)=round(obj.wpos(:,3));
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        % obj.wpos(i,3)*Electric_field_in_evpA
    end 
end


function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(1,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,2),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=pi/2+0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));      
 end
 gs.iniham=gs.ham+0;  
end

function []=spin_plot(obj,si)
figure('Color','white')
pos_cart=obj.atoms;
pos_cart=pos_cart*obj.a;
h=zeros(1,3);
% plot(obj.wpos(:,1),obj.wpos(:,2),'o')
lattice_plot(obj)
hold on
si=si.';
si=si./sqrt(2);
h(2)=quiver(pos_cart(:,1),pos_cart(:,2),real(si(:,1)),real(si(:,2)),0,'Color','#007EC9','LineWidth',2,'DisplayName','xy');
hold on
% h(3)=quiver(pos_cart(:,1),pos_cart(:,2),zeros(size(si(:,1))),real(si(:,3)),0,'Color','#FF7E81','LineWidth',2,'DisplayName','z');
% legend(h(2:3))
set(gca,'Fontsize',18,'FontName','Times New Roman','linewidth',0.8)

end

function lattice_plot(obj)
lattice=obj.a;
sublattice=obj.atoms;
figure('Color','white')
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} ) ;
drawArrow([0,lattice(1,1)],[0,lattice(1,2)],'linewidth',2,'color','k')
hold on
drawArrow([0,lattice(2,1)],[0,lattice(2,2)],'linewidth',2,'color','k')
hold on
sublattice(:,1:2)=sublattice(:,1:2)*(lattice(1:2,1:2));


for i=1:size(sublattice,1)
    if sublattice(i,end)==1

        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#4DBEEE','MarkerSize',40,'MarkerFaceColor','k')
    elseif sublattice(i,end)==2
        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#4DBEEE','MarkerSize',40,'MarkerFaceColor','k')
    else
        plot3(sublattice(i,1),sublattice(i,2),sublattice(i,3),'.','Color','#0072BD','MarkerSize',40,'MarkerFaceColor','k')
    end
    text(sublattice(i,1)+0.05,sublattice(i,2)+0.05,num2str(i))
end
axis equal
axis off


end