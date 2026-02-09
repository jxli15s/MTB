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

%%
%%
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);

Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.00;
gs=add_elec(gs,Electric_field_in_evpA);
%%
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;
bandindex=1:120;


% [Energy,Omega_k,kpath,kindex]=MTB.ham.get_bulk_bands_bcd(gs,hkpoints,nk,bandindex);
[Energy,Unk,kpath,kindex]=MTB.ham.get_bulk_bands_psi_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
weight=zeros(nbands,nbands,length(kpath));
weight=abs(Unk).*2;
weight=reshape(weight,8,15,nbands,301);
weight=sum(weight);
c=reshape(weight,15,120,301)
b=reshape(c(1,61,:),1,301)
%%
hold on;
% colormap(slanCM('RdBu'))
for i=1:15
    b=reshape(c(i,61,:),1,301);
    MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
    hold on;
    scatter(kpath,Energy(62,:),[],b);
    clim([0,2])
end


%%
weight=zeros(nbands,nbands,length(kpath));
weight=abs(Unk).*2;
weight=reshape(weight,8,15,nbands,301);
weight=sum(weight,2);
c=reshape(weight,8,120,301)
b=reshape(c(1,61,:),1,301)

%%
hold on;
% colormap(slanCM('RdBu'))
for i=1:8
    b=reshape(c(i,61,:),1,301);
    MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
    hold on;
    scatter(kpath,Energy(62,:),[],b);
    clim([0,2])
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
 gs.iniham=gs.ham+0;
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,1),a,Vamp);
 end
 % % for i=1:L/8
 % %     for j=1:8
 % %     gs.ham(8*(i-1)+j,8*(i-1)+j,onsite_index)=gs.ham(8*(i-1)+j,8*(i-1)+j,onsite_index)+moire(sub(j,1),a,Vamp);
 % %     end
 % % end
 
 function V=moire(x,a,Vamp)
       phi=0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));      
 end
 % gs.iniham=gs.ham+0;  
end