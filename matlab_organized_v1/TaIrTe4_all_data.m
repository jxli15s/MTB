clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.4;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);
nk=201;
%%
% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V04-E001-201.mat")
load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V04-E001-201.mat")
Electric_field_in_evpA=0.01;
%%
Enum=1001;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=40;
[Eaxis,bcd1,Omega_df1]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%[Eaxis,bcd1,Omega_df1]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);

% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V04-E001-201.mat")
load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V00-E001-201.mat")
Electric_field_in_evpA=0.01;

Enum=1001;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=40;
[Eaxis,bcd2,Omega_df2]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
figure('Color','White')
plot(Eaxis,bcd1,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on
plot(Eaxis,bcd2,'Linestyle','-','Color','red','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd1)-0.1,max(bcd1)+0.1]
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%filename='E-bcd-E001.dat'
%outlist=[Eaxis',bcd1',bcd2']
%writeoutput(filename,outlist)


%%
nk=201;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
dkx=gs.b(1,:)/nk;
dky=gs.b(2,:)/nk;
dS=norm(dkx)*norm(dky);


figure('Color','White')
% Omega_62=Omega_df1(:,:,62)+Omega_df1(:,:,62)
Omega_62=Omega_dk(:,:,61)+Omega_dk(:,:,61)
Omega_62=kron(ones(2,2),Omega_62)
% pcolor(Kx,Ky,Omega_62(150:450,150:450)./dS)
pcolor(Kx,Ky,Omega_62(101:301,101:301)./dS)
% pcolor(Kx,Ky,Omega_dk(:,:,61))
colormap(slanCM('RdBu'))
shading interp
% xlabel('E-E_f(eV)');
% ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
% colorbar; 
% caxis([-30 30]); 
% caxis([-8000 8000]); 
% caxis([-30000 30000]);
% caxis([-40000 40000]);
xlabel('$k_x$','Interpreter','latex','FontSize',20);
ylabel('$k_y$','Interpreter','latex','FontSize',20);
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
h=colorbar;
% h.Label.String='$\sigma_{xy} (\Omega\times cm)^-1$';
h.Label.String='$D_{xz}$';
h.Label.Interpreter='latex';
set(gca,'xtick',[])
set(gca,'ytick',[])
% set(gca,'xticklabel',[])
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])

% % % filename='BCD-V00-E001.dat'
% % % for i=1:size(Ky,1)
% % %     outlist=[Kx(:,i),Ky(:,i),Omega_62(150:450,149+i)./dS]
% % %     writeoutput(filename,outlist)
% % % end
%%

%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
%%
Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.01;
gs=add_elec(gs,Electric_field_in_evpA);

[nbands,~,nrpts]=size(gs.ham);
% % labels={'Y','\Gamma','X','Y'}; % labels for k
% % hkpoints={[0.0,0.5,0.0],...
% %           [0.0,0.0,0.0],...
% %           [0.5,0.0,0.0],...
% %           [0.0,0.5,0.0]};% hkpoints-high symmetry k points
labels={'Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=101;
bandindex=1:120;


[Energy,Omega_k,kpath,kindex]=MTB.ham.get_bulk_bands_bcd(gs,hkpoints,nk,bandindex);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
hold on;
colormap(slanCM('RdBu'))
scatter(kpath,Energy(62,:),[],-(Omega_k(:,61)+Omega_k(:,62)),"filled");
scatter(kpath,Energy(63,:),[],-(Omega_k(:,63)+Omega_k(:,64)),"filled");

scatter(kpath,Energy(59,:),[],-(Omega_k(:,59)+Omega_k(:,60)),"filled");
scatter(kpath,Energy(58,:),[],-(Omega_k(:,57)+Omega_k(:,58)),"filled");
% scatter(kpath,Energy(63,:),[],-(Omega_k(:,63)+Omega_k(:,64)),"filled");
caxis([-30 30]); 

%%
filename='Band-BC-V01-E001.dat'

for iband=1:nbands
    outlist=[kpath',Energy(iband,:)',Omega_k(:,iband)];
    writeoutput(filename,outlist)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cal High symmetry lines bands with Berry curvature breaking by onsite energy  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');

g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);

%%
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
pot=[0.02,0.04];
gs=add_site_potential(gs,pot);
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
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",0);
ylim([-0.4,0.4])
hold on;
colormap(slanCM('RdBu'))

scatter(kpath,Energy(61,:),[],(Omega_k(:,61)+Omega_k(:,62)),"filled");
scatter(kpath,Energy(63,:),[],(Omega_k(:,63)+Omega_k(:,64)),"filled");
% scatter(kpath,Energy(61,:),[],(Omega_k(:,61)+Omega_k(:,61)),"filled");
% scatter(kpath,Energy(62,:),[],(Omega_k(:,62)+Omega_k(:,62)),"filled");
% scatter(kpath,Energy(63,:),[],(Omega_k(:,63)+Omega_k(:,63)),"filled");
% scatter(kpath,Energy(64,:),[],(Omega_k(:,64)+Omega_k(:,64)),"filled");
%%
filename='Band-BC-V04-E001.dat'

for iband=1:nbands
    outlist=[kpath',Energy(iband,:)',Omega_k(:,iband)];
    writeoutput(filename,outlist)
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

function obj=add_site_potential(obj,pot)
    pot=repmat(pot,1,size(obj.ham,1)/2);
    dim_H=size(obj.ham,1);
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+pot(i);
    end 
end

function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(1,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,1),a,Vamp);
 end

 function V=moire(x,a,Vamp)
       phi=0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));
 end
end

function writeoutput(filename,list)
    file=fopen(filename,'a+');
    cloumns=size(list,2);
    raws=size(list,1);
    % fprintf(file,'Time on %s\n',datetime('today'));
    % fprintf(file,'raws %d cloums %d\n',raws,cloumns);
    for i=1:raws
        for j=1:cloumns
            fprintf(file,'%12.6f \t',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
 end