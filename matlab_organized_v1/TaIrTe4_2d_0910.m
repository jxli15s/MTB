clear;
clear all;
%p=parpool(8)
 g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

% g = MTB.read_poscar(g,"data/TaIrTe4/qe/gamma_low/POSCAR-TaIrTe4");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
%%
g.wpos=[];
g.wpos=g.atoms*g.a;

% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
%
Electric_field_in_evpA=0.1;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
[nbands,~,nrpts]=size(g.ham);
% labels={'R','Y','\Gamma','X'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0]};% hkpoints-high symmetry k points
labels={'Y1','Y2','\Gamma','X'}; % labels for k
hkpoints={[0.033,-0.5,0.0],...
          [0.033,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
%
Electric_field_in_evpA=0.1;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
%%
efermi=0.0
nk=801;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
Enum=100;
Emin=-0.25;
Emax= 0.25;
eps=(Emax-Emin)/Enum;
plottap=1;
Tem=100;
% [Eaxis,Dos,TDos]=MTB.ham.get_dos_FermiDirac(Enk-efermi,Tem,Enum,Emin,Emax,nk,plottap);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);

%%
% TDos=sum(Dos)*(Emax-Emin)/Enum

%% Plot the Dos and TDos
Dos_new=Dos
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('States/eV')
% xlim([-2,2])
%yticks([0 5 10])
% ylim([0,100])
%yticklabels({})
legend('DOS')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
print('data/TaIrTe4/dos/TaIrTe4-DOS','-dpng','-r300')

%%
figure('Color','white')
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16;
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
xlim([-0.2,0.2])
legend('Density')
%ylim([0,200])
% print('data/TaIrTe4/dos/TaIrTe4-Density','-dpng','-r300')

filename='forhaotian/dos_tdos.dat'
outlist=[Eaxis',Dos',TDos_new']
writeoutput(filename,outlist)
%% Berry Curvature Dipole
nk=1001;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.01;
nsband=1:8;
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
knum=nk;
dkx=g.b(1,:)./knum;
dkx=norm(dkx);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
 Kx_d=Kx_d+dkx;
%Kx_d=Kx_d+10^-3
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);

[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(g,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);

%%
Enum=1001;
Emin=-0.3;
Emax=0.3;
T=40;
% [Eaxis,bcd,~]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

figure('Color','White')
plot(Eaxis-efermi,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)
yrange=[min(bcd)-2,max(bcd)+2]
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%
%% Calculate bands
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
%%
g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=2;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);

[nbands,~,nrpts]=size(gs.ham);
% labels={'R','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.5,0.0]};% hkpoints-high symmetry k points

labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points

efermi=0.03;
nk=101;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.2,0.2])
%%
knum=81;
band1=1;
band2=62;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);

%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');
efermi=0.03;
g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.01;
gs=add_elec(gs,Electric_field_in_evpA);
nk=201;
%%
% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V04-E001-201.mat")
% load("data/TaIrTe4_2d_tb/1019/omegadata/Dip2-V00-E001-201.mat")
% load("data/TaIrTe4_2d_tb/1019/omegadata/Dip2-V01-E001-201.mat")
% load("data/TaIrTe4_2d_tb/1019/omegadata/Dip2-V00-E001-201.mat")
load("data/TaIrTe4_2d_tb/1019/omegadata/Dip2-V01-E002-201.mat")
Electric_field_in_evpA=0.01;
%%
Enum=1001;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=40;
[Eaxis,bcd1,Omega_df1]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

%%
figure('Color','White')
plot(Eaxis-efermi,bcd1,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on
yline(0, 'k--','LineWidth',1); 
% plot(Eaxis,bcd2,'Linestyle','-','Color','red','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd1)-1,max(bcd1)+1];
xrange=[-0.2,0.2];
xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%filename='E-bcd-E001.dat'
%outlist=[Eaxis',bcd1',bcd2']
%writeoutput(filename,outlist)

%%
%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

g.wpos=[];
g.wpos=g.atoms*g.a;
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
%%
Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
Electric_field_in_evpA=0.02;
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

efermi=0.03;
nk=101;
bandindex=1:120;


[Energy,Omega_k,kpath,kindex]=MTB.ham.get_bulk_bands_bcd(gs,hkpoints,nk,bandindex);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
hold on;
colormap(slanCM('RdBu'))
scatter(kpath,Energy(62,:)-efermi,[],-(Omega_k(:,61)+Omega_k(:,62)),"filled");
scatter(kpath,Energy(63,:)-efermi,[],-(Omega_k(:,63)+Omega_k(:,64)),"filled");

scatter(kpath,Energy(59,:)-efermi,[],-(Omega_k(:,59)+Omega_k(:,60)),"filled");
scatter(kpath,Energy(58,:)-efermi,[],-(Omega_k(:,57)+Omega_k(:,58)),"filled");
% scatter(kpath,Energy(63,:),[],-(Omega_k(:,63)+Omega_k(:,64)),"filled");
caxis([-40 40]); 

%%
function obj=add_elec(obj,Electric_field_in_evpA)
    obj.wpos(:,3)=round(obj.wpos(:,3));
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        obj.wpos(i,3)*Electric_field_in_evpA
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