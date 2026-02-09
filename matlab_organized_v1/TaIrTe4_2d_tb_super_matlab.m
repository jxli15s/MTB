clear;
clear all;
%p=parpool(8)
%delete(gcp('nocreate'))
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
%%
g.wpos=[];
g.wpos=g.atoms*g.a;

[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;


Electric_field_in_evpA=0.01;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
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
Omega_dk_bac=Omega_dk;
%%

Omega_dk=Omega_dk_bac;
Enum=1000;
Emin=-0.3;
Emax=0.3;
T=40;
% [Eaxis,bcd,~]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)
yrange=[min(bcd)-2,max(bcd)+2]
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


Eaxis_1=Eaxis;
bcd_1=bcd;
%%
nk=301;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
dkx=g.b(1,:)/nk;
dky=g.b(2,:)/nk;
dS=norm(dkx)*norm(dky);

Omega_62=Omega_df(:,:,6)%+Omega_df(:,:,61)
Omega_62=kron(ones(2,2),Omega_62)
pcolor(Kx,Ky,Omega_62(150:450,150:450)./dS)
colormap(slanCM('RdBu'))
shading interp
% xlabel('E-E_f(eV)');
% ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
% colorbar; 
% caxis([-800 800]);
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
%% Create the supercell and add cdw potential
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.4;
gs=moire_potential(g,gs,Vamp);

%%
% %% Nonlinear Hall Effect Berry Curvature Dipole
% nk=21;
% kxline=[0,1];
% kyline=[0,1];
% [Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
% 
% Electric_field_in_evpA=0.1;
% nsband=1:120;
% tic;
% [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
% [Kx_d,Ky_d,Kz_d] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
% Kx_d=Kx_d+10^-3;
% [~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
% [Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(gs,Kx,Ky,Kz,Enk,Unk,nsband,0.001,0);
% toc;
%% Here You could also load("Dip2-V00-E001.mat")
n1=15;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);

Vamp=0.4;
gs=moire_potential(g,gs,Vamp);
nk=301;
%load("data/TaIrTe4_2d_tb/omegadata/Dip2-V04-E001-301.mat")
load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V02-E001-201.mat")
Electric_field_in_evpA=0.01;
%%
Enum=1000;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=30;
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-0.1,max(bcd)+0.1]
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%
nk=301;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
dkx=gs.b(1,:)/nk;
dky=gs.b(2,:)/nk;
dS=norm(dkx)*norm(dky);
% 
%%
Omega_62=Omega_df(:,:,62)%+Omega_df(:,:,61)
Omega_62=kron(ones(2,2),Omega_62)
% pcolor(Kx,Ky,Omega_62(150:450,150:450)./dS)
pcolor(Kx,Ky,Omega_62(150:450,150:450)./dS)
colormap(slanCM('RdBu'))
shading interp
% xlabel('E-E_f(eV)');
% ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
% colorbar; 
caxis([-8000 8000]); 
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
            fprintf(file,'%12.6f',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
 end
