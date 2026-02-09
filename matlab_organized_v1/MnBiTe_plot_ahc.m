 clear;
 clear all;
 %p=parpool('local',8)
%delete(p)
%% Read POSCAR and Wannier
g = MTB.geometry("MnBiTe");
g = MTB.read_poscar(g,"data/MnBiTe-xu/POSCAR-MnBiTe");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/MnBiTe-xu/wannier90_hr_p1.dat","data/MnBiTe-xu/wannier90_hr_p2.dat");
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
efermi=0.0; %% set Fermi Level 0.0258

%% Load Omega_k_all and Enk_all
% save("Omega_k_all.mat","Omega_k_all","-v7.3")
% save("Enk_all.mat","Enk_all","-v7.3")
% h=6.626*10^-34;
% G=25812;%h/e^2 Om 
% coef=1/G*10^8 % sigmg(Om*Cm)^-1 sigma=e^2/h\sum(Chern_n)=coef*sum_n(Chern_n)

load("data/MnBiTe-xu/ahc/101-101/Omega_k_all-101-101.mat")
load("data/MnBiTe-xu/ahc/101-101/Enk_all_101-101.mat")
%% Calculate Hall conductivity along E
Enk=Enk_all(:,:,:,6);
Omega_k=Omega_k_all(:,:,:,6);
% Enk=flip(Enk,3); %
Omega_k=flip(Omega_k,3); %Because the sequence of Omega_k is 276:1 so it need to be fliped. But I have fixed this bug
Enum=101;
Emin=-1;
Emax=1;
T=0.001;
[Eaxis,sigma]=MTB.ham.get_ahc(Omega_k,Enk,Enum,Emin,Emax,T);
G=25812;
coef=1/G*10^8;
sigma=(sigma./coef-0.036).*coef %Because the AFM bands are not well symmetrical, so we need to set zero sigma at Ef 
%sigma(i)=sum(Omega_k(Enk<E))   
figure('Color','White')
plot(Eaxis,sigma,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$\sigma_{xy}$ ($\Omega\times cm)^{-1}$','Interpreter','latex','FontSize',24)
%yrange=[min(sigma)-20,max(sigma)+20]
%xrange=[-1,1]
%xlim(xrange)
%ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
print("data/MnBiTe-xu/ahc/101-101/Egate6_Sigma-flip.png",'-dpng','-r600')
%% Calculate Anomal Hall conductivity vs Electric field
step=11;
Electric_field_list=linspace(-0.003,0.003,step);
Enum=201;
Emin=-0.1;
Emax=0.1;
T=0.01;
G=25812;
coef=1/G*10^8;

% Eaxis=linspace(Emin,Emax,Enum);
sigma_all=zeros(step,Enum);
for i=1:step
    Enk=Enk_all(:,:,:,i);
    Omega_k=Omega_k_all(:,:,:,i);
    Omega_k=flip(Omega_k,3);
    [Eaxis,sigma]=MTB.ham.get_ahc(Omega_k,Enk,Enum,Emin,Emax,T);
    sigma=(sigma./coef-0.036).*coef; %Because the AFM bands are not well symmetrical, so we need to set zero sigma at Ef 
    sigma_all(i,:)=sigma;
end
sigma_all2=-flipud(sigma_all);
sigma_all=(sigma_all+sigma_all2)/2;

[Kx,Ky]=meshgrid(Eaxis,Electric_field_list*10000);
figure('Color','white')
[Kqx,Kqy]=meshgrid(linspace(Emin,Emax,500),linspace(-0.003,0.003,500)*10000);
sigma_all_dense=interp2(Kx,Ky,sigma_all,Kqx,Kqy);
% pcolor(Kx,Ky,sigma_all)
pcolor(Kqx,Kqy,sigma_all_dense)
colormap(slanCM('RdBu'))
shading interp
xlabel('E-E_f(eV)');
ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
h=colorbar;
% h.Label.String='$\sigma_{xy} (\Omega\times cm)^-1$';
h.Label.String='$\sigma_{xy} (\Omega)^-1$';
h.Label.Interpreter='latex';
% title(h,'$O$')

print("data/MnBiTe-xu/ahc/101-101/new-unit/EvsSigma-flip.png",'-dpng','-r600')
%% Calculate anc along E
Enk=Enk_all(:,:,:,1);
Omega_k=Omega_k_all(:,:,:,1);
% Enk=flip(Enk,3); %
Omega_k=flip(Omega_k,3); %Because the sequence of Omega_k is 276:1 so it need to be fliped. But I have fixed this bug
Enum=201;
Emin=-1;
Emax=1;
T=300;
[Eaxis,sigma]=MTB.ham.get_anc(Omega_k,Enk,Enum,Emin,Emax,T);
G=25812;
e=1.6*10^-19;
coef=1/G*10^10;
% sigma=(sigma./coef.*T-6.3*10^-5).*coef./T %Because the AFM bands are not well symmetrical, so we need to set zero sigma at Ef 
%sigma(i)=sum(Omega_k(Enk<E))   
figure('Color','White')
plot(Eaxis,sigma,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$\alpha_{xy} (A/m/K)$','Interpreter','latex','FontSize',20)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
yrange=[min(sigma)-0.2,max(sigma)+0.2]
xrange=[-1,1]
% xlim(xrange)
% ylim(yrange)

%% Calculate Anomal Nernst conductivity vs Electric field
step=11;
Electric_field_list=linspace(-0.003,0.003,step);
Enum=201;
Emin=-0.1;
Emax=0.1;
T=100;
G=25812;
coef=1/G*10^8;

% Eaxis=linspace(Emin,Emax,Enum);
sigma_all=zeros(step,Enum);
for i=1:step
    Enk=Enk_all(:,:,:,i);
    Omega_k=Omega_k_all(:,:,:,i);
    Omega_k=flip(Omega_k,3);
    [Eaxis,sigma]=MTB.ham.get_anc(Omega_k,Enk,Enum,Emin,Emax,T); 
    sigma_all(i,:)=sigma;
end
sigma_all2=-flipud(sigma_all);
sigma_all=(sigma_all+sigma_all2)/2;

[Kx,Ky]=meshgrid(Eaxis,Electric_field_list*10000);
figure('Color','white')
pcolor(Kx,Ky,sigma_all)
colormap(slanCM('RdBu'))
shading interp
xlabel('E-E_f(eV)');
ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
h=colorbar;
h.Label.String='$\alpha_{xy} (A/K)$';
h.Label.Interpreter='latex';
print("data/MnBiTe-xu/ahc/101-101/EvsAlpha-flip-100K-1meV.png",'-dpng','-r600')
% title(h,'$O$')




%% Calculate Hall Conductivity
step=11;
Electric_field_list=linspace(-0.003,0.003,step);
Enum=101;
Emin=-1;
Emax=1;
Eaxis=linspace(Emin,Emax,Enum);
sigma_all=zeros(step,Enum);
kb=8.61733*10^-5
tem=0.10
e=-5
ef=0
phase=get_phase(e,ef,kb,tem)
[Kx,Ky]=meshgrid(Eaxis,Electric_field_list)


for i=1:step
    Enk=Enk_all(:,:,:,i);
    Omega_k=Omega_k_all(:,:,:,i);
    sigma=zeros(1,Enum);
    for j=1:Enum
        E=Eaxis(j);
        % sigma(j)=sum(Omega_k(Enk<E));
        sigma(j)=get_ahc(Omega_k,Enk,E,kb,tem);
    end
    sigma_all(i,:)=sigma;
end
sigma_all=sigma_all*coef
% sigma_all_2=flipud(sigma_all)
% sigma_all(1:5,:)=-sigma_all_2(1:5,:);
figure('Color','white')
pcolor(Kx,Ky,sigma_all)
shading interp


%% Calculate Anomal Nernst  vs Electric field
step=11;
Electric_field_list=linspace(-0.003,0.003,step);
Enum=101;
Emin=-1;
Emax=1;
T=0.01;
G=25812;
coef=1/G*10^8;

% Eaxis=linspace(Emin,Emax,Enum);
sigma_all=zeros(step,Enum);
for i=1:step
    Enk=Enk_all(:,:,:,i);
    Omega_k=Omega_k_all(:,:,:,i);
    Omega_k=flip(Omega_k,3);
    [Eaxis,sigma]=MTB.ham.get_ahc(Omega_k,Enk,Enum,Emin,Emax,T);
    sigma=(sigma./coef-0.036).*coef; %Because the AFM bands are not well symmetrical, so we need to set zero sigma at Ef 
    sigma_all(i,:)=sigma;
end
sigma_all2=-flipud(sigma_all);
sigma_all=(sigma_all+sigma_all2)/2;

[Kx,Ky]=meshgrid(Eaxis,Electric_field_list*10000);
figure('Color','white')
pcolor(Kx,Ky,sigma_all)
colormap(slanCM('RdBu'))
shading interp
xlabel('E-E_f(eV)');
ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
h=colorbar;
h.Label.String='$\sigma_{xy} (\Omega\times cm)^-1$';
h.Label.Interpreter='latex';
% title(h,'$O$')

print("data/MnBiTe-xu/ahc/101-101/EvsSigma-flip.png",'-dpng','-r600')


%% 
Enk=Enk_all(:,:,:,6);
Omega_k=Omega_k_all(:,:,:,6)

Enum=101;
Emin=-1;
Emax=1;
Eaxis=linspace(Emin,Emax,Enum);
sigma=zeros(1,Enum);
for i=1:Enum
    E=Eaxis(i);
    sigma(i)=sum(Omega_k(Enk<E));
end

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

%%
kb=8.61733*10^-5
tem=0.01
e=-5
ef=0
phase=get_phase(e,ef,kb,tem)
sigma=get_ahc(Omega_k,Enk,ef,kb,tem)

Enum=301;
Emin=-10;
Emax=10;
Eaxis=linspace(Emin,Emax,Enum);
sigma2=zeros(1,Enum);
for i=1:Enum
    E=Eaxis(i);
    sigm2a(i)=sum(Omega_k(Enk<E));
end

function phase=get_phase(e,ef,kb,tem)
    f=1.0/(exp((e-ef)/kb/tem)+1.0)
    phase=((e-ef)*f+kb*tem*log(1+exp(-((e-ef)/kb/tem))))
end

% function sigma=get_ahc(Omega_k,Enk,ef,kb,tem)
%     fi=1.0/(exp((Enk-ef)/kb/tem)+1.0);
%     sigma=sum(Omega_k.*fi,'all');
% end