clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/wannier90_hr_p1.dat','data/TaIrTe4_2d/wannier90_hr_p2.dat');
[nbands,~,nrpts]=size(g.ham);
efermi=-0.4423;
nk=101;

g.wpos=[]
g.wpos=g.atoms*g.a
orbital_num=[18,18,12,12,8,8,8,8,8,8,8,8]
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
        repmat(g.wpos(3,:),[orbital_num(3),1]);...
        repmat(g.wpos(4,:),[orbital_num(4),1]);...
        repmat(g.wpos(5,:),[orbital_num(5),1]);...
        repmat(g.wpos(6,:),[orbital_num(6),1]);...
        repmat(g.wpos(7,:),[orbital_num(7),1]);...
        repmat(g.wpos(8,:),[orbital_num(8),1]);...
        repmat(g.wpos(9,:),[orbital_num(9),1]);...
        repmat(g.wpos(10,:),[orbital_num(10),1]);...
        repmat(g.wpos(11,:),[orbital_num(11),1]);...
        repmat(g.wpos(12,:),[orbital_num(12),1]);...
    ]
g.wpos=g.wpos*5

%% Berry Curvature Dipole
nk=31;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
% Electric_field_list=[0.0,0.005,0.01,0.02,0.03,0.04,0.05];
Electric_field_list=[0.0,0.05]
Enum=101;
Emin=-1;
Emax=1;
T=0.01;
nsband=1:100;
Enk_all=zeros(nk,nk,nbands,length(Electric_field_list));
Omega_dk_all=zeros(nk,nk,size(nsband,2),length(Electric_field_list));
bcd_all=zeros(length(Electric_field_list),Enum);
for i=1:length(Electric_field_list)
    tic;
	Electric_field_in_evpA=Electric_field_list(i);
	[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
	Enk_all(:,:,:,i)=Enk-efermi;
	[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip(g,Kx,Ky,Kz,Enk,Unk,nsband,0.03,0);
	Omega_dk_all(:,:,:,i)=Omega_dk;
	[Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk,Enk(:,:,nsband)-efermi,Enum,Emin,Emax,T);
	bcd_all(i,:)=bcd;
    toc;
end

% save("Omega_k_all-401.mat","Omega_dk_all","-v7.3");
% save("Enk_all-401.mat","Enk_all","-v7.3");
% save("bcd.mat","bcd","Eaxis","-v7.3");

%%
load("data/TaIrTe4_2d/51-51/bcd-51.mat");
load("data/TaIrTe4_2d/51-51/Enk_all-51.mat");
load("data/TaIrTe4_2d/51-51/Omega_k_all-51.mat");

%%
nk=51;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_list=[0.0,0.005,0.01,0.02,0.03,0.04,0.05];
Enum=501;
Emin=-1;
Emax=1;
T=200;
nsband=1:90;
bcd_all=zeros(length(Electric_field_list),Enum);
for i=1:length(Electric_field_list)
    tic;
    Omega_dk=Omega_dk_all(:,:,:,i);
    Enk=Enk_all(:,:,:,i);
    [Eaxis,bcd]=MTB.ham.get_bcd(Omega_dk(:,:,nsband),Enk(:,:,nsband),Enum,Emin,Emax,T);
    bcd_all(i,:)=bcd;
    toc;
end
%%
figure('Color','White')
% plot(Eaxis,bcd_all(:,:),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,bcd_all(:,:),'Linestyle','-','LineWidth',2)
legend('E=0.000','E=0.005','E=0.010','E=0.020','E=0.030','E=0.040','E=0.050','Location','best','NumColumns',1)
legend('boxon')
xlabel('E-E_f(eV)')
ylabel('$D_{xz} (\AA)$','Interpreter','latex','FontSize',24)
% yrange=[-6,6]
%xrange=[-1,1]
xlim([-0.2,0.2])
% ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)



