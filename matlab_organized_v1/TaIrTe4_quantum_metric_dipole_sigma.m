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


Electric_field_in_evpA=0.1;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Add zeeman and displacemant field at the same time  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clear all;
%p=parpool(8)
%delete(gcp('nocreate'))
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;

[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
% labels={'Y','\Gamma','Y','R','X','\Gamma','-X'}; % labels for k
% hkpoints={[0.0,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,-0.5,0.0],...
%           [0.5,-0.5,0.0],...
%           [0.5,-0.0,0.0],...
%           [0,0,0],...
%           [-0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;

s3=[1  0
    0  -1];

Zeeman=kron(0.001*s3,eye(4));
g.add_zeeman(Zeeman)


Electric_field_in_evpA=0.00;
g=add_elec(g,Electric_field_in_evpA);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

%%
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
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

Electric_field_in_evpA=0.0;
g=add_elec(g,Electric_field_in_evpA);
%%
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
%% Berry Curvature Dipole
nk=201;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

Electric_field_in_evpA=0.0;
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
%% plot Berry Curvature
Omega=Omega_dk(:,:,1)+Omega_dk(:,:,2)+Omega_dk(:,:,3)+Omega_dk(:,:,4);
pcolor(Kx,Ky,Omega*40000);
colormap(slanCM('RdBu'))
shading interp
% xlabel('E-E_f(eV)');
% ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
% colorbar; 
% maxval=max(Omega_62./dS,[],'all');
% clim([-maxval maxval]); 

% caxis([-100 100]);
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
%%
Omega_dk=Omega_dk_bac;
Enum=100;
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate the QMD                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
g.wpos=[];
g.wpos=g.atoms*g.a;
n1=1;
n2=1;
gs = MTB.ham.get_supercell(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
%%
Electric_field_in_evpA=0.0;
gs=add_elec(gs,Electric_field_in_evpA);

%% Add Zeeman
s3=[1  0
    0  -1];
% % for hfmf
% Zeeman=kron(eye(60),0.0008*s3);
% gs.add_zeeman(Zeeman)

Zeeman_bulk=kron(0.01*s3,eye(4));
Zeeman_cdw=kron(eye(n1),Zeeman_bulk);
gs.add_zeeman(Zeeman_cdw)
%%
[nbands,~,nrpts]=size(gs.ham);
labels={'R','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points

labels={'-X','\Gamma','X','\Gamma'}; % labels for k
hkpoints={[-0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [1.0,0.0,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','Y','R','X','\Gamma','-X'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,-0.5,0.0],...
          [0.5,-0.5,0.0],...
          [0.5,-0.0,0.0],...
          [0,0,0],...
          [-0.5,0.0,0.0]};% hkpoints-high symmetry k points

% efermi=get_ef()
% efermi= 4.5302;
%%
nk=251;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

%%
% ===== 已有的部分 =====
knum  = 500;
kxline = [-0.2,0.2];
kyline = [-0.5,0.5];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]  = MTB.ham.get_bulk_plane_bands(g, Kx,Ky,Kz);
[Hamk,Unk,Enk]  = MTB.ham.get_bulk_plane_bands_with_Ham(gs, Kx,Ky,Kz);
%
efermi=calculate_ef(Enk(:),0.5);
%
Enk=Enk-efermi;
%%
plottap=2;
nk=knum;
Enum=3000;
Emin=-0.1;
Emax=0.3;
eps=(Emax-Emin)/Enum*1;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,nk,plottap);
% Dos=sum(Dos)*(Emax-Emin)/Enum;
Dos=Dos*(Emax-Emin)/Enum;

%% Plot the Dos and TDos
Dos_new=Dos;
TDos_new=TDos/norm(cross(gs.a(1,:),gs.a(2,:)))*10^16;
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos*10^15,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
% plot(TDos_new,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
plot(Eaxis,TDos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(TDos_new,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
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
% print('data/TaIrTe4/dos/TaIrTe4-DOS','-dpng','-r300')

%%
figure('Color','white')
% TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16;
TDos_new=TDos/norm(cross(gs.a(1,:),gs.a(2,:)))*10^16;
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
% xlim([-2,2])
legend('Density')
%ylim([0,200])
% print('data/TaIrTe4/dos/TaIrTe4-Density','-dpng','-r300')


%%

% 2. dk_vecs: 对应两个方向的 k 步长, 比如沿 b1, b2
% ===== 定义 dk_list (2D 情况) =====
b1 = gs.b(1,:);    % 你的 reciprocal vectors，如果是别的字段，就替换
b2 = gs.b(2,:);
dk_vecs = [b1/knum;   % 对应 "kx 方向"的单位向量
           b2/knum];  % 对应 "ky 方向"的单位向量


% 3. 其它参数
% band_list =55:66;               % 比如 valence top + conduction bottom…
band_list =1:8;  
Emin=-0.1;
Emax=0.3;
NEF=3000;
Ef_list   = linspace(Emin,Emax,NEF);
% AreaBZ    = norm(cross(b1,b2));
AreaBZ = norm(gs.b(1,:))*norm(gs.b(2,:));
eta       = 1e-3;
weights   = [];               % 让函数内部设置均匀 AreaBZ/(Nkx*Nky)
deltaE_reg = 1e-5;
TK=30;

sigma_abc = MTB.ham.get_sigma_quantum_metric_dipole( ...
    Hamk, Unk, Enk, band_list, dk_vecs, Ef_list, ...
    AreaBZ, eta, weights, deltaE_reg,TK);

[Dk_core,gk] = MTB.ham.get_Dk_qmd_plainD_core(Hamk, Unk, Enk, band_list, dk_vecs, deltaE_reg);
%%
dk121=squeeze(Dk_core(:,:,1,2,2,6));
% dk121=sign(dk121).*log(abs(dk121)+100)
s = prctile(abs(dk121(:)), 80);   % 也可 95/99
dk121 = asinh(dk121 / s);
% dk121=dk121(1:71,31:71)
% dk121=squeeze(Dk_core(:,:,2,1,1,8));
 figure()
 surf(Kx,Ky,dk121,'EdgeColor','none')
  % surf(dk121,'EdgeColor','none')
  view(2)
  % figure()
  % hold on;
  % surf(squeeze(Enk(:,:,60)),'EdgeColor','none')
 %%
 Ef=0.0;
 eta       = 1e-3;
 band_list=55:66;
 Enk_sel = Enk(:,:,band_list);
 Nb_sel=length(band_list);
 delta_n = (1/pi)*eta ./ ((Enk_sel - Ef).^2 + eta^2); % Nkx x Nky x Nb_sel
 [Nkx,Nky]=size(Kx);
Dmap = zeros(Nkx,Nky);
for in=1:Nb_sel
    Dmap = Dmap + squeeze(Dk_core(:,:,1,2,1,in)) .* delta_n(:,:,in) / (2*pi)^2;
end

% Dmap=squeeze(Dk_core(:,:,1,2,1,6)) .* delta_n(:,:,6) / (2*pi)^2;

D_yxx = sum(Dmap,'all');

figure()
surf(Kx,Ky,Dmap,'EdgeColor','none')
colorbar
% figure()
% surf(Kx,Ky,squeeze(delta_n(:,:,6)))
% figure()
% surf(Kx,Ky,squeeze(Enk_sel(:,:,6)))
%%
% 例如 σ_{yxx}(E_F)
% load("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d/qmd/sigma_qmd_V01_E001_B1e-4_k101.mat")
% Ef_list   = linspace(-0.15,0.15,601);
a = 1; b = 2; c =1;
sigma_yxx = squeeze(sigma_abc(a,b,c,:));
figure()
plot(Ef_list,1000*sigma_yxx)
%%
% load("/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K.mat")
load('/Volumes/T9/work/tb/matlab/data/TaIrTe4_2d_tb/QMD/QMD_B_1meV_30K_band.mat')
figure()
hold on;
for x=1:2
    for y=1:2
        for z=1:2
            % plot(Ef_list,1000*squeeze(sigma_abc(2,1,1,:)))
            plot(TDos_new,1000*squeeze(sigma_abc(x,y,z,:))*pi)
            % filename="sigma_abc_"+int2str(x)+int2str(y)+int2str(z)+".dat";
            % outlist=[Ef_list.',reshape(sigma_abc(x,y,z,:),300,1)];
            % writeoutput(filename,outlist)
        end
    end
end
% plot(Ef_list,1000*sigma_yxx,'ko')
ylabel('mA/V^2')
% xlabel('Energy(eV)')
xlabel('Carrier Density n (1/cm^2)')
% xlim([-0.1 0.1])

%%
abc=[2,1,1];
a=abc(1);b=abc(2);c=abc(3);
knum=500;
band1=zeros(knum/2,4);
banddk1=zeros(knum/2,4);
bandgxx1=zeros(knum/2,4);
bandgyx1=zeros(knum/2,4);

idx=0;
for i=knum:-1:knum/2+1
    idx=idx+1;
    for j=59:62
        band1(idx,j-58)=Enk(i,i,j);
        banddk1(idx,j-58)=Dk_core(i,i,a,b,c,j-54);
        bandgxx1(idx,j-58)=gk(i,i,1,1,j-54);
        bandgxx1(idx,j-58)=gk(i,i,2,1,j-54);
    end
end

band2=zeros(knum/2,4);
banddk2=zeros(knum/2,4);
bandgxx2=zeros(knum/2,4);
bandgyx2=zeros(knum/2,4);
idx=0
for i=knum/2+1:knum
    idx=idx+1;
    for j=59:62
        band2(idx,j-58)=Enk(knum/2+1,i,j);
        banddk2(idx,j-58)=Dk_core(knum/2+1,i,a,b,c,j-54);
        bandgxx2(idx,j-58)=gk(knum/2+1,i,1,1,j-54);
        bandgyx2(idx,j-58)=gk(knum/2+1,i,2,1,j-54);
    end
end

band3=zeros(knum/2,4);
banddk3=zeros(knum/2,4);
bandgxx3=zeros(knum/2,4);
bandgyx3=zeros(knum/2,4);
idx=0
for i=knum:-1:knum/2+1
    idx=idx+1
    for j=59:62
        band3(idx,j-58)=Enk(knum/2*3+1-i,i,j);
        banddk3(idx,j-58)=Dk_core(knum/2*3+1-i,i,a,b,c,j-54);
        bandgxx3(idx,j-58)=gk(knum/2*3+1-i,i,1,1,j-54);
        bandgyx3(idx,j-58)=gk(knum/2*3+1-i,i,2,1,j-54);
    end
end

band=[band1;band2;band3];
banddk=[banddk1;banddk2;banddk3];
bandgxx=[bandgxx1;bandgxx2;bandgxx3];
bandgyx=[bandgyx1;bandgyx2;bandgyx3];
%%

figure()
kdist=kpath(1:3*knum/2);
msize=10;
Dmax=max(abs(banddk),[],'all');
hold on;
for i=1:4
% scatter(kdist,band(:,i),msize,banddk(:,i)/Dmax,'filled')
scatter(kdist,band(:,i),msize,bandgxx(:,i),'filled')
% plot(squeeze(Enk(100,:,i+54)))
end

kk=kindex;
linesize=1;
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',2)
for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[-0.1 0.1],'--k','LineWidth',linesize)
end

grid off
box on
colorbar
colormap(slanCM('RdBu'))
colormap(flipud(colormap));
shading interp


% clim([-0.5,0.5])

% colorbar
% colormap(sky)
% clim([0,5000])
%%
figure()
surf(squeeze(Enk(:,:,60)),'EdgeColor','none')

xlabel('x')
%%
figure()
hold on;
% plot(Enk(151,:,60))
% plot(Enk(:,151,60))
% in_sel=5:8
% % ---------- 插值：D(k) （只取你关心的 in_sel） ----------
% D2D = squeeze(Dk_core(:,:,a,b,c,in_sel));   % Nkx x Nky
% D_line = interp2(Kx, Ky, D2D, kx_path, ky_path, 'linear');
% D_line_abs = abs(D_line);

%%
% ---------- 画图：scatter 带颜色 ----------
figure; hold on;

msize = 12; % marker size
for ib = 1:nb_plot
    scatter(kdist, E_line(:,ib), msize, D_line_abs, 'filled');
end

% 高对称点竖线
for t = xticks(:).'
    xline(t, '--k', 'LineWidth', 0.8);
end

colormap(parula);
cb = colorbar;
cb.Label.String = sprintf('|D_{%d%d%d}(k)|', a,b,c);

if ~isempty(clim)
    caxis(clim);
end

xlabel('k-path');
ylabel('Energy (eV)');
set(gca,'XTick',xticks,'XTickLabel',klabel);
title(title_str);
box on; grid on;



%%
filename='nlh-V01-E001.dat'
outlist=[Eaxis.',bcd.'];
writeoutput(filename,outlist)

%% Here You could also load("Dip2-V00-E001.mat")
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
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
% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V04-E001-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-E_neg_001-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-E001-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-E002-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-neg-E002-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-E002-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-E003-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V04-neg-E003-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V00-E001-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/xiaojun/Dip2-V00-E003-201.mat")

load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V01-E001-201.mat")

% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V00-E001-201.mat")
% load("data/TaIrTe4_2d_tb/omegadata/new/Dip2-V03-E001-201.mat")
%%
Enum=101;
Emin=-0.3;
Emax=0.3;
nsband=1:120;
T=40;
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(gs,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
figure('Color','White')
plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-0.1,max(bcd)+0.1];
%xrange=[-1,1]
%xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%
filename='nlh-V01-E001.dat'
outlist=[Eaxis.',bcd.'];
writeoutput(filename,outlist)

%%
%%
nk=201;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,nk);
dkx=gs.b(1,:)/nk;
dky=gs.b(2,:)/nk;
dS=norm(dkx)*norm(dky);
%%

%%
figure;
Omega_62=Omega_df(:,:,62)%+Omega_df(:,:,62);
Omega_62=kron(ones(2,2),Omega_62);
% pcolor(Kx,Ky,Omega_62(150:450,150:450)./dS)
pcolor(Kx,Ky,Omega_62(101:301,101:301)./dS);
colormap(slanCM('RdBu'))
shading interp
% xlabel('E-E_f(eV)');
% ylabel('$V_g$(V/nm)','Interpreter','latex','FontSize',20);
% set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8);
% colorbar; 
maxval=max(Omega_62./dS,[],'all');
clim([-maxval maxval]); 
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
% 
filename='BCD-V00-E003.dat'
for i=1:size(Ky,1)
    outlist=[Kx(:,i),Ky(:,i),Omega_62(101:301,100+i)./dS];
    writeoutput(filename,outlist)
end
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


function ef=get_ef(g,pars,model,nbands,knum)
    kxline=[-0.1,0.1];
    kyline=[-0.1,0.1];
    u=0.5;
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_kp(pars,model,nbands,Kx,Ky,Kz);
    ef=calculate_ef(Enk(:), u);
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end


function [kx_path, ky_path, kdist, xticks] = build_kpath(kpts, npts_seg)
% kpts: Np x 2
Np = size(kpts,1);
kx_path = [];
ky_path = [];
kdist = [];
xticks = zeros(Np,1);

s = 0;
for ip = 1:(Np-1)
    p0 = kpts(ip,:);
    p1 = kpts(ip+1,:);
    t = linspace(0,1,npts_seg);
    kx_seg = p0(1) + (p1(1)-p0(1))*t;
    ky_seg = p0(2) + (p1(2)-p0(2))*t;

    % 距离（如果你有真实 reciprocal vectors，可改成用笛卡尔距离）
    ds = sqrt((kx_seg - kx_seg(1)).^2 + (ky_seg - ky_seg(1)).^2);
    kdist_seg = s + ds;

    if ip > 1
        % 去掉段首重复点
        kx_seg = kx_seg(2:end);
        ky_seg = ky_seg(2:end);
        kdist_seg = kdist_seg(2:end);
    end

    kx_path = [kx_path; kx_seg(:)];
    ky_path = [ky_path; ky_seg(:)];
    kdist   = [kdist;   kdist_seg(:)];

    s = kdist(end);
    xticks(ip+1) = s;
end
end