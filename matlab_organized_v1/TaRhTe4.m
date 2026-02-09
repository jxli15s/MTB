clc;
clear;
parpool('local',4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplot  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("TaRhTe4-fplo-1s");

for i=1:size(g.wpos,1)
    if g.wpos(i,1)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(1,:);
    end
    if g.wpos(i,2)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(2,:);
    end
    if g.wpos(i,3)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(3,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate the band structure           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','X','\Gamma','Y','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0;
nk=101;

Electric_field_in_evpA=0.0*0.529177;


[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-fplo-1s",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(75,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(76,:)-efermi,"Color",'red','LineWidth',2);

filename='./data/TaRhTe4/data/TaRhTe4-bulk-band.dat';
for iband=1:nbands;
    outlist=[kpath',Energy(iband,:)'];
    writeoutput(filename,outlist);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   Calculate DOS                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nk=51;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
%%
%% save("./data/TaRhTe4/data/bulk-plane-201-GammaCenter.mat","Enk","-v7.3")
%%
% load("./data/TaRhTe4/data/bulk-plane-201-GammaCenter.mat")
% nk=201;
load("./data/TaRhTe4/data/TaIrTe4_bulk-400x400.mat")
nk=401;
plottap=1;
Enum=1501;
Emin=-0.2;
Emax=0.1;
eps=0.001;
efermi=0;

[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);

% filename='./data/TaRhTe4/data/TaRhTe4-bulk-dos-201.dat';
% outlist=[Eaxis',Dos',TDos'];
% writeoutput(filename,outlist);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Plot the 3D bands for VHS               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("./data/TaRhTe4/data/bulk-plane-201-GammaCenter.mat")
load("./data/TaRhTe4/data/TaIrTe4_bulk-400x400.mat")

nk=400;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);

figure('Color','white')
surf(Kx,Ky,Enk(:,:,76))
hold on;
surf(Kx,Ky,Enk(:,:,77))
colorbar;
% shading flat;
shading interp;
lighting gouraud;
 light('Position', [0 0 1], 'Style','infinite');
 material shiny;
colormap(slanCM('RdBu'))
hold on;

numContours = 1000;  % 
contourLevels = linspace(-0.3, 0.3, numContours);  %

contour3(Kx,Ky,Enk(:,:,77),contourLevels,'LineColor',[0.5 0.5 0.5]); 
contour3(Kx,Ky,Enk(:,:,76), contourLevels, 'LineColor',[0.5 0.5 0.5]); 
% contour3(Kx,Ky,Enk(:,:,76), 30, 'LineColor',[0.5 0.5 0.5]); 
zlim([-0.3 0.3]);
ylim([-0.2 0.2]);
caxis([-0.4 0.3])
axis 
view(-92,18)
% valence position 0.0 0.0506 0.0  -0.1407
% conduc 0.056 0.0177 0.0  0.01152
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate the Willson loop of bulk     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
knum=101;
band1=1;
band2=76;
[wx,~]=MTB.ham.get_wilsonloop(g,knum,band1,band2);

kx=linspace(0,1,knum);

% filename='./data/TaRhTe4/data/TaRhTe4-bulk-wilsonloop.dat';
% for iband=band1:band2
%     outlist=[kx',wx(:,iband)];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     Surface states                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=read_fplo("TaRhTe4-fplo-1s");

for i=1:size(gs.wpos,1)
    if gs.wpos(i,1)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(1,:)
    end
    if gs.wpos(i,2)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(2,:)
    end
    if gs.wpos(i,3)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(3,:)
    end
end
% MillerIndices=[1,0,0];
MillerIndices=[0,-1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
efermi=0; %% set Fermi Level
labels={'Y','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
% labels={'X','\Gamma','X'}; % labels for k
% hkpoints={[0.5,0.0000000000],...
%           [0.0000000000,0.0000000000],...
%           [0.5,0.0000000000]};% hkpoints-high symmetry k points


%%
nk=31;
nslab=20;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaRhTe4-slab")

% save("./data/TaRhTe4/data/ene.mat","Energy","kpath")
% filename="./data/TaRhTe4/data/TaRhTe4-slab-50s-010.dat";
% for iband=1:nbands*nslab
%     outlist=[kpath',Energy(iband,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Surface states by Greenfunction           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=read_fplo("TaRhTe4-fplo-1s");
for i=1:size(gs.wpos,1)
    if gs.wpos(i,1)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(1,:);
    end
    if gs.wpos(i,2)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(2,:);
    end
    if gs.wpos(i,3)<0
        gs.wpos(i,:)=gs.wpos(i,:)+gs.a(3,:);
    end
end
% MillerIndices=[1,0,0];
MillerIndices=[1,0,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
efermi=0; %% set Fermi Level
labels={'Y','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
%%
% labels={'X','\Gamma','X'}; % labels for k
% hkpoints={[0.0,0.0000000000],...
%           [0.5000000000,0.0000000000],...
%           [0.0,0.0000000000]};% hkpoints-high symmetry k points

%%
Np=1;
omegamin=-0.5;
omegamax=0.5;
omeganum=100;
omegas=linspace(omegamin,omegamax,omeganum);
nk=51;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);

%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l);
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_r)
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_bulk)
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)


%% For the next paper %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                       CDW states                                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=1;
n2=7;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','X','\Gamma','Y','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=21;
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%%
% filename='./data/TaRhTe4/data/TaRhTe4-17cdw-0V-band.dat';
% for iband=1:nbands;
%    outlist=[kpath',Energy(iband,:)'];
%    writeoutput(filename,outlist);
% end
%%
hold on;
plot(kpath,Energy(76*n2,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(76*n2-1,:)-efermi,"Color",'red','LineWidth',2);
%% Calculate the Wilson loop
knum=21;
band1=1;
band2=76*n2+2;
[wx,~]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);

kx=linspace(0,1,knum);
%%
n1=1;
n2=17;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.1;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','X','\Gamma','Y','R'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=21;
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])
%
%filename='./data/TaRhTe4/data/TaRhTe4-17cdw-0V-band.dat';
%for iband=1:nbands;
%    outlist=[kpath',Energy(iband,:)'];
%    writeoutput(filename,outlist);
%end






%%             %%%%%%%%%%%%%%%%%%%%%%%%
%%             %% For the next paper %%
%%             %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Nonlinear Hall Effect Berry Curvature Dipole along ky              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear Hall Effect Berry Curvature Dipole along ky
nk=21;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
%Kx_d=Kx_d+10^-3;

dky=g.b(2,:)./nk;
dky=norm(dky);
Ky_d=Ky_d+dky;

Electric_field_in_evpA=0.05*0.529177;
nsband=1:82;
tic;
fprintf("Processing on plane Enk and Unk\n");
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
toc;

tic;
fprintf("Processing on Omega_dk\n");
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(g,Kx,Ky,Kz,Enk,Unk,nsband,0.003,0);
toc;
%% Here You could also load("Dip2-V04-E01.mat")
%%
% save("Dip2-V02-E001-201.mat","Omega_dk","Enk","Enk_d","-v7.3");
%%
efermi=-0.05
Enum=1000;
Emin=-0.3;
Emax=0.3;
nsband=1:82;
T=30;
%[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd2(Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd_ky(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);
%%
figure('Color','White')
plot(Eaxis-efermi,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-0.1,max(bcd)+0.1]
xrange=[-0.2,0.2]
xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


%% Nonlinear Hall Effect Berry Curvature Dipole along kx
nk=51;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[Kx_d,Ky_d,Kz_d] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
%Kx_d=Kx_d+10^-3;

dkx=g.b(1,:)./nk;
dkx=norm(dkx);
Kx_d=Kx_d+dkx;

Electric_field_in_evpA=0.01*0.529177;
nsband=1:82;
tic;
fprintf("Processing on plane Enk and Unk\n");
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
[~,Enk_d]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx_d,Ky_d,Kz_d);
toc;

tic;
fprintf("Processing on Omega_dk\n");
[Omega_dk,KX,KY,KZ]=MTB.ham.get_Berrycurvature_dip2(g,Kx,Ky,Kz,Enk,Unk,nsband,0.003,0);
toc;


efermi=-0.05;
Enum=1000;
Emin=-0.3;
Emax=0.3;
nsband=1:82;
T=30;
[Eaxis,bcd,Omega_df]=MTB.ham.get_bcd(g,Omega_dk,Enk(:,:,nsband),Enk_d(:,:,nsband),Enum,Emin,Emax,T);

%%
figure('Color','White')
plot(Eaxis-efermi,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% plot(Eaxis,bcd,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
xlabel('E-E_f(eV)')
ylabel('$D_{xz}(\AA)$','Interpreter','latex','FontSize',24)

yrange=[min(bcd)-0.1,max(bcd)+0.1];
xrange=[-0.2,0.2];
xlim(xrange)
ylim(yrange)
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)


%%
clc;
clear;
g = MTB.geometry("TaIrTe4-fplo-1s");
g = MTB.read_poscar(g,"data/TaIrTe4/fplo/1s/POSCAR");
ham=[];
pos=textread("data/TaIrTe4/fplo/1s/wpos");
g.wpos=pos;
ham=textread("data/TaIrTe4/fplo/1s/mydata-p1");
orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
orbital_index=find(orbital_index_logical);
orbital=ham(orbital_index_logical,1:2);
orbital_num=sqrt(size(orbital_index,1));

%%
a=[];
ham2=ham;
tic;
for i=1:size(orbital_index,1)-1
    fprintf("%d\n",i)
    if (orbital_index(i+1)-orbital_index(i))>1
        a=[a;i];
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
    end
end
i=size(orbital_index,1);
ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)/g.a;
toc;

%
tic;
dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
[x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                 round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                 round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
hopr=round([x(:),y(:),z(:)]);
toc;
%
tic;
for i=1:size(orbital_index,1)-1
    fprintf("%d\n",i)
    if (orbital_index(i+1)-orbital_index(i))>1
        for j=orbital_index(i)+1:orbital_index(i+1)-1
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
        end
    end
end


i=size(orbital_index,1)-1
for j=orbital_index(i+1)+1:size(ham2,1)
    index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
    hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
end

toc;
g.ham=hambac;
g.hopr=hopr;


% g.wpos(1:10,3)=g.wpos(1:10,3)+0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)-0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)+0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)-0.2013;
%%
for i=1:size(g.wpos,1)
    if g.wpos(i,1)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(1,:)
    end
    if g.wpos(i,2)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(2,:)
    end
    if g.wpos(i,3)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(3,:)
    end
end
%%
n1=1;
n2=15;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.3;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
labels={'R','Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.0;
nk=21;
%%
Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-0.4,0.4])

%%
hold on;
plot(kpath,Energy(1140,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(1139,:),"Color",'red','LineWidth',2);
%%
knum=21;
band1=1;
band2=76*n2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(2,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,2),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=pi/2+0.01;
       V=Vamp.*(cos(2*pi/7/a*x+phi));      
 end
end


function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/TaRhTe4/fplo/1s/POSCAR");
        pos=textread("data/TaRhTe4/fplo/1s/wpos");
        g.wpos=pos;
        ham=textread("data/TaRhTe4/fplo/1s/mydata-p1");
        orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
        orbital_index=find(orbital_index_logical);
        orbital=ham(orbital_index_logical,1:2);
        orbital_num=sqrt(size(orbital_index,1));

        ham2=ham;
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
            end
        end
        i=size(orbital_index,1);
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)/g.a;
        toc;

                %
        tic;
        dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
        hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
        [x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                         round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                         round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
        hopr=round([x(:),y(:),z(:)]);
        toc;
        %
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                for j=orbital_index(i)+1:orbital_index(i+1)-1
                    index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
                    hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
                end
            end
        end


        i=size(orbital_index,1)-1
        for j=orbital_index(i+1)+1:size(ham2,1)
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
        end

        toc;
        g.ham=hambac;
        g.hopr=hopr;
        fplo=g;
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


