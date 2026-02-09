clear;
clear all;
%parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/qe/gamma_low/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
% 
% g = MTB.read_poscar(g,"data/TaIrTe4/terminal2/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/terminal2/wannier90_hr_p1.dat','data/TaIrTe4/terminal2/wannier90_hr_p2.dat');
%%
%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.000000],...
          [0.0000000000,0.5000000000,0.000000],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points

% labels={'-X','\Gamma','X'}; % labels for k
% hkpoints={[0.5000000000,0.4000000000,0.000000],...
%           [0.0,0.4,0.0],...
%           [0.5,0.4,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=1.242;
% efermi=-0.4805
%% plot bulk bands
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-bulk",0)
hold on;

plot(kpath,Energy(28,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(60,:)-efermi,"Color",'red','LineWidth',2);
plot(kpath,Energy(88,:)-efermi,"Color",'blue','LineWidth',2);
%%
filename='TaIrTe4_bulk_band.dat';
for i=1:size(Energy,1)
    outlist=[kpath',Energy(i,:)'];
    writeoutput(filename,outlist);
end

%%
%% Calculate Wilson loop
knum=201;
band1=85;
band2=88;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%% Calculate plane bands

knum=501;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
% [Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);

%%
% save("TaIrTe4_1sEnk-501x501.mat","Enk","-v7.3")
load("TaIrTe4_1sEnk-501x501.mat")
Occ=60;
efermi=1.242;
% Occ=88;
% filename='TaIrTe4_1sEnk-501x501_new.dat';
% writeEnk(Enk,Kx,Ky,Occ,filename)

%% plot the plane 
figure('Color','White')
hold on;
% V = [-1, 0, 1];
% V = [0.02, 0.04761,0.04761*2,0.14761*3,0.14761*4,0.14761*5,0.14761*6,0.14761*7,0.14761*8];
% V = 0:0.05/5:0.3;
V = 0:0.0467/2:0.3;
surf(Kx,Ky,Enk(:,:,62)-efermi);
contour3(Kx, Ky, Enk(:,:,62)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1);
% surf(Kx,Ky,Enk(:,:,89)-efermi);
% contour3(Kx, Ky, Enk(:,:,89)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
% V=V+0.005;
% surf(Kx,Ky,Enk(:,:,61)-efermi+0.005);
% contour3(Kx, Ky, Enk(:,:,61)-efermi+0.005, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
V=0:-0.063/2:-0.3;
surf(Kx,Ky,Enk(:,:,60)-efermi);
contour3(Kx, Ky, Enk(:,:,60)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1);
% surf(Kx,Ky,Enk(:,:,88)-efermi);
% contour3(Kx, Ky, Enk(:,:,88)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
% V=V-0.005;
% surf(Kx,Ky,Enk(:,:,60)-efermi-0.005);
% contour3(Kx, Ky, Enk(:,:,60)-efermi-0.005, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
% V = [-0.02, -0.062,-0.1,-0.14761,-0.14761-0.062];

% 设置真实比例
% axis equal;
% 
% colormap(slanCM('bwr'))
% colormap(slanCM('bjy'))
colormap(slanCM('RdBu'))
shading interp

% camlight    
lighting gouraud
% grid off;
% material dull
% lightangle(-45,0)
% lightangle(45,80)
% lightangle(90,45)
% colormap(coolwarm);
% 设置手动轴范围
% xlim([-15 15]); % 限制 X 轴范围
xlim([-0.3 0.3]);   % 限制 Y 轴范围
ylim([-0.2 0.2]);   % 限制 Y 轴范围
zlim([-0.3,0.3])
clim([-0.29 0.29]); 
view(-7,10)
grid on;
colorbar
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
xlabel("Kx")
ylabel("Ky")
zlabel("Energy")
print('-depsc2', '-r600', 'example_figure.eps'); % EPS (vector)
print('-dpdf', '-r600', 'example_figure.pdf');   % PDF (vector)
print('-dpng', '-r600', 'example_figure.png');   % PNG (raster)
% exportgraphics(gcf, 'figure_name.eps', 'ContentType', 'vector');
% exportgraphics(gcf, 'figure_name.pdf', 'ContentType', 'vector');
% exportgraphics(gcf, 'figure_name.png', 'Resolution', 300);
%% Calculate Dos
nk=knum;
Enum=120;
Emin=-0.25;
Emax= 0.25;
eps=(Emax-Emin)/Enum;
plottap=1;
Tem=50;
% [Eaxis,Dos,TDos]=MTB.ham.get_dos_FermiDirac(Enk-efermi,Tem,Enum,Emin,Emax,nk,plottap);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);
%%
figure()
plot(Eaxis,Dos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
plot(Eaxis,TDos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
%plot(Eaxis,TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16,'Linestyle','-','Color','#4DA1D7','LineWidth',2)

%%
%% Enk(nk,nk,nband) Kx(nk,nk) Ky(nk,nk)
filename='TaIrTe4_Enk-400x400.dat'
writeEnk(Enk,Kx,Ky,Occ,filename)


%%
filename="TaIrTe4_qe_dos_.dat";
outlist=[Eaxis',Dos'];
writeoutput(filename,outlist)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/qe/gamma_low/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
MillerIndices=[1,0,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=51;
nslab=101;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=1.241; %% set Fermi Level

MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/qe/gamma_low/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
MillerIndices=[1,0,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
labels={'X','\Gamma','X'};
Np=1;
efermi=1.241;
omegamin=-0.5-efermi;
omegamax=0.5-efermi;
omeganum=100;
omegas=linspace(omegamin,omegamax,omeganum);
nk=201;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(g.ham,g.hopr2,nbands,nrpts,hkpoints,nk,Np,g.a2,g.b2',omegamax,omegamin,omeganum);

%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('inferno'))

shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_r)
colormap(slanCM('hot'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_bulk)
colormap(slanCM('heat'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)


%%%%%%%%%%%%%%%%%
%% FPLO VERSION
%%%%%%%%%%%%%%%%%
clc;
clear;

g=read_fplo("TaIrTe4-fplo-1s");

% g.wpos(1:10,3)=g.wpos(1:10,3)+0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)-0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)+0.2013;
% g.wpos(1:10,3)=g.wpos(1:10,3)-0.2013;
%%
for i=1:size(g.wpos,1)
    if g.wpos(i,3)<0
        g.wpos(i,:)=g.wpos(i,:)+g.a(3,:);
    end
end

%%

[nbands,~,nrpts]=size(g.ham);
labels={'R','Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.05;
nk=101;

Electric_field_in_evpA=0.0*0.529177;


[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-fplo-1s",Electric_field_in_evpA*10000);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-fplo-1s",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(75,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(76,:)-efermi,"Color",'red','LineWidth',2);


%%
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

function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/TaIrTe4/fplo/1s/POSCAR");
        pos=textread("data/TaIrTe4/fplo/1s/wpos");
        g.wpos=pos;
        ham=textread("data/TaIrTe4/fplo/1s/mydata-p1");
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

%% Enk(nk,nk,nband) Kx(nk,nk) Ky(nk,nk)
function writeEnk(Enk,Kx,Ky,Occ,filename)
    file=fopen(filename, 'w');
    fprintf(file, '%s\n', '# kx     ky     kz     Ev5     Ev4     Ev3     Ev2     Ev1     Ec1     Ec2     Ec3     Ec4     Ec5');
    for i=1:size(Enk,1)
        for j=1:size(Enk,2)
            fprintf(file, [repmat('%12.6f',1,13) '\n'],...
                Kx(i,j), Ky(i,j), 0.0, ...
                Enk(i,j,Occ-4), Enk(i,j,Occ-3), Enk(i,j,Occ-2), Enk(i,j,Occ-1), Enk(i,j,Occ), ...
                Enk(i,j,Occ+1), Enk(i,j,Occ+2), Enk(i,j,Occ+3), Enk(i,j,Occ+4), Enk(i,j,Occ+5));
        end
        fprintf(file,'\n');
    end
    fclose(file);
end