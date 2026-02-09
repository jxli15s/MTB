clear;
clear all;
%parpool(8)
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4/POSCAR-TaIrTe4");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/wannier90_hr_p1.dat','data/TaIrTe4/wannier90_hr_p2.dat');

g = MTB.read_poscar(g,"data/TaIrTe4/terminal2/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/terminal2/wannier90_hr_p1.dat','data/TaIrTe4/terminal2/wannier90_hr_p2.dat');

%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
nk=201;
efermi=7.3742;
efermi=-0.4805; 
%% plot bulk bands
tic;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
toc;
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-bulk",0)
hold on;
% plot(kpath,Energy(4*n1+1,:)-efermi,'Color','magenta','LineWidth',2);
% plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'red','LineWidth',2);
% print(bandname,'-dpng','-r300')
%% Calculate Wilson loop
knum=101;
band1=85;
band2=88;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%% Calculate Dos
nk=151;
knum=nk;
kxline = [0, 1]; kyline = [0, 1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
Enk=reshape(Enk,[knum^2,nbands]);
%%
Enum=200;
Emin=-0.15;
Emax= 0.15;
eps=(Emax-Emin)/Enum;
plottap=1;
Tem=50;
[Eaxis,Dos,TDos]=MTB.ham.get_dos_FermiDirac(Enk-efermi,Tem,Enum,Emin,Emax,nk,plottap);
% [Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);
%%
figure()
S=norm(g.a(1,:))*norm(g.a(2,:));
plot(Eaxis,Dos/S*1e16,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
% S=norm(g.a(1,:)*0.529177)*norm(g.a(2,:)*0.529177)
% plot(Eaxis,TDos/S*1e16,'Linestyle','-','Color','red','LineWidth',2)
xlim([-0.2,0.2])
ylabel('n (eV^{-1}cm^{-2})','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
%
exportgraphics(gca, 'forjian/dos.png', 'Resolution', 300)
filename='forjian/dos.dat';
outlist=[Eaxis',Dos'/S*1e16];
writeoutput(filename,outlist)
%%
%%
figure()
% plot(Eaxis,Dos/S*1e16/2,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
% hold on;
S=norm(g.a(1,:))*norm(g.a(2,:))
plot(Eaxis,TDos/S*1e16,'Linestyle','-','Color','red','LineWidth',2)
xlim([-0.2,0.2])
ylabel('n (cm^{-2})','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
%
exportgraphics(gca, 'forjian/t-dos.png', 'Resolution', 300)
filename='forjian/t-dos.dat';
outlist=[Eaxis',Dos'/S*1e16];
writeoutput(filename,outlist)

%% plot slab bands
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
labels={'Y','\Gamma','X'}; % labels for k
hkpoints={[0.0,0.5],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
%%
nk=51;
nslab=1;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-slab",0)


%% calculate slab plane bands

knum=81;
nslab=3;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky] = g.get_Slab2Dkmesh(kxline,kyline,knum);
%% 
[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);
%% Calculate the Dos 
plottap=1
nk=knum;
Enum=101;
Emin=-10;
Emax=10;
eps=0.01;
efermi=7.3742;
Nband=size(Enk,3);
[Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);
sum(Dos)*(Emax-Emin)/Enum
%% Plot the Dos and TDos
Dos_new=Dos
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_new,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('DOS')
xlim([-2,2])
%yticks([0 5 10])
yticklabels({})
legend('DOS')
print('TaIrTe4-DOS','-dpng','-r300')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%%
figure('Color','white')
TDos_new=TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16
plot(Eaxis,TDos_new,'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('n (cm^{-2})')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
xlim([-2,2])
legend('Density')
%ylim([0,200])
print('TaIrTe4-Density','-dpng','-r300')
%%
myCluster = parcluster('Processes')
delete(myCluster.Jobs)

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

