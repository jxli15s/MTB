clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Write the wannier90_hr.dat from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename="data/SrSnO/data_all/wtool/wannier90_hr.dat";
% write_hr(g,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.5,0.5],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;

% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(6,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(7,:),"Color",'red','LineWidth',2);

% filename='./data/SrSnO/data_all/SrSnO_bulk_band.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=31;
nslab=9;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level

MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")

% filename='./data/SrSnO/data_all/SrSnO_slab_band_001.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("SrSnO");
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
omegamin=-1;
omegamax=1;
omeganum=500;
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
[Kx,Ky]=meshgrid(kpath,omegas)
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
[Kx,Ky]=meshgrid(kpath,omegas)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the Mirrorx symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("SrSnO");
[nbands,~,nrpts]=size(g.ham);
%%
kpoint1=[0.0,0.2,0.4]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.0,0.2,0.4]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);

kpoint=kpoint1*g.b;
m1=[1,0,0,0,0;...
    0,-1,0,0,0;...
    0,0,-1,0,0;...
    0,0,0,1,0;...
    0,0,0,0,-1];
% m1=eye(5)
m1=kron(m1,1j*[0,1;1,0]);
m2=m1*1;
m3=m1*1;
m4=[-1,0,0;...
    0,-1,0;...
    0,0,1;];
% m4=eye(3)
m4=kron(m4,1j*[0,1;1,0]);
% m=[m1,zeros(size(m1)),zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
%    zeros(size(m1)),m2,zeros(size(m1)),zeros(size(m1,1),size(m4,2));...
%    zeros(size(m1)),zeros(size(m1)),m3,zeros(size(m1,1),size(m4,2));...
%    zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),zeros(size(m4,1),size(m1,2)),m4];
m=blkdiag(m1,m2,m3,m4);
% 
for i =1:size(g.wpos,1)
    m(i,:)=exp(-2j*(g.wpos(i,1))*kpoint(1)).*m(i,:);
     % m(i,:)=exp(-2j*(g.wpos(i,:))*kpoint').*m(i,:);
end

c=m*hk1*inv(m)-hk2;
max(c,[],'all')


[v,e]=eig(m);

d1=inv(v)*hk1*v;

d2=inv(v)*m*v.*1j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=1001;
band1=1;
band2=3;
[wx,unk]=MTB.ham.get_wilsonloop_mirror(g,v,knum,band1,band2);
%%
% kx=linspace(-0.17*pi,0.17*pi,knum);
% kx=linspace(-pi,pi,knum);
% filename='./data/SrSnO/data_all/SrSnO_bulk_wilsonloop_mx_all.dat';
% outlist=[kx',wx];
% writeoutput(filename,outlist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=401;
band1=1;
band2=6;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check the Mirrorxy symmetry          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint1=[0.2,0.2,0.4]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
% [Energy,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b)
kpoint2=[0.2,0.2,0.4]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
% [Energy,Psik,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b)

kpoint=kpoint1*g.b;


m12=[0,0,0,0,0,-1,0,0,0,0;...
     0,0,0,0,0,0,0,0,-1,0;...
     0,0,0,0,0,0,0,-1,0,0;...
     0,0,0,0,0,0,-1,0,0,0;...
     0,0,0,0,0,0,0,0,0,1;...
     -1,0,0,0,0,0,0,0,0,0;...
     0,0,0,-1,0,0,0,0,0,0;...
     0,0,-1,0,0,0,0,0,0,0;...
     0,-1,0,0,0,0,0,0,0,0;...
     0,0,0,0,1,0,0,0,0,0];
m12=kron(m12,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);

m3=[-1,0,0,0,0;...
    0,0,0,-1,0;...
    0,0,-1,0,0;...
    0,-1,0,0,0;...
    0,0,0,0,1];
m3=kron(m3,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0]);
% m3(1:10,:)=m3(1:10,:).*exp(1j*(g.wpos(21,:)-[g.wpos(21,2),g.wpos(21,1),g.wpos(21,3)])*kpoint');

m4=[0,0,-1;...
    0,-1,0;...
    -1,0,0;];
m4=kron(m4,1j*[0,exp(-1j*pi*3/4);exp(1j*pi*3/4),0])%.*exp(2j*(g.wpos(32,:))*kpoint');
m4(1:6,:)=m4(1:6,:).*exp(1j*(g.wpos(31,:)-[g.wpos(31,2),g.wpos(31,1),g.wpos(31,3)])*kpoint');

m=blkdiag(m12,m3,m4);


c=inv(m)*hk1*m-hk2;
max(c,[],'all')

%%
[v,e]=eig(m);

d1=inv(v)*hk1*v;

d2=inv(v)*m*v.*1j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=1001;
band1=1;
band2=3;
[wx,unk]=MTB.ham.get_wilsonloop_mirror_xy(g,v,knum,band1,band2);
%%
% kx=linspace(-0.08*pi,0.08*pi,knum);
% kx=linspace(-0.5*pi,0.5*pi,knum);
% filename='./data/SrSnO/data_all/SrSnO_bulk_wilsonloop_mxy_all_scal008.dat';
% outlist=[kx',wx];
% writeoutput(filename,outlist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Get wilsonloop of bulk states        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

knum=301;
band1=1;
band2=6;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Plot the mus-vs-ky--gap           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
% % % absEnergys=abs(Energys);
% % % minEnergys=min(absEnergys);
% % % minEnergys=reshape(minEnergys(:,:,:),100,100)
% % % [KX,KY]=meshgrid(ky,mus)
% % % figure
% % % % scatter(KX,KY,50,log(2*minEnergys))
% % % surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
% % % colormap(slanCM('inferno'))
% % % clim([-10,-7])
% % % % shading interp
% % % f=find(ismember(mus,-min(abs(mus))));
% % % figure('Color','White')
% % % E=Energys(:,63,:);
% % % for i=1:length(Energy(:,1))
% % %     plot(ky,E(i,:)-0,'Color','black','LineWidth',2);
% % %     hold on
% % % end

%% all
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky_scal1.mat")
load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
absEnergys=abs(Energys)
minEnergys=min(absEnergys)
% minEnergys=reshape(minEnergys(:,:,:),101,101)
 minEnergys=reshape(minEnergys(:,:,:),100,100)
[KX,KY]=meshgrid(ky,mus)
figure
% scatter(KX,KY,50,log(2*minEnergys))
surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
% colormap(slanCM('inferno'))
colormap(slanCM('RdBu'))
clim([-11,0])
ylim([-0.2,0.1])
% xlim([0.058,0.1334])
% ylim([-0.2,-0.15])
% xlim([0.043,0.128])
% ylim([-0.049,0.076])
ylabel('$\mu$','FontSize',20,'Interpreter','latex')
xlabel('$k_y$','FontSize',20,'Interpreter','latex')
% colormap(slanCM('heat'))
shading interp
% Energy-mus-500s-0015_ky_scal1.mat
%% ky_scal1
load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky_scal1.mat")
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
absEnergys=abs(Energys)
minEnergys=min(absEnergys)
minEnergys=reshape(minEnergys(:,:,:),101,101)
 % minEnergys=reshape(minEnergys(:,:,:),100,100)
[KX,KY]=meshgrid(ky,mus)
figure
% scatter(KX,KY,50,log(2*minEnergys))
surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
% colormap(slanCM('inferno'))
colormap(slanCM('RdBu'))
clim([-11,0])
% clim([-10,-7])
% ylim([-0.2,0.1])
% xlim([0.058,0.1334])
% ylim([-0.2,-0.15])
xlim([0.043,0.128])
ylim([-0.049,0.076])
ylabel('$\mu$','FontSize',20,'Interpreter','latex')
xlabel('$k_y$','FontSize',20,'Interpreter','latex')
% colormap(slanCM('heat'))
shading interp
% Energy-mus-500s-0015_ky_scal1.mat

%% ky_scal2
load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky_scal2.mat")
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
absEnergys=abs(Energys)
minEnergys=min(absEnergys)
minEnergys=reshape(minEnergys(:,:,:),101,101)
 % minEnergys=reshape(minEnergys(:,:,:),100,100)
[KX,KY]=meshgrid(ky,mus)
figure
% scatter(KX,KY,50,log(2*minEnergys))
surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
% colormap(slanCM('inferno'))
colormap(slanCM('RdBu'))
clim([-11,0])
% clim([-10,-7])
% ylim([-0.2,0.1])
xlim([0.058,0.1334])
ylim([-0.2,-0.15])
% xlim([0.043,0.128])
% ylim([-0.049,0.076])
ylabel('$\mu$','FontSize',20,'Interpreter','latex')
xlabel('$k_y$','FontSize',20,'Interpreter','latex')
% colormap(slanCM('heat'))
shading interp
% Energy-mus-500s-0015_ky_scal1.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Plot the E-mu diagram           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load("data/SrSnO/1000s/isaac/001/Energy-mus-1000s-001.mat")
load("data/SrSnO/1000s/isaac/001/Energy-mus-500s-001.mat")
%%
figure()
mus=linspace(-0.5,0.5,801); 
x=repmat(mus,50,1);
plot(x',Energy','.','Color','#007EC9')
ylim([-0.0104,0.0104])
ylabel('Energy(eV)')
xlabel('mu')
title('mu=0.01 layers=1000')
%%
filename="data/SrSnO/data_all/E-mu/Energy-mus-1000s-001.dat";
for i=1:size(Energy,1)
    outlist=[mus',Energy(i,:)'];
    writeoutput(filename,outlist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Plot the Wilsonloop for BdG           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kx=linspace(0,0.12,k_mesh_dim(1));
load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx.mat")
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del.mat")
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del_s1.mat")
% load("data/SrSnO/wilsonloop/60s/wilsonloop_60s_kx_del_s2.mat")
figure('Color','white')
% kx=linspace(0,1,101);
% kx=linspace(0,0.12,201);
plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
% xticks([0,1/2,1])
% xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

%%
filename="data/SrSnO/data_all/wilsonloop/wilsonloop_001_60s_kx_del.dat";
for i=1:size(wx,2)
    outlist=[kx',wx(:,i)];
    writeoutput(filename,outlist);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Functions                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/SrSnO/fplo/POSCAR");
        pos=textread("data/SrSnO/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/SrSnO/fplo/mydata-p1");
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


 function write_hr(g,filename)
    % open a file to write wannier hoppings
    fileID=fopen(filename,'w');
    
    [numBands,~,numRpts]=size(g.ham);

    % Get the current date and time
    currentTime = datetime('now');
    formattedTime = datestr(currentTime, 'mm/dd/yyyy at HH:MM:SS');

    % Write the header with the current data and time to the file
    fprintf(fileID,' write on %s\n', formattedTime);
    % Write the orbital numbers to the file
    fprintf(fileID,'\t %d\n', numBands);
    % Write the sites numbers to the file
    fprintf(fileID,'\t %d\n', numRpts);
    
    % numbers per line
    numsPerLine = 15;
    
    % get the site numbers of hoppings
    len=size(g.hopr,1);

    degeneracy = ones(1,len);

    % Write the data into the file with 15 numbers per line
    for i = 1:numsPerLine:len
        % Determine the index range of data for the current line
        endIdx = min(i+numsPerLine-1,len);
        fprintf(fileID, '%5d', degeneracy(i:endIdx));
        fprintf(fileID, '\n');
    end

    for i = 1:numRpts
        for j = 1:numBands
            for k = 1:numBands
                fprintf(fileID, '%5d %5d %5d %5d %5d %12.6f %12.6f\n',g.hopr(i,:),j,k,real(g.ham(j,k,i)),imag(g.ham(j,k,i)));
            end
        end
    end

    
 end
