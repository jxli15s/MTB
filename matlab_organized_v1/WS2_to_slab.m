clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplot  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("WS2");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Create 3s Slab Ham                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roate the lattice
MillerIndices=[1,-1,0]; %for a axis
% MillerIndices=[1,1,0]; %for b axis 
% MillerIndices=[0,0,1]; %for c axis 
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
g.hopr=g.hopr2;

% Calculate the bands
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=51;
nslab=3;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level

MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")
%%
% Create the supercell along z
n1=1;
n2=1;
n3=3;
gs=MTB.ham.get_supercell_wannier_3d(g,n1,n2,n3)

% rule out the z hoppings
temp=find(ismember(gs.hopr(:,3),0,'rows'));
gs.hopr=gs.hopr(temp,:);
gs.ham=gs.ham(:,:,temp);
%%
% redefine the new lattice constants and b vectors
gs.a=[6.0957,0,0;...
      0.0000,10.7920,0.0;...
      0.0000,0.00000,100]
gs.b=(2*pi*inv(gs.a))'

gs.wpos=gs.wpos(:,[2,3,1])
gs.atoms=gs.wpos*inv(gs.a)
%%
% Calculate the bands from the real space Ham
[nbands,~,nrpts]=size(gs.ham);
labels={'N','\Gamma','N1'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]
          };% hkpoints-high symmetry k points

nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(gs.ham,gs.hopr,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"2M-WS2-fplo")
% hold on;
% plot(kpath,Energy(29,:),'Color','magenta','LineWidth',2);
% plot(kpath,Energy(28,:),"Color",'red','LineWidth',2);

%%
filename="data/WS2/3-slab/wannier90_hr_3s.dat";
MTB.write_hr(gs,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Surface states by Greenfunction           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=read_fplo("WS2");

% for i=1:size(gs.wpos,1)
%     if gs.wpos(i,1)<0
%         gs.wpos(i,:)=gs.wpos(i,:)+gs.a(1,:);
%     end
%     if gs.wpos(i,2)<0
%         gs.wpos(i,:)=gs.wpos(i,:)+gs.a(2,:);
%     end
%     if gs.wpos(i,3)<0
%         gs.wpos(i,:)=gs.wpos(i,:)+gs.a(3,:);
%     end
% end
%%
% MillerIndices=[0,0,1];
MillerIndices=[1,-1,0];
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = gs.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
efermi=0; %% set Fermi Level
%%
labels={'C','\Gamma','C'}; % labels for k
hkpoints={[0.2,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.2,0.0000000000]};% hkpoints-high symmetry k points
%%
Np=1;
omegamin=-0.3;
omegamax=0.3;
omeganum=100;
omegas=linspace(omegamin,omegamax,omeganum);
nk=101;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(gs.ham,gs.hopr2,nbands,nrpts,hkpoints,nk,Np,gs.a2,gs.b2,omegamax,omegamin,omeganum);

%%
%Plot surface states
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l)
colormap(slanCM('plasma')); %magma plasma inferno cividis inferno hot heat

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
colormap(slanCM('ice'))
shading interp
% caxis([1, 50])
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

%%
% MillerIndices=[1,-1,0]; %for a axis
% MillerIndices=[1,1,0]; %for b axis 
MillerIndices=[0,0,1]; %for c axis 
Umatrix = gs.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(gs.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")
% load("data/WS2/data/E-u/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/Energy-mus-200s-001_v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001-v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-1000s-001.mat")
 % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Gamma.mat")
  load("data/WS2/data/E-u/viper/Energy-mus-300s-1-10.mat")
   % load("data/WS2/data/E-u/viper/Energy-mus-500s-1-10.mat")
  %%
figure()
% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
mus=linspace(-0.5,0.3,301)
% mus=linspace(-0.5,0.3,201)
x=repmat(mus,50,1)
% plot(x',Energy','or',...
%     'LineWidth',1,...
%     'MarkerSize',4,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r')
plot(x',Energy','o')
hold on;
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
% ylim([-0.1 0.1])
% xlim([])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
%x=repmat(mus,30,1)
%load("Energy-mus.mat")
%plot(x',Energy','*-','Color','red')

filename="data/WS2/data/final_data/E-mus/1-10-300s.dat";
for i=1:size(Energy,1)
    outlist=[mus',Energy(i,:)'];
    writeoutput(filename,outlist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Energy-mus                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/WS2/data/1000s/Energy-mus01.mat")
% load("data/WS2/data/E-u/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/Energy-mus-200s-001_v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-300s-001-v2.mat")
% load("data/WS2/data/E-u/issac/Energy-mus-1000s-001.mat")
 load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Gamma.mat")
 Energy_G=Energy;
  % load("data/WS2/data/E-u/viper/Energy-mus-500s-001-M-201mus.mat")
% load("data/WS2/data/E-u/viper/Energy-mus-500s-001-X.mat")
load("data/WS2/data/E-u/viper/Energy-mus-500s-001-Y.mat")
Energy_Y=Energy;
figure()
% mus=linspace(6.4464,7.0464,201);%mu for E_f-mu to E_f+mu of 100 points
% x=repmat(mus,20,1)
% mus=linspace(-0.5,0.3,501)
mus=linspace(-0.5,0.3,201)
x=repmat(mus,50,1)
plot(x',Energy_G','or',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
hold on;
plot(x',Energy_Y','ob',...
    'LineWidth',1,...
    'MarkerSize',4,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
ylabel('E$_{gap}$','FontSize',24,'Interpreter','latex')
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
% ylim([-0.1 0.1])
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
%x=repmat(mus,30,1)
%load("Energy-mus.mat")
%plot(x',Energy','*-','Color','red')

% 
% filename="data/WS2/data/final_data/E-mus/Gamma.dat";
% for i=1:size(Energy_G,1)
%     outlist=[mus',Energy_G(i,:)'];
%     writeoutput(filename,outlist);
% end
% filename="data/WS2/data/final_data/E-mus/Y.dat";
% for i=1:size(Energy_Y,1)
%     outlist=[mus',Energy_Y(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Energy-mus-k_{x,y}                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky_scal1.mat")
% load("data/SrSnO/fplo/finite-ky/Energy-mus-500s-0015_ky.mat")
load("data/WS2/data/finite-ky/Energy-mus-500s-001_kx.mat")
% load("data/WS2/data/finite-ky/Energy-mus-500s-001_ky.mat")
absEnergys=abs(Energys)
minEnergys=min(absEnergys)
% minEnergys=reshape(minEnergys(:,:,:),101,101)
 minEnergys=reshape(minEnergys(:,:,:),201,201)
[KX,KY]=meshgrid(ky,mus)
figure
% scatter(KX,KY,50,log(2*minEnergys))
surface(KX,KY,log(2*minEnergys),'edgecolor','none');colorbar; shading flat;
colormap(slanCM('RdBu'))
% clim([-10,-7])
% ylim([-0.2,0.1])
% xlim([0.058,0.1334])
% ylim([-0.2,-0.15])
% xlim([0.043,0.128])
% ylim([-0.049,0.076])
ylabel('$\mu$','FontSize',20,'Interpreter','latex')
xlabel('$k_x$','FontSize',20,'Interpreter','latex')
% colormap(slanCM('heat'))
shading interp
% Energy-mus-500s-0015_ky_scal1.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Wilson loop                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% load("./data/WS2/data/wloop/wilsonloop_50s_kx.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_kx_scal.mat")
% plot(x',Energy','*-','Color','red')
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky_scal1.mat")
load("./data/WS2/data/wloop/001/wilsonloop_001_80s_ky.mat")
figure('Color','white')
% plot(kx,wx)
knum=101;
band1=1;
band2=50*44;

% kx=linspace(0,1,knum);

plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
% xticks([0,1/2,1])
% xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

% filename="data/WS2/data/wloop/001/wilsonloop_001_80s_ky.dat";
% for i=1:size(wx,2)
%     outlist=[kx',wx(:,i)];
%     writeoutput(filename,outlist);
% end



%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/WS2/fplo3/POSCAR");
        pos=textread("data/WS2/fplo3/wpos");
        g.wpos=pos;
        ham=textread("data/WS2/fplo3/mydata-p1");
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

    
