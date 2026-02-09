clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("LaNiO");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','X','M','\Gamma','Z','R','A','Z'}; % labels for k
hkpoints={[0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5],...
          [0.5,0.0,0.5],...
          [0.5,0.5,0.5],...
          [0.0,0.0,0.5],...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;

% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"LaNiO-bulk")
hold on;
% plot(kpath,Energy(6,:),'Color','magenta','LineWidth',2);
% plot(kpath,Energy(7,:),"Color",'red','LineWidth',2);

% filename='./data/SrSnO/data_all/SrSnO_bulk_band.dat';
% for i=1:size(Energy,1)
%     outlist=[kpath',Energy(i,:)'];
%     writeoutput(filename,outlist);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("LaNiO");
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'X','\Gamma','M','X'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.5000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
nslab=81;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"LaNiO-slab")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Calculate Surface states             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=read_fplo("LaNiO");
MillerIndices=[1,0,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'M','\Gamma','R','A','M'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.0,0.5000000000],...
          [0.5,0.5],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
nslab=81;

[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"LaNiO-slab")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Check the Time reversal symmetry       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

g=read_fplo("LaNiO");
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(23),i*s2);
h1=T*conj(hk1)*inv(T)-hk2;
max(h1,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the 2pi periodictivity        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kpoint=[0.0,0.0,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[1.0,-0.0,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

c=hk1-hk2;
max(c,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check the Inversion symmetry at single point     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpoint1=[0.1,0.4,0.2]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=[-0.1,-0.4,-0.2]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint2,g.a,g.b);
P=eye(46);
kpoint=kpoint1*g.b;
for i =1:size(g.wpos)
    P(i,i)=exp(-2j*(g.wpos(i,:))*kpoint');
end
P(1:14,1:14)=-P(1:14,1:14);
P(35:46,35:46)=-P(35:46,35:46);

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
% b=eig(Psik(:,1:6)'*P*Psik(:,1:6))
% sum(b(1:6))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          BdG Domain Wall E-mu diagram             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-100.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-001.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-1000s-delta001-001.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-001.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-100.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-001_X.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-001_M.mat")
load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-100-A.mat")
% load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-100-M.mat")
 % load("data/LaNiO/112/E-mu/Energy-mus-500s-delta001-100-R.mat")
figure()
% mus=linspace(3.3616,3.9616,501);
x=repmat(mus,50,1)
plot(x',Energy','o')



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
        g = MTB.read_poscar(g,"data/LaNiO/112/POSCAR");
        pos=textread("data/LaNiO/112/wpos");
        g.wpos=pos;
        ham=textread("data/LaNiO/112/mydata-p1");
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
