clc;
clear;

%%
g=read_fplo("WS2")
%%
MillerIndices=[0,0,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham);
nslab=20;
zshift=g.a*inv(Urot);
zvalue=zeros(nbands*nslab*2,3);
zpos=1:nslab*2;
zpos=zpos-1;
zpos=kron(zpos,ones(1,nbands));
wpos=repmat(g.wpos,nslab*2,1);
wpos(:,3)=wpos(:,3)+zshift(3,3)*zpos';
g.wpos=wpos;

%%

[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points

efermi=0.0; %% set Fermi Level

knum=11;
k_mesh_dim=[8,11]
band1=1;
band2=nslab*nbands;
mu=0;
delta=0.03;

%%
tic;
% wx=MTB.ham.get_wilsonloop_slab_BdG_v2(g,knum,nslab,band1,band2,mu,delta)
wx=MTB.ham.get_wilsonloop_slab_BdG_v3(g,k_mesh_dim,nslab,band1,band2,mu,delta)
toc;


%% Plot      wloop issac for 1-10           
% load("./data/WS2/data/wloop/wilsonloop_50s_kx.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_kx_scal.mat")
% load("./data/WS2/data/wloop/issac/new/wilsonloop_50s_kx_scal2.mat")
% plot(x',Energy','*-','Color','red')
load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky.mat")
% load("./data/WS2/data/wloop/issac/wilsonloop_50s_ky_scal1.mat")
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
% 
% filename="data/WS2/data/wloop/issac/wilsonloop_50_ky.dat";
% for i=1:size(wx,2)
%     outlist=[kx',wx(:,i)];
%     writeoutput(filename,outlist);
% end

%% plot E-u 1-10
load('./data/WS2/data/E-u/Energy-mus.mat')
%
figure()
mus=linspace(-0.3,0.3,301);
x=repmat(mus,30,1);
plot(x',Energy','o')

%% plot E-u 001
load('./data/WS2/data/E-u/Energy-mus-300s-001.mat')
%
figure()
mus=linspace(-0.3,0.3,301);
x=repmat(mus,30,1);
plot(x',Energy','o')
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
