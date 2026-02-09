clc;
clear;
g = MTB.geometry("CaN2");
g = MTB.read_poscar(g,"data/CaN2/fplo/POSCAR");
ham=[];
% fid=fopen('data/SrSnO/fplo/mydata-p1');
% tline=fgetl(fid);
% while ischar(tline)
%     disp(tline);
%     tline=fgetl(fid);
% end
pos=textread("data/CaN2/fplo/wpos");
g.wpos=pos;
ham=textread("data/CaN2/fplo/mydata-p1");
orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
orbital_index=find(orbital_index_logical);
orbital=ham(orbital_index_logical,1:2);
orbital_num=sqrt(size(orbital_index,1));


%%
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
%%
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=51;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"2M-WS2-fplo")
hold on;
plot(kpath,Energy(9,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(10,:),"Color",'red','LineWidth',2);

%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
MillerIndices=[0,1,0];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
labels={'X','\Gamma','X'}; % labels for k
hkpoints={[-0.5,0.0],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0]};% hkpoints-high symmetry k points
%%
nk=51;
nslab=40;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"CaN2-slab",0)

%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"CaN2-slab",0)
plot(kpath,Energy(397,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(398,:),"Color",'red','LineWidth',2);
%%
nk=51;
nslab=40;
kpoint=[0,0]
hk=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint,g.a2);
hk=(hk+hk')./2.0;
[V,D]=eig(hk);
[E,ind]=sort(diag(D));
Psik=V(:,ind);
x=1:40*86
plot(x,abs(Psik(:,397).^2))
%%

a=kron(diag(ones(1,40)),ones(1,86));
c=a*abs(Psik(:,399).^2)
x=1:40
plot(x,c)
%%
clc;
clear;
g = MTB.geometry("CaN2");
g = MTB.read_poscar(g,"data/CaN2/fplo/POSCAR");
ham=[];
% fid=fopen('data/SrSnO/fplo/mydata-p1');
% tline=fgetl(fid);
% while ischar(tline)
%     disp(tline);
%     tline=fgetl(fid);
% end
pos=textread("data/CaN2/fplo/wpos");
g.wpos=pos;
ham=textread("data/CaN2/fplo/mydata-p1");
orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
orbital_index=find(orbital_index_logical);
orbital=ham(orbital_index_logical,1:2);
orbital_num=sqrt(size(orbital_index,1));


%%
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

n1=10;
n2=10;
gs = MTB.ham.get_supercell(g,n1,n2);
%%
vals=sort(real(eig(gs.ham(:,:,5))));
Energy_obs=vals;
x=1:8600;
%%
figure()
plot(x,real(Energy_obs),'*')