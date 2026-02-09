clc;
clear;
%%
g=read_fplo("WS2");
%%
g = MTB.geometry("WS2");
g = MTB.read_poscar(g,"data/WS2/fplo3/POSCAR");
ham=[];
% fid=fopen('data/SrSnO/fplo/mydata-p1');
% tline=fgetl(fid);
% while ischar(tline)
%     disp(tline);
%     tline=fgetl(fid);
% end
pos=textread("data/WS2/fplo3/wpos");
g.wpos=pos;
ham=textread("data/WS2/fplo3/mydata-p1");
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
%%
%% Set K-path
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
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"2M-WS2-fplo")
hold on;
plot(kpath,Energy(6,:),'Color','magenta','LineWidth',2);
plot(kpath,Energy(7,:),"Color",'red','LineWidth',2);
%% Check Time Reversal Symmetry
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(22),i*s2);

h1=inv(T)*conj(hk1)*T-hk2;
max(h1,[],'all')
%% Check Inversion Symmetry
s1=[0  1
    1  0];
kpoint=[0.5,0.0,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.5,-0.0,-0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
P=blkdiag(d1,p1,p2);
pos=eye(size(g.wpos,1));
kpoint=kpoint*g.b;
for i=1:size(g.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    fprintf('%.6f',pos(i,i))
end
P=P*pos;
c=P*hk1*inv(P)-hk2;
max(c,[],'all')
a=Psik'*P*Psik;
diag(a)

%% Check Inversion Symmetry for slab
% MillerIndices=[1,-1,0];
MillerIndices=[0,0,1];
%%
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nslab=20;
s1=[0  1
    1  0];

kpoint1=[0.5,0.0]
kpoint1=kpoint1*g.b2;
hk1=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint1,g.a2)
[V,D]=eig(hk1);
[Energy1,ind]=sort(diag(D));
Psik1=V(:,ind);

kpoint2=[-0.5,0.0]
kpoint2=kpoint2*g.b2
hk2=MTB.ham.get_slab_hk(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint2,g.a2)
[V,D]=eig(hk2);
[Energy2,ind]=sort(diag(D));
Psik2=V(:,ind);

d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab))
P=kron(orbital_site,p_bulk)

c=P*hk1*inv(P)-hk2;
max(c,[],'all')
%%
a=eig(Psik1'*P*Psik1);
diag(a)
%% Parity of slab BdG
[nbands,~,nrpts]=size(g.ham);
nslab=20
mu=0.0
delta=0.03

kpoint=[0.0,0.0]
tic;
[Energy,Psik]=MTB.ham.get_slab_bands_BdG_at_q(g.ham,g.hopr2,nslab,nbands,nrpts,kpoint,g.a2,g.b2,mu,delta);
toc;
d1=kron(s1,diag(repmat([1],1,10)));
p1=kron(s1,diag(repmat([-1],1,6)));
p2=p1;
p_bulk=blkdiag(d1,p1,p2);

%p_bulk=eye(44)
orbital_site=flip(eye(nslab));
P=kron(orbital_site,p_bulk);

P=[P,zeros(size(P));zeros(size(P)),-P];

a4=eig(Psik(:,1:nslab*nbands)'*P*Psik(:,1:nslab*nbands))

%%
load("data/WS2/data/Unk/001/Unk_13201.mat")
% load("data/WS2/data/Unk/001/Unk_13200.mat")
%%
x=kron(diag(ones(1,600)),ones(1,44));
y=x*abs(a.^2);
x=1:600;
plot(x,y)
hold on;


%% Calculate the surface states

MillerIndices=[1,-1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
%%
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
%%
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=51;
nslab=20;


%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=0; %% set Fermi Level
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")


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

