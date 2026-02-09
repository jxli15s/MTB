clc;
clear;
g = MTB.geometry("SrSnO");
g = MTB.read_poscar(g,"data/SrSnO/fplo/POSCAR");
%%
ham=[];
% fid=fopen('data/SrSnO/fplo/mydata-p1');
% tline=fgetl(fid);
% while ischar(tline)
%     disp(tline);
%     tline=fgetl(fid);
% end
pos=textread("data/SrSnO/fplo/wpos");
g.wpos=pos;
ham=textread("data/SrSnO/fplo/mydata-p1");
orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
orbital_index=find(orbital_index_logical);
orbital=ham(orbital_index_logical,1:2);
orbital_num=sqrt(size(orbital_index,1));
%%
a=[];
ham2=ham;

for i=1:size(orbital_index,1)-1
    if (orbital_index(i+1)-orbital_index(i))>2
        a=[a;i];
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
    elseif i==size(orbital_index,1)-1
        ham2(orbital_index(i+1)+1:end,1:3)=ham2(orbital_index(i+1)+1:end,1:3)+(g.wpos(orbital(i+1,1),:)-g.wpos(orbital(i+1,2),:));
        ham2(orbital_index(i+1)+1:end,1:3)=ham2(orbital_index(i+1)+1:end,1:3)/g.a;
    end
end


dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
[x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                 round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                 round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
hopr=round([x(:),y(:),z(:)]);

for i=1:size(orbital_index,1)-1
    if (orbital_index(i+1)-orbital_index(i))>2
        for j=orbital_index(i)+1:orbital_index(i+1)-1
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
        end
    elseif i==size(orbital_index,1)-1
        for j=orbital_index(i+1)+1:size(ham2,1)
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
        end
    end
end
g.ham=hambac;
g.hopr=hopr;
%%
%% Set K-path
[nbands,~,nrpts]=size(g.ham);
%%
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

%%
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
%%
kpoint=[0.5,0.0,0.0]
[Energy,Psik,hk]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint,g.a,g.b);
P=eye(36);
kpoint=kpoint*g.b;
for i =1:size(g.wpos)
    % P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,-4.8960])*kpoint');
    % P(i,i)=exp(2j*(g.wpos(i,:)-[-4.8960,4.8960,-4.8960])*kpoint');
    P(i,i)=exp(2j*g.wpos(i,:)*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

% P=eye(36);
% P(31:end,:)=-P(31:end,:);
a=eig(Psik'*P*Psik)


%%
kpoint=[0.5,0,0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.5,0,0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
% P=eye(36);
% P(31:end,:)=-P(31:end,:);

c=inv(P)*hk1*P-hk2
max(c,[],"all")

%%
kpoint=[0.2,0.1,0.3]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.2,-0.1,-0.3]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
s2=[0  -1i
    1i  0];
T=kron(eye(18),i*s2);

h1=inv(T)*conj(hk1)*T-hk2;
max(h1,[],"all")
%%
knum=101;
band1=1;
band2=36;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
%%
[nbands,~,nrpts]=size(g.ham);
labels={'R','\Gamma','X','M','\Gamma'}; % labels for k
hkpoints={[0.5,0.5,0.5],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]...
          };% hkpoints-high symmetry k points
nk=401;
efermi=0;
mu=0.3;
delta=0.1;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_BdG(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b,mu,delta);
Egap=min(Energy(37,:))-max(Energy(36,:))
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")

%%
kpoint=[0.5,-0.1,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.5,0.1,0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
P=eye(36);
kpoint=kpoint*g.b;
for i =1:size(g.wpos)
    P(i,i)=exp(2j*(g.wpos(i,:)-[-4.8960,4.8960,-4.8960])*kpoint');
    fprintf('%.6f\n',P(i,i))
    % P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,4.8960])*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

% c=inv(P)*hk1*P-hk2;
c=P*hk1*inv(P)-hk2;
max(c,[],'all')
%%
% P=eye(36);
% P(31:end,:)=-P(31:end,:);
a=Psik'*P*Psik
diag(a)

%%
kpoint=[0.1,-0.0,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[-0.1,0.0,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint,g.a,g.b);
P=eye(36);
kpoint=kpoint*g.b;

for i =1:size(g.wpos,1)
    P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,4.8960])*kpoint');
    % P(i,i)=exp(-2j*(g.wpos(i,:)-[4.8960,4.8960,4.8960])*kpoint');
    if i>30
        P(i,i)=-P(i,i);
    end
end

c=inv(P)*hk1*P-hk2;
max(c,[],'all')
%%
% P=eye(36);
% P(31:end,:)=-P(31:end,:);
a=Psik'*P*Psik
diag(a)

