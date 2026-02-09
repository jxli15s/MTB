%%
for n1=10:11 
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
T=zeros(size(g.ham,1),size(g.ham,2));
for i=1:size(g.ham,1)
    if mod(i,2)==1
    T(i,ceil(i/2))=1;
    else
    T(i,4+i/2)=1;
    end
end

for i=1:size(g.ham,3)
    g.ham(:,:,i)=T*g.ham(:,:,i)*inv(T);
end

g.wpos=[];
g.wpos=g.atoms*g.a;
% orbital_num=[4,4];
% g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
%         repmat(g.wpos(2,:),[orbital_num(2),1])
%     ]
g.wpos=[g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:);...
        g.wpos(1,:);...
        g.wpos(1,:);...
        g.wpos(2,:);...
        g.wpos(2,:)];
% g.wpos=g.wpos

n2=1;
gs = MTB.ham.get_supercell_wannier(g,n1,n2);
Vamp=0.0;
gs=moire_potential(g,gs,Vamp);
[nbands,~,nrpts]=size(gs.ham);
% labels={'R','Y','\Gamma','X','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.5,0.0,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','X','Y'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=0.00;
nk=21;

stepmax=1000;
minstepmax=11;
step=1;
knum=21;
u=(4*n1+0.3)/nbands;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
% U=[0.1,0.25,0.5,1,2];
Us=[1,2];

for uindex=1:size(Us,2)
U=Us(uindex);
xinitial=zeros(nbands,nbands);
nsite=nbands/8;
nelectrons1=4*n1/2+6+0.3;
nelectrons2=4*n1/2-6;
for i = nsite*3+1:nsite*4
    nup=nelectrons2/nsite/2/2; 
    ndn=nelectrons2/nsite/2/2;  
    sx=0.5;
    sy=0.5;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U;
end
for i = nsite*1+1:nsite*2
    nup=nelectrons1/nsite/2/2; 
    ndn=nelectrons1/nsite/2/2;  
    sx=0.5;
    sy=0.5;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U*1.2;
end
for i = nsite*2+1:nsite*3
    nup=nelectrons2/nsite/2/2; 
    ndn=nelectrons2/nsite/2/2; 
    sx=0.5;
    sy=0.5;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U;
end
for i = 1:nsite
    nup=nelectrons1/nsite/2/2; 
    ndn=nelectrons1/nsite/2/2;  
    sx=0.5;
    sy=0.5;
    % xinitial=[xinitial,nup,sx,sy,ndn];
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U*1.2;
end
%save("tmp.dat","tmp")
 % load("tmp.mat")
 % xinitial=tmp;

% [xinitial,T_energy,ni,si]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);
[xinitial,T_energy,ni,si,efermi]=runhartreev2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);
gs.onsite_modify(xinitial);

Electric_field_in_evpA=0.00*0.529177;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
bandname="TIT_band_U"+int2str(uindex)+"_n"+int2str(n1);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1+1,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'red','LineWidth',2);
print(bandname,'-dpng','-r600')
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);



knum2=101;
band1=1;
band2=4*n1;
[wx1,unk]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
wxname="TIT_wx1_U"+int2str(uindex)+"_n"+int2str(n1);
print(wxname,'-dpng','-r600')

knum2=101;
band1=4*n1+1;
band2=4*n1+2;
[wx2,unk]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
wxname="TIT_wx2_U"+int2str(uindex)+"_n"+int2str(n1);
print(wxname,'-dpng','-r600')

knum2=101;
band1=1;
band2=4*n1+2;
[wx3,unk]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
wxname="TIT_wx3_U"+int2str(uindex)+"_n"+int2str(n1);
print(wxname,'-dpng','-r600')

filename="TIT_data_U"+int2str(uindex)+"_n"+int2str(n1)+".mat";
save(filename,"xinitial","T_energy","ni","si","efermi","wx1","wx2","wx3");
end
end
%%
function [xinitial,T_energy,ni,si,efermi]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end


        tap=sum(abs(xnew-xinitial),"all");

    
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end
        xinitial=xnew*0.8+0.2*xinitial;
        step=step+1;
    end
end

function [xinitial,T_energy,ni,si,efermi]=runhartreev2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/8;
    for i=1:stepmax
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        xnew=zeros(size(xinitial));
        for k = nsite*3+1:nsite*4
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U*2;
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
        end
        for k = nsite*1+1:nsite*2
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U/2;
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U*1.2;       
        end
        for k = nsite*2+1:nsite*3
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
        end
        for k = 1:nsite
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U*1.2;
        end


        tap=sum(abs(xnew-xinitial),"all");

    
        if abs(tap)>critial && step>1
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        elseif  abs(tap)>critial && step<2
            disp("Initial Order: ")
        elseif abs(tap)<critial   && step>minstepmax
            disp("coverged ")
            break
        else
            disp("iter step "+num2str(step)+" diff order = "+ num2str(tap))
        end
        xinitial=xnew*0.8+0.2*xinitial;
        step=step+1;
    end
end

function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
    sitenum=size(Unk,2)/2;
    paulix=[0,1;1,0];pauliy=[0,-1i;1i,0];pauliz=[1,0;0,-1];
    parfor i=1:size(kindex,1)
        nki=abs(Unk(:,bandindex(i),kindex(i))).^2+nki;
        psik=reshape(Unk(:,bandindex(i),kindex(i)),[2,sitenum]);
        sxki=psik'*paulix*psik./2+sxki;
        syki=psik'*pauliy*psik./2+syki;
        szki=psik'*pauliz*psik./2+szki;
    end
    knum=size(Unk,3);
    ni=nki/knum;
    ni=reshape(ni,[2,size(ni,1)/2]);%first row spin up, second row spin down
    si=real([diag(sxki).';diag(syki).';diag(szki).']./knum);
end

function [Etot,kindex,bandindex,efermi]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
    efermi=max(a,[],"all");
end

%%

function obj=add_elec(obj,Electric_field_in_evpA)
    % obj.wpos(:,3)=round(obj.wpos(:,3));
    dim_H=size(obj.ham,1);
    hke=zeros(dim_H,dim_H);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        % obj.wpos(i,3)*Electric_field_in_evpA
    end 
end


function gs=moire_potential(g,gs,Vamp)
 a=norm(g.a(1,:));
 sub=gs.wpos;
 L=size(sub,1);
 onsite_index=find(ismember(gs.hopr,[0,0,0],'rows'));
 for i=1:L
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,2),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=pi/2+0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));      
 end
 gs.iniham=gs.ham+0;  
end