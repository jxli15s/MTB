%%
clc;
clear;

g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
%[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

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

n1=1;
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

efermi=-0.00;


stepmax=1000;
minstepmax=11;
step=1;
knum=201;
u=(4*n1+0.1)/nbands;
% u=1/3;
Electric_field_in_evpA=0;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = gs.get_Bulk2Dkmesh(kxline,kyline,knum);
critial=10^-10;
U=3;
V=1;

xinitial=zeros(nbands,nbands);
nsite=nbands/8;
nelectrons1=4*n1/2+2+2;
nelectrons2=4*n1/2-2;
% nelectrons1=0;
% nelectrons2=0;
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
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U*1.0;
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
    xinitial((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+2)=[ndn,-(sx-1j*sy);-(sx+1j*sy),nup].*U*1.0;
end
%save("tmp.dat","tmp")
 % load("tmp.mat")
 % xinitial=tmp;
    
load('data.mat')
%%
[xinitial,T_energy,ni,si,efermi]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u);
%[xinitial,T_energy,ni,si,efermi]=runhartreev2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,V,u);
%[xinitial,ni,si,efermi]=runhartreev3(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u);
%[xinitial,ni,si,efermi]=runhartreev4(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u);
%%
save('data.mat', 'xinitial');
%%
xinitial=zeros(nbands,nbands);
gs.onsite_modify(xinitial);
%%

Electric_field_in_evpA=0.00*0.529177;
nk=101;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
% ylim([-0.4,0.4])

hold on;
% plot(kpath,Energy(4*n1,:)-efermi,'Color','magenta','LineWidth',2);
% plot(kpath,Energy(4*n1+1,:)-efermi,"Color",'red','LineWidth',2);
% plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:),"Color",'blue','LineWidth',2);
% Energy_ori=Energy;
% save('energyori.mat',"Energy_ori")

load('energyori.mat')
for i=1:size(Energy_ori,1)
    plot(kpath,Energy_ori(i,:)-0.1562,'Color','red','LineWidth',2);
    % hold on
end
%% 0.0931 for n=0.2
%% 0.0747 for n=0.1
%% 0.1212 for n=0.3
%% 0.1562 for n=0.4
%% 0.3908 for n =0.6
%% 0.6378 for n =1

% Energy_ori=Energy;

load('energyori.mat')
for i=1:size(Energy_ori,1)
    plot(kpath,Energy_ori(i,:)-0.678,'Color','red','LineWidth',2);
    % hold on
end
%%
knum=101;
band1=1;
band2=4*n1;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);

knum=101;
band1=4*n1+1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
function [xinitial,T_energy,ni,si,efermi]=runhartree(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,u)
    step=1;
    T_energy=[];
    nbands=size(gs.ham,1);
    nsite=nbands/2;
    xorders=zeros(size(xinitial,1),size(xinitial,2),5);

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
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
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

        % % % step_mod=mod(step-1,4)+1; 
        % % % xorders(:,:,step_mod)=xnew;
        % % % if step<5
        % % %     xinitial=xnew*0.8+0.2*xinitial;
        % % % else
        % % %     xinitial=(xinitial+xorders(:,:,1)+xorders(:,:,2)+xorders(:,:,3)+xorders(:,:,4)+xorders(:,:,5))/6;
        % % % end
        step=step+1;
    end
end

function [xinitial,T_energy,ni,si,efermi]=runhartreev2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,minstepmax,critial,U,V,u)
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
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U*1.2;
        end

        Vfock=caldensity(Unk,kindex,bandindex);
        % size(Vfock)
        for k=1:nsite/4
            Vfock(8*(k-1)+1,8*(k-1)+1)=ni(2,4*(k-1)+1+2);
            Vfock(8*(k-1)+2,8*(k-1)+2)=ni(1,4*(k-1)+1+2);
            Vfock(8*(k-1)+3,8*(k-1)+3)=ni(2,4*(k-1)+1+3);
            Vfock(8*(k-1)+4,8*(k-1)+4)=ni(1,4*(k-1)+1+3);
            Vfock(8*(k-1)+5,8*(k-1)+5)=ni(2,4*(k-1)+1);
            Vfock(8*(k-1)+6,8*(k-1)+6)=ni(1,4*(k-1)+1);            
            Vfock(8*(k-1)+7,8*(k-1)+7)=ni(2,4*(k-1)+1+1);
            Vfock(8*(k-1)+8,8*(k-1)+8)=ni(1,4*(k-1)+1+1);
        end
        % size(Vfock)
        
        xnew=xnew+Vfock.*V;

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

function [xinitial,ni,si,efermi]=runhartreev3(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u)
        nbands=size(gs.ham,1);
        xinitial = [real(xinitial(:));imag(xinitial(:))];
        % 优化选项 for fminunc
        options = optimoptions('fminunc', ...
            'Algorithm', 'quasi-newton', ... % 使用拟牛顿方法
            'Display', 'iter', ...           % 显示迭代信息
            'OptimalityTolerance', critial, ...
            'MaxIterations', stepmax);
        objective = @(xinitial) one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        [xinitial_opt,fval] = fminunc(objective,xinitial,options);


        % % swarm_size=50;
        % % lb = -1;
        % % ub = 1;
        % % num_variables=2*nbands^2;
        % % lower_bound=lb*ones(num_variables,1);
        % % upper_bound=ub*ones(num_variables,1);
        % % options = optimoptions('particleswarm', ...
        % % 'SwarmSize', swarm_size, ...        % 粒子数量
        % % 'MaxIterations', stepmax, ...% 最大迭代次数
        % % 'Display', 'iter', ...              % 显示每次迭代的信息
        % % 'PlotFcn', 'pswplotbestf');         % 绘制优化过程
        % % objective = @(xinitial) one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        % % [xinitial_opt,fval] = particleswarm(objective,num_variables,lower_bound,upper_bound,options);



        xinitial_real=reshape(xinitial_opt(1:nbands^2),nbands,nbands);
        xinitial_imag=reshape(xinitial_opt(1+nbands^2:end),nbands,nbands);
        xinitial=xinitial_real+1j*xinitial_imag;


        disp('Relaxed xinitial:');
        disp(xinitial);

        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

%%
function [xinitial,ni,si,efermi]=runhartreev4(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,stepmax,critial,U,u)
        nbands=size(gs.ham,1);
        xinitial = xinitial(:);
        objective = @(xinitial) one_step_hf_v2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial,U,u);
        xinitial = quasi_newton(objective, xinitial);

        % disp('Relaxed xinitial:');
        % disp(xinitial);
        xinitial=reshape(xinitial,nbands,nbands);

        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end
%%
function fval = one_step_hf(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial_0,U,u)
        nbands=size(gs.ham,1);
        xinitial_real=reshape(xinitial_0(1:nbands^2),nbands,nbands);
        xinitial_imag=reshape(xinitial_0(1+nbands^2:end),nbands,nbands);
        xinitial=xinitial_real+1j*xinitial_imag;
        gs.onsite_modify(xinitial);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        % T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        nsite=nbands/2;
        xnew=zeros(size(xinitial));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end       
        xnew=xnew*0.8+0.2*xinitial;
        fval = norm(xnew - xinitial, 'fro')^2;
        % fval=sum(abs(xnew-xinitial),"all");
end

function xnew = one_step_hf_v2(gs,knum,Kx,Ky,Kz,Electric_field_in_evpA,xinitial_0,U,u)
        nbands=size(gs.ham,1);
        xinitial_0=reshape(xinitial_0,nbands,nbands);
        gs.onsite_modify(xinitial_0);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        % fprintf("Unk_one")
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        Enk=reshape(Enk,[knum^2,nbands]);
        [T_e,kindex,bandindex,efermi]=Total_energy(Enk,u);
        [ni,si]=calonsite(Unk,kindex,bandindex);
        % T_energy=[T_energy,T_e+U*ni(1,:)*ni(2,:).'];
        nsite=nbands/2;
        xnew=zeros(size(xinitial_0));
        for k = 1:nsite
            xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k));-(si(1,k)+1j*si(2,k)),ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k)/2+ni(1,k)/2,-(si(1,k)-1j*si(2,k))/8;-(si(1,k)+1j*si(2,k))/8,ni(1,k)/2+ni(2,k)/2].*U;
            % xnew((k-1)*2+1:(k-1)*2+2,(k-1)*2+1:(k-1)*2+2)=[ni(2,k),0;0,ni(1,k)].*U;
        end       
        xnew=xnew*0.8+0.2*xinitial_0;
        xnew=xnew(:);

end



%%

function [ni,si]=calonsite(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    sitenum=size(Unk,2)/2;
    sxki=zeros(sitenum,sitenum);
    syki=zeros(sitenum,sitenum);
    szki=zeros(sitenum,sitenum);
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

function Vfock=caldensity(Unk,kindex,bandindex)
    nki=zeros(size(Unk,1),1);
    Vfock=zeros(size(Unk,1),size(Unk,1));
    knum=size(Unk,3);
    for i =1:size(Unk,1)/8
        for j=1:2
            t11=zeros(size(Unk,1),size(Unk,1));
            t11(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5)=1;
            % t11_2=t11';
            t12=zeros(size(Unk,1),size(Unk,1));
            t12(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5)=1;
            % t12_2=t11';
            t21=zeros(size(Unk,1),size(Unk,1));
            t21(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5+1)=1;
            % t21_2=t11';
            t22=zeros(size(Unk,1),size(Unk,1));
            t22(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5+1)=1;
            % t22_2=t11';
            h11=0;h12=0;h21=0;h22=0;
            parfor k=1:size(kindex,1)
                unk=Unk(:,bandindex(k),kindex(k))
                h11=h11+unk'*t11*unk
                h12=h12+unk'*t12*unk
                h21=h21+unk'*t21*unk
                h22=h22+unk'*t22*unk
            end
            hh=[h11,h12;h21,h22]./knum;
            Vfock(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5)=hh(1,1);
            Vfock(8*(i-1)+2*(j-1)+1,8*(i-1)+2*(j-1)+5+1)=hh(1,2);
            Vfock(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5)=hh(2,1);
            Vfock(8*(i-1)+2*(j-1)+1+1,8*(i-1)+2*(j-1)+5+1)=hh(2,2);
        end
    end
end

function [Etot,kindex,bandindex,efermi]=Total_energy(Enk,u)
    %u: filling factor
    % tag='ele';
    [knum,~]=size(Enk);
    % if tag=="hole"
    % [a,b]=maxk(Enk(:),ceil(size(Enk(:),1)*u));
    % else
    [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    efermi = max(a,[],"all");
    % rule out the bands below the fermi level
    a=a(ceil(size(Enk(:),1)*0.5)+1:end);
    b=b(ceil(size(Enk(:),1)*0.5)+1:end);

    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;

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
     gs.ham(i,i,onsite_index)=gs.ham(i,i,onsite_index)+moire(sub(i,1),a,Vamp);
 end
 
 function V=moire(x,a,Vamp)
       phi=0;
       V=Vamp.*(cos(2*pi/15/a*x+phi));      
 end
 gs.iniham=gs.ham+0;  
end

function xn = quasi_newton(func, x0, q, tol, history)
    if nargin < 3
        q = 10;
    end
    if nargin < 4
        tol = 1e-7;
    end
    if nargin < 5
        history = {};
    end

    f0 = func(x0);
    n = length(f0);

    alpha = 1e-3;
    U = zeros(n, q);
    V = zeros(n, q);

    xn = x0;
    fn = f0;

    disp('Starting quasi-Newton...');

    % Warm-up: Build U and V
    for i = 1:q
        fn = func(xn);
        xn_1 = xn - alpha * fn;

        % Check if vector can be reshaped into a square matrix
        N = sqrt(length(xn_1));

        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        U(:, i) = fn - xn;
        V(:, i) = func(fn) - fn;
        xn = xn_1;
    end

    disp('Warmed up.');

    % Main loop
    for i = 1:10000
        fn = func(xn);

        % Compute matrix C
        C = transpose(U) * U - transpose(U) * V; % 使用普通转置

        b = transpose(U) * (xn - fn); % 使用普通转置

        % Regularization if singular
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        delta = V * (C \ b);
        xn_1 = fn - delta;

        % Check if vector can be reshaped into a square matrix
        N = sqrt(length(xn_1));

        xn_1 = reshape(xn_1, [N, N]);
        xn_1(1:N+1:end) = max(0.0, diag(xn_1));
        xn_1 = reshape(xn_1, [], 1);

        history{end+1} = xn_1;

        fn_1 = func(xn_1);
        U = [U(:, 2:end), fn_1 - xn_1];
        V = [V(:, 2:end), func(fn_1) - fn_1];

        % dx = max(abs(xn_1 - xn));
        dx = sum(abs(xn_1 - xn));
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Finished at iteration %d\n', i);
            return;
        end

        xn = xn_1;
    end

    disp('DIDN''T CONVERGE!!!!');
end

