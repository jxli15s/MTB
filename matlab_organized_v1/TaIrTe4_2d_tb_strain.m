%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% load('xinitial_good.mat')
% Step 1: Initialize Geometry and Hamiltonian
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
% g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);

% Step 4: Construct Supercell
n1 = 1; n2 = 1; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

gs.iniham = gs.ham + 0; % Store the initial Ham
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

% Step 5: Neighbor Search
[result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Calculate the plane               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate plane bands
knum=501;
kxline=[-0.5,0.5];
kyline=[-0.5,0.5];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
%%
Occ=4;
filename='TaIrTe4_1sEnk-501x501.dat';
writeEnk(Enk,Kx,Ky,Occ,filename)
%%
efermi=0
figure()
hold on;
% V = [-1, 0, 1];
% V = [0.02, 0.04761,0.04761*2,0.14761*3,0.14761*4,0.14761*5,0.14761*6,0.14761*7,0.14761*8];
V = 0:0.04667/50:0.3;
surf(Kx,Ky,Enk(:,:,5)-efermi);
contour3(Kx, Ky, Enk(:,:,5)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
% V=V+0.005;
% surf(Kx,Ky,Enk(:,:,61)-efermi+0.005);
% contour3(Kx, Ky, Enk(:,:,61)-efermi+0.005, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);
V=0:-0.063/50:-0.3;
surf(Kx,Ky,Enk(:,:,4)-efermi);
contour3(Kx, Ky, Enk(:,:,4)-efermi, V, 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5);

colormap(slanCM('RdBu'))
shading interp

xlim([-0.3 0.3]);   % 限制 Y 轴范围
ylim([-0.18 0.18]);   % 限制 Y 轴范围
zlim([-0.3,0.3])
clim([-0.29 0.29]); 
view(-6,9)
%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Check Time Reversal Symmetry            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.3,0.2,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.3,-0.2,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
s2=[0  -1i
    1i  0];
T=kron(eye(4),i*s2);
h1=T*conj(hk1)*inv(T)-hk2; %% TH=HT T=i*sigma_y*K  THT^-1=H 
max(h1,[],'all')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Check glide mirror Mx  Symmetry            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.5,0.3,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.5,0.3,0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
% sigma1=[1  0
%         0  1];
% t0=[0 1
%     1 0];
sigma1=[1  0
        0  1];
t0=[1 0
    0 1];
s1=[0 1
    1 0];
M=kron(sigma1,t0)
% M=kron(t0,sigma1)
M=1j*kron(M,s1)

h1=M*hk1*inv(M)-hk2; %% TH=HT T=i*sigma_y*K  THT^-1=H 
max(h1,[],'all')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Check Inversion Symmetry         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];
kpoint=[0.5,0.0,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.5,0.0,0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
Pblk=blkdiag(d1,d2);
P=kron(fliplr(eye(n1)),Pblk);
%%
pos=eye(size(g.wpos,1));
kpoint=kpoint*g.b;
for i=1:size(g.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    % fprintf('%.6f\n',pos(i,i));
end
P=P*pos;
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
a=Psik'*P*Psik; %% P \psi=e \psi    e=\psi.T.conj P \psi
c=eig(a);
%%
a=eig(Psik(:,5:6)'*P*Psik(:,5:6));
sum(a)
%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(4,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
PnkZ2 = zeros(4,size(kpoints,1));
Z2 = zeros(size(kpoints,1),1);
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];
d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
P=blkdiag(d1,d2);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:)
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:4)'*P*Psik(:,1:4));
    Z4(i)=sum(Pnk(:,i));
    % PnkZ2=eig(Psik(:,1:2)'*P*Psik(:,1:2))
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check Inversion Symmetry for CDW     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(gs.ham);
s1=[0  1
    1  0];
kpoint=[0.3,0.2,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);%[Energy,Psik,hk]
kpoint=[-0.3,-0.2,0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);

d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
Pblk=blkdiag(d1,d2);
P=kron(fliplr(eye(n1)),Pblk);
%%
pos=eye(size(gs.wpos,1));
kpoint=kpoint*gs.b;
for i=1:size(gs.wpos)
    pos(i,i)=1;%exp(-2j*(g.wpos(i,:)-[1.0,0.0,0.5]*g.a)*kpoint');
    % fprintf('%.6f\n',pos(i,i));
end
P=P*pos;
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(nbands/2,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
% PnkZ2 = zeros(1,size(kpoints,1));
Z2 = zeros(size(kpoints,1),1);
s1=[0  1
    1  0];
d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
Pblk=blkdiag(d1,d2);
P=kron(fliplr(eye(n1)),Pblk);

for i=1:size(kpoints,1)
    kpoint=kpoints(i,:);
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:nbands/2)'*P*Psik(:,1:nbands/2));
    Z4(i)=sum(Pnk(:,i));
    Z2(i)=(-1)^(sum(abs(Pnk(:,i) + 1) <= 0.1)/2);
    % PnkZ2=eig(Psik(:,1:2)'*P*Psik(:,1:2))
end

%%

g2=initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
% g2=initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");
% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g2.ham = transformBasis(g2.ham, T);
s=1.02;
strain=[1/s,0,0;
        0.0,s,0.0;
        0,0,1];
strain=[1/s*1.00,0,0;
        0.0,s,0.0;
        0,0,1];

g2.a=strain*g2.a;
g2.wpos = setWannierPosition(g2);
[result_matrices2, pairsU0, pairsU, pairsV] = findNeighbors(g2);

[nbands,~,nrpts]=size(gs.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
% t_ij'=t_ij*e^{-\beta*(\delta a / a)}
% for i=2:size(result_matrices2,2)-1
for i=2:size(result_matrices2,2)-1
        delta_ij=result_matrices2{i}(:,6) - result_matrices{i}(:,6);
        t_ij=exp(-1.5 * delta_ij./result_matrices{i}(:,6));
        result_matrices2{i}=[result_matrices2{i},t_ij];
        num=0;
        for j=1:length(result_matrices2{i})
            bandindex_i=result_matrices2{i}(j,1);
            bandindex_j=result_matrices2{i}(j,2);
            scale=result_matrices2{i}(j,7);
            index=find(ismember(g2.hopr,result_matrices2{i}(j,3:5),"rows"));
            if ~isempty(index)
                num=num+1;
                g2.ham(2*bandindex_i-1:2*bandindex_i,2*bandindex_j-1:2*bandindex_j,index)=g2.ham(2*bandindex_i-1:2*bandindex_i,2*bandindex_j-1:2*bandindex_j,index)*scale;
            end
        end
end

g2.iniham = g2.ham + 0; 
efermi=0.022;
Electric_field_in_evpA=0.00*0.529177; nk=51;
[nbands,~,nrpts]=size(g2.ham);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g2.ham,g2.hopr,g2.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g2.a,g2.b);
% bandname="tit_hf_data/3x1/HFscanBand_"+int2str(n1)+"n_";
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
% plot(kpath,Energy(4*n1+1,:)-efermi,'Color','red','LineWidth',2);
% plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')

% % function obj=offsite_modify(obj,hopr,ham)
% % index=find(ismember(obj.hopr,hopr,'rows'));
% % ham_new=obj.ham(:,:,index)+ham;
% % obj.ham(:,:,index)=ham_new(:,:);
% % end


knum=501;
band1=1;
band2=4;
[wx,unk]=MTB.ham.get_wilsonloop(g2,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g2,knum,band1,band2);
%%
%%
kpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.0,0.5,0.0;0.5,0.5,0.0];
Unk = zeros(nbands,nbands,size(kpoints,1));
Pnk = zeros(4,size(kpoints,1));
Z4  = zeros(size(kpoints,1),1);
[nbands,~,nrpts]=size(g.ham);
s1=[0  1
    1  0];
d1=kron(s1,diag(repmat([1],1,2)));
d2=kron(s1,diag(repmat([1],1,2)));
P=blkdiag(d1,d2);
for i=1:size(kpoints,1)
    kpoint=kpoints(i,:)
    [Energy,Psik,hk]=MTB.ham.get_parity_singleK(g2.ham,g2.hopr,nbands,nrpts,kpoint,g2.a,g2.b);
    toc;
    Unk(:,:,i)=Psik;
    Pnk(:,i)=eig(Psik(:,1:4)'*P*Psik(:,1:4));
    Z4(i)=sum(Pnk(:,i));
end
%%
s_values = 0.98:0.01:1.02;  % 变化范围
Electric_field_in_evpA = 0.00 * 0.529177;  % 施加电场
nk = 51;  % k 点数量
efermi = 0.022;  % 费米能级
[nbands,~,nrpts]=size(gs.ham);
labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
% strain_list=0.98:0.01:1.05;
strain_list=0.98:0.01:1.05;
gap=strain_list;
n=0;
% 遍历应变参数 s
for s = 1:length(strain_list)
    n=n+1;
    g2=initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
    % g2=initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");
    % Step 2: Basis Transformation
    T = getBasisTransformMatrix();
    g2.ham = transformBasis(g2.ham, T);
    strain=[1/strain_list(s)*1.00,0,0;
        0.0,strain_list(s),0.0;
        0,0,1];
    g2.a=strain*g2.a;
    g2.wpos = setWannierPosition(g2);

    % 计算邻近矩阵
    [result_matrices2, pairsU0, pairsU, pairsV] = findNeighbors(g2);

    % 计算修正后的跳跃积分
    for i = 2:size(result_matrices2,2)-1
        delta_ij = result_matrices2{i}(:,6) - result_matrices{i}(:,6); % 确保 result_matrices2{i} 正确
        t_ij = exp(-1.5 * delta_ij ./ result_matrices{i}(:,6)); % 这里 beta = 3, 可调
        
        % 添加修正因子到矩阵
        result_matrices2{i} = [result_matrices2{i}, t_ij];

        % 遍历所有邻接项
        for j = 1:length(result_matrices2{i})
            bandindex_i = result_matrices2{i}(j,1);
            bandindex_j = result_matrices2{i}(j,2);
            scale = result_matrices2{i}(j,7);

            % 找到匹配的跳跃项
            index = find(ismember(g2.hopr, result_matrices2{i}(j,3:5), "rows"));
            if ~isempty(index)
                g2.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) = ...
                    g2.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) * scale;
            end
        end
    end

    % 记录初始哈密顿量
    g2.iniham = g2.ham + 0;

    % 计算能带
    [nbands,~,nrpts] = size(g2.ham);
    [Energy, kpath, kindex] = MTB.ham.get_bulk_bands_add_electric(...
        g2.ham, g2.hopr, g2.wpos, Electric_field_in_evpA, nbands, nrpts, hkpoints, nk, g2.a, g2.b);
    delta=min(Energy(5,:))-max(Energy(4,:));
    gap(s)=delta;
    % 绘制能带
    % hold on;
    % MTB.plot.plot_bands(Energy, nbands, efermi, kpath, labels, kindex, ...
        % strcat("TaIrTe4-2d (s=", num2str(s), ")"), Electric_field_in_evpA * 10000);
end

figure;
plot(strain_list,gap,'r-o')
% 
% % 运行计算并绘制不同应变下的能带
% strain_band_plot("TaIrTe4", "data/TaIrTe4_2d_tb/1019", s_values, Electric_field_in_evpA, nk, efermi, labels);

%%
filename="gap_strain.dat";
list=[strain_list'-1,gap']
writeoutput(filename,list)
%%
efermi=0.022;
Electric_field_in_evpA=0.00*0.529177; nk=51;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
% bandname="tit_hf_data/3x1/HFscanBand_"+int2str(n1)+"n_";
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1+1,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     DFT starin str_10             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_13/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_13/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_13/wannier90_hr_p2.dat');
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_8/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_8/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_8/wannier90_hr_p2.dat');
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_18/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_18/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_18/wannier90_hr_p2.dat');

g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/strain_x_relax_y/str_10/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/strain_x_relax_y/str_10/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/strain_x_relax_y/str_10/wannier90_hr_p2.dat');
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.45;
nk=101;

g.wpos=[]
g.wpos=g.atoms*g.a
orbital_num=[18,18,12,12,8,8,8,8,8,8,8,8]
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
        repmat(g.wpos(3,:),[orbital_num(3),1]);...
        repmat(g.wpos(4,:),[orbital_num(4),1]);...
        repmat(g.wpos(5,:),[orbital_num(5),1]);...
        repmat(g.wpos(6,:),[orbital_num(6),1]);...
        repmat(g.wpos(7,:),[orbital_num(7),1]);...
        repmat(g.wpos(8,:),[orbital_num(8),1]);...
        repmat(g.wpos(9,:),[orbital_num(9),1]);...
        repmat(g.wpos(10,:),[orbital_num(10),1]);...
        repmat(g.wpos(11,:),[orbital_num(11),1]);...
        repmat(g.wpos(12,:),[orbital_num(12),1]);...
    ]
g.wpos=g.wpos;
%%

Electric_field_in_evpA=0.00;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

hold on;
plot(kpath,Energy(28,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(89,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')
%%
knum=101;
band1=85;
band2=88;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
%%
knum=51;
band1=1;
band2=88;
[wx2,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
% [wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);

%%
kx=linspace(0,1,knum);
% filename="data/tit_hf/15x1/10nmd/U-epsilon-f/"+int2str(12)+"-15x1/wtool/r"+int2str(11)+"/wannier90_s1_eps"+int2str(13)+".00_r"+int2str(11)+"_wloopky_all.dat";
filename="data/TaIrTe4_2d/strain/strain_x_relax_y/str_10/wilsonloop/"+"wloopky.dat";
for iband=1:size(wx2,2)
   outlist=[kx',wx2(:,iband)];
   writeoutput(filename,outlist)
end

%%
knum=101;
band1=29;
band2=88;
% [wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%
knum=101;
band1=85;
band2=88;
% [wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     DFT starin str_16             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_13/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_13/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_13/wannier90_hr_p2.dat');
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_8/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_8/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_8/wannier90_hr_p2.dat');
% g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/str_18/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/str_18/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/str_18/wannier90_hr_p2.dat');

g = MTB.read_poscar(g,"data/TaIrTe4_2d/strain/strain_x_relax_y/str_16/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/strain/strain_x_relax_y/str_16/wannier90_hr_p1.dat','data/TaIrTe4_2d/strain/strain_x_relax_y/str_16/wannier90_hr_p2.dat');
%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Y'}; % labels for k
% hkpoints={[0.5,0.5,0.0],...
%           [0.0,0.0,0.0],...
%           [0.0,0.5,0.0],...
%           [0.5,0.0,0.0]};% hkpoints-high symmetry k points

labels={'R','Y','\Gamma','X'}; % labels for k
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=-0.65;
nk=20;

g.wpos=[]
g.wpos=g.atoms*g.a
orbital_num=[18,18,12,12,8,8,8,8,8,8,8,8]
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
        repmat(g.wpos(3,:),[orbital_num(3),1]);...
        repmat(g.wpos(4,:),[orbital_num(4),1]);...
        repmat(g.wpos(5,:),[orbital_num(5),1]);...
        repmat(g.wpos(6,:),[orbital_num(6),1]);...
        repmat(g.wpos(7,:),[orbital_num(7),1]);...
        repmat(g.wpos(8,:),[orbital_num(8),1]);...
        repmat(g.wpos(9,:),[orbital_num(9),1]);...
        repmat(g.wpos(10,:),[orbital_num(10),1]);...
        repmat(g.wpos(11,:),[orbital_num(11),1]);...
        repmat(g.wpos(12,:),[orbital_num(12),1]);...
    ]
g.wpos=g.wpos;
%%

Electric_field_in_evpA=0.00;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

hold on;
plot(kpath,Energy(28,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(89,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')
%%
%%
% filename="data/TaIrTe4_2d/strain/strain_x_relax_y/str_16/bands.dat";
% for iband=1:nbands
%    outlist=[kpath',Energy(iband,:)'];
%    writeoutput(filename,outlist)
% end
%%
knum=51;
band1=85;
band2=88;
[wx2,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%
knum=101;
band1=1;
band2=88;
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%
knum=101;
band1=29;
band2=88;
% [wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%
knum=101;
band1=85;
band2=88;
% [wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%
%%
kx=linspace(0,1,knum);
% filename="data/tit_hf/15x1/10nmd/U-epsilon-f/"+int2str(12)+"-15x1/wtool/r"+int2str(11)+"/wannier90_s1_eps"+int2str(13)+".00_r"+int2str(11)+"_wloopky_all.dat";
filename="data/TaIrTe4_2d/strain/strain_x_relax_y/str_16/wilsonloop/"+"wloopky.dat";
for iband=1:size(wx2,2)
   outlist=[kx',wx2(:,iband)];
   writeoutput(filename,outlist)
end
%%
%U0 onsite orbital 1 2 3 4
%          pairsU{1}~pairsU{12}
%U onsite_nn 13&31 24&42 offsite_nn 13&31 24&42 offsite_nnn 12&21 14&41 32&32 34&43 offsite_nnnn 12&21 14&41 %32&32 34&43
%V onsite_nn 13&31 24&42 offsite_nn 13&31 24&42 offsite_nnn 12&21 14&41 32&32 34&43 offsite_nnnn 12&21 14&41 %32&32 34&43
% Generate random U and V
% % U = rand(1, 4); % Random U values (4 elements)
% % V = rand(1, 3) * -1; % Random V values (3 elements, negative values)
% Generate random U0 with constraints
U0 = rand()*0.8 + 0.3; % U1 is the largest value, within [0.3, 1.1]
perturbation_U = 1-(rand(1, 4) - 0.5) * 0.3; %[+- 0.15%]
U0 = [U0,U0,U0,U0].*perturbation_U*2;
% Generate random U1 with constraints
shared_value_onsite_nn = rand() * 0.6 + 0.05; % U2-U4 are close, within [0.05, 0.66]
pertu1=1-(rand(1, 2) - 0.5) * 0.3; %[%[+- 0.15%]]
shared_value_onsite_nn=[shared_value_onsite_nn,shared_value_onsite_nn].*pertu1;
shared_value_nn = rand() * 0.5 + 0.05; % U2-U4 are close, within [0.05, 0.55]
pertu2=1-(rand(1, 2) - 0.5) * 0.3; %[%[+- 0.15%]]
shared_value_nn=[shared_value_nn,shared_value_nn].*pertu2;
shared_value_nnn = rand() * 0.4 ; % U2-U4 are close, within [0.0, 0.4]
pertu3=1-(rand(1, 4) - 0.5) * 0.3; %[%[+- 0.15%]]
shared_value_nnn=[shared_value_nnn,shared_value_nnn,shared_value_nnn,shared_value_nnn].*pertu3;
shared_value_nnnn = rand() * 0.3; % U2-U4 are close, within [0.0, 0.3]
pertu4=1-(rand(1, 4) - 0.5) * 0.3; %[%[+- 0.15%]]
shared_value_nnnn = [shared_value_nnnn,shared_value_nnnn,shared_value_nnnn,shared_value_nnnn].*pertu4;

U = [shared_value_onsite_nn,...
     shared_value_nn,...
     shared_value_nnn,...
     shared_value_nnnn]*1.2; % U values

% Generate V based on U with perturbation
perturbation = (rand(1, length(U)) - 0.5) * 0.3; % Small perturbation in range [-0.1, 0.1]
V = -(abs(U) + perturbation); % V values are close to -U(2:4)
%%
% % U0=[ 0.7207,   0.6844];
% % % U=[ 0.3255,    0.3137,   0.3520,   0.320,    0.289,  0.3539];
% % % V=[-0.3459  , -0.2663,  -0.3328,  -0.398,   -0.3085, -0.3440];
% % U=[ 0.3,    0.3, 0.3520,   0.320,    0.289,  0.3539];
% % V=[-0.3, -0.2, -0.3328,  -0.398,   -0.3085, -0.3440];
% xinitial_0 = diag(kron(ones(1, n1), [nec1 - 0.4, nec1 - 0.4, nec2 + 0.4, nec2 + 0.4, ...
%                                     nec1 + 0.4, nec1 + 0.4, nec2 - 0.4, nec2 - 0.4]));
%%
% Initialize Hartree-Fock States

[xinitial_U, xinitial_V] = initializeStates(U, V, pairsU, pairsV, nbands);
% xinitial = {xinitial_0, xinitial_U, xinitial_V};
% % xinitial_U{1}=result.xinitial{2}{1};
% % xinitial_V{1}=result.xinitial{3}{1};
% % xinitial_U{2}=result.xinitial{2}{2};
% % xinitial_V{2}=result.xinitial{3}{2};
% % xinitial_0=result.xinitial{1};
xinitial = {xinitial_0, xinitial_U, xinitial_V};
% Load and Save
%save('xinitial_good.mat',"xinitial","efermi","-v7.3")
% load('xinitial_good.mat')
%%

% Run Hartree-Fock
% u1 = 1e-10; u2 = (4 * n1 + 1.9) / nbands;
u1 = 1e-10; u2 = (4 * n1 * n2 + 2) / nbands;
knum=30;
stepmax=100;
parpool('Threads')
[xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, 0, ...
    xinitial, stepmax, critial, U0, U, V, u1, u2, pairsU0, pairsU, pairsV);



%%

[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,0,Kx,Ky,Kz);
Enk=reshape(Enk,[knum^2,nbands]);
Unk=reshape(Unk,[nbands,nbands,knum^2]);
[~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
[ni,si]=calonsite(Unk,kindex,bandindex);
%%
% efermi=0
Electric_field_in_evpA=0.00*0.529177; nk=51;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
% bandname="tit_hf_data/3x1/HFscanBand_"+int2str(n1)+"n_";
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1+1,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')

%% Reload and plt data
%load('RandomUV_HartreeFock_results.mat')
result=load("/Users/jxli/work/test/sparse/15x1/4nmd/11-15x1/result_1_U4.mat")
%%
idx=22
U0=results(idx).U0
U=results(idx).U
V=results(idx).V
xinitial=results(idx).xinitial
efermi=results(idx).efermi  

%% Rerun HF
d=10*10^(-9); %10nm
epsilon=10;
n_max=10^7;
% ra=1.895*10^(-10); % r=0.7
ra=0.85*10^(-10); % r=0.7
Ua=dual_gate_potential(ra, d, epsilon, n_max)
rnn=3.77*10^(-10); % r=0.7A
v1=dual_gate_potential(rnn, d, epsilon, n_max)
rnnn=4.4*10^(-10); % r=0.7A
v2=dual_gate_potential(rnnn, d, epsilon, n_max)
rnnnn=6.91*10^(-10); % r=0.7A
v3=dual_gate_potential(rnnnn, d, epsilon, n_max)
rnnnnn=7.54*10^(-10); % r=0.7A
v4=dual_gate_potential(rnnnnn, d, epsilon, n_max)
r5=8.66*10^(-10);
v5=dual_gate_potential(r5, d, epsilon, n_max)
r6=10.17*10^(-10);
v6=dual_gate_potential(r6, d, epsilon, n_max)
%%
Electric_field_in_evpA=0.00;
Ua=Ua;
Ub=Ua;
% Ua=
U0=[ Ua Ua Ub Ub]; 

v_intra_13=Ua;
v_intra_24=v_intra_13;
% v1=0
v_inter_nn_13=v1;
v_inter_nn_24=v1;
% v2=0
v_inter_nnn_12=v2;
v_inter_nnn_14=v2;
v_inter_nnn_32=v_inter_nnn_14;
v_inter_nnn_34=v2;
% v3=0;
v_inter_nnnn_12=v3;
v_inter_nnnn_14=v3;
v_inter_nnnn_32=v_inter_nnnn_14;
v_inter_nnnn_34=v3;
% v4=0
v_inter_nnnnn_13=v4;
v_inter_nnnnn_24=v4;
% 1-3 2-4 1-2 1-4 3-2 3-4 1-2 1-4 3-2 3-4
V=[v_intra_13,v_intra_24,...
    v_inter_nn_13,v_inter_nn_24,...
    v_inter_nnn_12,v_inter_nnn_14,v_inter_nnn_32,v_inter_nnn_34,...
    v_inter_nnnn_12,v_inter_nnnn_14,v_inter_nnnn_32,v_inter_nnnn_34,...
    v_inter_nnnnn_13,v_inter_nnnnn_24,...
    v5,v6];

% V=[ 0.28      0.28    0.20    0.22    0.3091    0.2835    0.0856    0.0837    0.0938    0.0812] % 1-3 2-4 1-2 1-4 3-2 3-4 1-2 1-4 3-2 3-4
%U=[ 0.4834     0.4625    0.4422    0.4965    0.4214    0.4834    0.3462    0.3618    0.3716    0.3309]*0.5
%V=[-0.4446    -0.4435   -0.5235   -0.3551   -0.4465   -0.4126   -0.3102   -0.5033   -0.4112   -0.2158]*1


%%
%
% xinitial=result.xinitial
% Run Hartree-Fock
% xinitial_0 = diag(kron(ones(1, n1), [nec1 - 0.46, nec1 - 0.46, nec2 + 0.46, nec2 + 0.46, ...
%                                      nec1 + 0.46, nec1 + 0.46, nec2 - 0.46, nec2 - 0.46]))*0.5;
xinitial_0 = diag(kron(ones(1, 4), rand(1,60)));
[xinitial_V1, xinitial_V2] = initializeStates(V, V, pairsU, pairsV, nbands,xinitial_0);
xinitial_0=xinitial_0.*U0(1);
% xinitial_0 = diag(kron(ones(1, n1), rand(1,8)));
% xinitial_0 = diag(kron(ones(1, n1), 1));
xinitial = {xinitial_0, xinitial_V1, xinitial_V2};
%%
%%
% profile on
Electric_field_in_evpA=0;
[xinitial, metaData] = flattenNestedCell(xinitial);        % 展平操作
xinitial = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, V, V, u1, u2, pairsU0, pairsU, pairsV);
xinitial = restoreNestedCell(xinitial, metaData);  % 复原操作
modifyHam(gs, xinitial, V, V, pairsU, pairsV)  %修改Ham
nbands=size(gs.ham,1);

%%
[Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
%%
knum=200
Enk=reshape(Enk,[knum^2,nbands]);
Unk=reshape(Unk,[nbands,nbands,knum^2]);
[~,kindex,bandindex,efermi]=Total_energy(Enk,0.37);
[ni,si]=calonsite(Unk,kindex,bandindex);
%%
%%
knum=15;
stepmax=60;
critial=1e-7;
% parpool('Threads')
u1 = 1e-10; u2 = (4 * n1 + 2) / nbands;

[xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, 0, ...
    xinitial, stepmax, critial, U0, V, V, u1, u2, pairsU0, pairsU, pairsV);

%%
% efermi=0;
% modifyHam(gs, xinitial, V, V, pairsU, pairsV)
Electric_field_in_evpA=0.00*0.529177; nk=51;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
% bandname="tit_hf_data/3x1/HFscanBand_"+int2str(n1)+"n_";
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
plot(kpath,Energy(4*n1*n2+1,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(4*n1*n2+2,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')

%%
    result=struct();
    result.U0=U0;
    result.U=V;
    result.V=V;
    result.xinitial=xinitial;
    result.ni=ni;
    result.si=si;
    result.efermi=efermi;
    result.pairsU=pairsU;
    result.pairsV=pairsV;
save('tit_hf_data/target/result_scf_v1.mat',"result")
%%
load('tit_hf_data/target/result_scf_v1.mat')
%%
load('/Users/jxli/work/tb/matlab/tit_hf_data/target/result_3.mat')
%%
    % result.U0=U0;
    % result.U=U;
    % result.V=V;
    % result.xinitial=xinitial;
    % result.ni=ni;
    % result.si=si;
    % result.efermi=efermi;
    % save('tit_hf_data/1-15x1/result_452_new.mat',"result")
%%
% result=load("/Users/jxli/work/test/sparse/15x1/4nmd/11-15x1/result_1_U4.mat");
clc;
clear;
result=load("/Users/jxli/work/test/sparse/15x1/10nmd/6-15x1/result_1_U1.mat");
U0=result.U0
% U=result.
V=result.V
U=V
xinitial=result.xinitial
efermi=result.efermi
%%
knum = 60;
kxline = [0, 1]; kyline = [0, 1];
[Kx, Ky, Kz] = gs.get_Bulk2Dkmesh(kxline, kyline, knum);
kpoints = [Kx(:), Ky(:), Kz(:)];
stepmax = 50; critial = 1e-9; % HF parameters
u1 = 1e-10; u2 = (4 * n1 + 2) / nbands;

[xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, 0, ...
    xinitial, stepmax, critial, U0, U, V, u1, u2, pairsU0, pairsU, pairsV);
%%
% result=load("tit_hf_data/1-15x1/result_452.mat");
% result=load("/Users/jxli/work/test/sequence/test/1-15x1/result_1.mat")
% result=load("/Users/jxli/work/test/sparse/1-15x1/result_1.mat")
 % result=load("/Users/jxli/work/test/sparse/1-15x1/result_205.mat")
% load("/Users/jxli/work/test/sparse/15x1/4nmd/11-15x1/result_1_U4.mat")
% result=load("/Users/jxli/work/test/sparse/15x1/10nmd/15-15x1/result_1_U14.mat")
result=load("/Users/jxli/work/test/sparse/15x1/10nmd/pargram/u-epsilon/newbasis/diff/1-15x1/result_s2_eps2.00_r9.mat")
U0=result.U0
U=result.V
V=result.V
xinitial=result.xinitial;
pairsU=result.pairsU;
pairsV=result.pairsV;
efermi=result.efermi;
%
% xinitial{2}{1}=xinitial{2}{1}*3;
%%
xinitial{2}{1}=xinitial{2}{1}*1.0;
xinitial{2}{2}=xinitial{2}{2}*1.0;
xinitial{2}{3}=xinitial{2}{3}*0.0;
xinitial{2}{4}=xinitial{2}{4}*0.0;
xinitial{2}{5}=xinitial{2}{5}*0.0;
xinitial{2}{6}=xinitial{2}{6}*0.0;
xinitial{2}{7}=xinitial{2}{7}*0.0;
xinitial{2}{8}=xinitial{2}{8}*0.0;
xinitial{2}{9}=xinitial{2}{9}*0.0;
xinitial{2}{10}=xinitial{2}{10}*0.0;

xinitial{3}{1}=xinitial{3}{1}*1.0;
xinitial{3}{2}=xinitial{3}{2}*1.0;
xinitial{3}{3}=xinitial{3}{3}*0.0;
xinitial{3}{4}=xinitial{3}{4}*0.0;
xinitial{3}{5}=xinitial{3}{5}*0.0;
xinitial{3}{6}=xinitial{3}{6}*0.0;
xinitial{3}{7}=xinitial{3}{7}*0.0;
xinitial{3}{8}=xinitial{3}{8}*0.0;
xinitial{3}{9}=xinitial{3}{9}*0.0;
xinitial{3}{10}=xinitial{3}{10}*0.0;
%%
% modifyHam(gs, xinitial, V, V, pairsU, pairsV)
% 
% % xinitial_0=zeros(nbands,nbands)
% % % xinitial_0=diag(kron(ones(1,n1),kron(U0,ones(1,2))));
% % gs.onsite_modify(xinitial_0)
% % efermi=0
    % % 
    % % [a,b]=mink(Enk(:),ceil(size(Enk(:),1)*u));
    % % efermi = max(a,[],"all");


labels={'Y','\Gamma','X','R','\Gamma','X','Y','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.0,0.5,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

labels={'X','\Gamma','R','Y'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points

labels={'Y','\Gamma','X','R','\Gamma'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points

Electric_field_in_evpA=0.00*0.529177;
nk=101;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
ylim([-1,1])
hold on;
% plot(kpath,Energy(4*n1,:)-efermi,'Color','magenta','LineWidth',2);
plot(kpath,Energy(4*n1+1,:)-efermi,"Color",'red','LineWidth',2);
plot(kpath,Energy(4*n1+2,:)-efermi,"Color",'blue','LineWidth',2);
% plot(kpath,Energy(4*n1+3,:)-efermi,"Color",'yellow','LineWidth',2);
%%
knum=201;
band1=1;
band2=4*n1;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
% [wx,unk]=MTB.ham.get_wilsonloop_ky(gs,knum,band1,band2);
%%
%% Calculate Dos
nk=knum;
Enum=100;
Emin=-0.1;
Emax= 0.1;
eps=(Emax-Emin)/Enum;
plottap=1;
Tem=8;

knum = 200;
kxline = [0, 1]; kyline = [0, 1];
[Kx, Ky, Kz] = gs.get_Bulk2Dkmesh(kxline, kyline, knum);
%%
[~,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
% Enk=reshape(Enk,[knum^2,nbands]);
% Unk=reshape(Unk,[nbands,nbands,knum^2]);
%%
[Eaxis,Dos,TDos]=MTB.ham.get_dos_FermiDirac(Enk-efermi,Tem,Enum,Emin,Emax,nk,plottap);
% [Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk-efermi,eps,Enum,Emin,Emax,nk,plottap);
%%
figure()
plot((Eaxis-0.015)*1000,Dos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
xlim([-0.08 0.08]*1000)
% ylim([0 30])
xlabel("meV")
ylabel("DOS")
% plot(Eaxis+0.014,TDos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
%plot(Eaxis,TDos/norm(cross(g.a(1,:),g.a(2,:)))*10^16,'Linestyle','-','Color','#4DA1D7','LineWidth',2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Write the wannier90_hr.dat       %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
filename="/Users/jxli/work/test/sparse/15x1/wannier90_15cdw_newbasis.dat";
MTB.write_hr(gs,filename)
%%
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,0,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,0.5);


%%
knum=30;
band1=1;
band2=4*n1;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
knum=30;
band1=4*n1+1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);
%%
knum=100;
band1=1;
band2=4*n1+2;
[wx,unk]=MTB.ham.get_wilsonloop(gs,knum,band1,band2);


%% SCAN
clc;
clear;
% profile on; 
% Step 1: Initialize Geometry and Hamiltonian
g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb/1019");
% g = initializeGeometry("TaIrTe4", "data/TaIrTe4_2d_tb");

% Step 2: Basis Transformation
T = getBasisTransformMatrix();
g.ham = transformBasis(g.ham, T);

% Step 3: Set Wannier Position
g.wpos = setWannierPosition(g);

% Step 4: Construct Supercell
n1 = 15; n2 = 1; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

gs.iniham = gs.ham + 0; % Store the initial Ham
[nbands,~,nrpts]=size(gs.ham);
labels={'Y','\Gamma','X','R','\Gamma'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0],...
          [0.0,0.0,0.0]};% hkpoints-high symmetry k points

% Step 5: Neighbor Search
[result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs);
%
% Step 6: Initialize fixed Hartree-Fock States parameters
nec1 = (4 * n1 * n2 + 4) / n1 / 2 / 2 / 2;
nec2 = (4 * n1 * n2) / n1 / 2 / 2 / 2;
% xinitial_0 = diag(kron(ones(1, n1), [nec1 - 0.4, nec1 - 0.4, nec2 + 0.4, nec2 + 0.4, ...
%                                     nec1 + 0.4, nec1 + 0.4, nec2 - 0.4, nec2 - 0.4]));


knum = 20;
kxline = [0, 1]; kyline = [0, 1];
[Kx, Ky, Kz] = gs.get_Bulk2Dkmesh(kxline, kyline, knum);
kpoints = [Kx(:) Ky(:), Kz(:)];
stepmax = 50; critial = 1e-5; % HF parameters
u1 = 1e-10; u2 = (4 * n1 *n2 + 2) / nbands;

% Parameters for scanning
num_scans = 1; % Number of rlansdom combinations to scan
results = struct(); % Structure to store results

% Start total time
total_time = tic;
    s=0.22; %%0.22
% Loop over random combinations of U and V
for idx = 1:num_scans
    fprintf('Scan %d started.\n', idx);

    % Start timer for this scan
    scan_start = tic;
    result=struct();
    %% Initial From dual gate potential

    Ua=0.75*s;
    v1=0.3720*s;
    v2=0.3173*s;
    v3=0.1984*s;
    v4=0.1810*s;
    v5=0.1563*s;
    v6=0.1316*s;
    U0=[ Ua Ua Ua Ua];
    v_intra_13=Ua;
    v_intra_24=v_intra_13;
    % v1=0
    v_inter_nn_13=v1;
    v_inter_nn_24=v1;
    % v2=0
    v_inter_nnn_12=v2;
    v_inter_nnn_14=v2;
    v_inter_nnn_32=v2;
    v_inter_nnn_34=v2;
    % v3=0;
    v_inter_nnnn_12=v3;
    v_inter_nnnn_14=v3;
    v_inter_nnnn_32=v3;
    v_inter_nnnn_34=v3;
    % v4=0
    v_inter_nnnnn_13=v4;
    v_inter_nnnnn_24=v4;
    % 1-3 2-4 1-2 1-4 3-2 3-4 1-2 1-4 3-2 3-4
    V=[v_intra_13,v_intra_24,...
        v_inter_nn_13,v_inter_nn_24,...
        v_inter_nnn_12,v_inter_nnn_14,v_inter_nnn_32,v_inter_nnn_34,...
        v_inter_nnnn_12,v_inter_nnnn_14,v_inter_nnnn_32,v_inter_nnnn_34,...
        v_inter_nnnnn_13,v_inter_nnnnn_24,...
        v5,v6];

    % % % % Initialize Hartree-Fock States
    xinitial_0 = diag(kron(ones(1, 1), rand(1,120))); %for 15x1
    % xinitial_0 = diag(kron(ones(1, 7), rand(1,8))); %for 7x1
    % xinitial_0 = diag(kron(ones(1, n1), rand(1,8))); %for 3x1
    [xinitial_V1, xinitial_V2] = initializeStates(V, V, pairsU, pairsV, nbands,xinitial_0);
    xinitial_0=xinitial_0.*U0(1);
    xinitial = {xinitial_0, xinitial_V1, xinitial_V2};
    % load("tit_hf_data/3x1/result_1.mat")


    % Run Hartree-Fock
    [xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, 0, ...
                                              xinitial, stepmax, critial, U0, V, V, u1, u2, pairsU0, pairsU, pairsV);

    Electric_field_in_evpA=0.00*0.529177; nk=51;
    [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(gs.ham,gs.hopr,gs.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,gs.a,gs.b);
    bandname="tit_hf_data/3x1/HFscanBand_"+int2str(n1)+"n_"+int2str(idx);
    MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
    hold on;
    plot(kpath,Energy(4*n1*n2+1,:)-efermi,'Color','magenta','LineWidth',2);
    plot(kpath,Energy(4*n1*n2+2,:)-efermi,"Color",'red','LineWidth',2);
    print(bandname,'-dpng','-r600')
    % 
    knum2=21;
    band1=1;
    band2=4*n1*n2;
    [wx1,~]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
    wxname="tit_hf_data/3x1/HFscanZ2_VB"+int2str(n1)+"n_"+int2str(idx);
    print(wxname,'-dpng','-r600')
    % 
    % knum2=51;
    % band1=4*n1+1;
    % band2=4*n1+2;
    % [wx2,~]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
    % wxname="tit_hf_data/3x1/HFscanZ2_CB"+int2str(n1)+"n_"+int2str(idx);
    % print(wxname,'-dpng','-r600')
    % 
    knum2=21;
    band1=1;
    band2=4*n1*n2+2;
    [wx3,~]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
    wxname="tit_hf_data/3x1/HFscanZ2_VCB"+int2str(n1)+"n_"+int2str(idx);
    print(wxname,'-dpng','-r600')

    knum2=21;
    band1=4*n1*n2+1;
    band2=4*n1*n2+2;
    [wx3,~]=MTB.ham.get_wilsonloop(gs,knum2,band1,band2);
    wxname="tit_hf_data/3x1/HFscanZ2_VCB"+int2str(n1)+"n_"+int2str(idx);
    print(wxname,'-dpng','-r600')
    % 
    % Store results for this random U and V
    % % % results(idx).U0 = U0;
    % % % results(idx).U = U;
    % % % results(idx).V = V;
    % % % results(idx).xinitial = xinitial;
    % % % results(idx).ni = ni;
    % % % results(idx).si = si;
    % % % results(idx).efermi = efermi;
    % results(idx).wx1 = wx1;
    % results(idx).wx2 = wx2;
    % results(idx).wx3 = wx3;

    % Record and display the time for this scan

    % save 
    result.U0=U0;
    % result.V1=V;
    result.V=V;
    result.xinitial=xinitial;
    result.ni=ni;
    result.si=si;
    result.efermi=efermi;
    result.pairsU=pairsU;
    result.pairsV=pairsV;

    % pathn="tit_hf_data/3x1/result"+int2str(n1)+"n_"+int2str(idx)+".mat";
    % save('pathn',"-fromstruct",result)
    save(sprintf("tit_hf_data/3x1/result_%d.mat",idx),"-fromstruct",result);

    scan_time = toc(scan_start);
    fprintf('Scan %d completed in %.2f seconds.\n', idx, scan_time);
end


% Record and display total runtime
total_runtime = toc(total_time);
fprintf('Total runtime for all scans: %.2f seconds.\n', total_runtime);

% Save results to a .mat file for later analysis
% save('RandomUV_HartreeFock_results.mat', 'results','-v7.3');

disp('Random U and V Hartree-Fock calculations completed.');

% profile viewer;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Check Time Reversal Symmetry            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k
kpoint=[0.3,0.1,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
kpoint=[-0.3,-0.1,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(gs.ham,gs.hopr,nbands,nrpts,kpoint,gs.a,gs.b);
s2=[0  -1i
    1i  0];
T=kron(eye(60),i*s2);

h1=T*conj(hk1)*inv(T)-hk2; %% TH=HT T=i*sigma_y*K  THT^-1=H 
max(h1,[],'all')
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                              初始几何和Ham                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = initializeGeometry(materialName, dataPath)
    % initializeGeometry Initializes the material geometry and Hamiltonian
    %
    % Inputs:
    %   materialName - String, the name of the material (e.g., "TaIrTe4")
    %   dataPath     - String, the path to the directory containing input files
    %
    % Outputs:
    %   g            - Struct, contains the geometry and Hamiltonian information
    %
    % This function reads the material geometry from a POSCAR file and
    % constructs the Hamiltonian using Wannier90 HR files.

    % Create geometry object
    g = MTB.geometry(materialName);
    % Read geometry from POSCAR file
    g = MTB.read_poscar(g, fullfile(dataPath, "POSCAR"));
    % Read Wannier Hamiltonian data
    [g.ham, g.hopr] = MTB.wannier.read_hr(...
        fullfile(dataPath, "wannier90_hr_p1.dat"), ...
        fullfile(dataPath, "wannier90_hr_p2.dat"));
end

function T = getBasisTransformMatrix()
    % getBasisTransformMatrix Returns the basis transformation matrix
    %
    % Outputs:
    %   T - Matrix, the basis transformation matrix for Hamiltonian
    %
    % The transformation matrix `T` is used to change the basis of the
    % Hamiltonian matrices.
    % % 
    % % T = [1,0,0,0,0,0,0,0; ...
    % %      0,0,0,0,0,1,0,0; ...
    % %      0,0,1,0,0,0,0,0; ...
    % %      0,0,0,0,0,0,0,1; ...
    % %      0,1,0,0,0,0,0,0; ...
    % %      0,0,0,0,1,0,0,0; ...
    % %      0,0,0,1,0,0,0,0; ...
    % %      0,0,0,0,0,0,1,0];

    T = [1,0,0,0,0,0,0,0; ...
         0,0,0,0,1,0,0,0; ...
         0,1,0,0,0,0,0,0; ...
         0,0,0,0,0,1,0,0; ...
         0,0,1,0,0,0,0,0; ...
         0,0,0,0,0,0,1,0; ...
         0,0,0,1,0,0,0,0; ...
         0,0,0,0,0,0,0,1];
end

function ham = transformBasis(ham, T)
    % transformBasis Applies a basis transformation to the Hamiltonian
    %
    % Inputs:
    %   ham - 3D array, the original Hamiltonian in the initial basis
    %   T   - Matrix, the basis transformation matrix
    %
    % Outputs:
    %   ham - 3D array, the Hamiltonian in the transformed basis

    % Apply the transformation for each Hamiltonian slice
    for i = 1:size(ham, 3)
        ham(:, :, i) = T * ham(:, :, i) * inv(T);
    end
end

function wpos = setWannierPosition(g)
    % setWannierPosition Sets the Wannier positions for the geometry
    %
    % Inputs:
    %   g    - Struct, contains the geometry and atomic positions
    %
    % Outputs:
    %   wpos - Array, the Wannier positions for each orbital
    %
    % This function computes the Wannier positions based on the atomic
    % positions and lattice vectors.

    % Compute Wannier positions based on atomic positions
    wpos = g.atoms * g.a;
    wpos = [wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :); ...
            wpos(1, :); wpos(1, :); wpos(2, :); wpos(2, :)];
end

function gs = constructSupercell(g, n1, n2)
    % constructSupercell Constructs the supercell Hamiltonian
    %
    % Inputs:
    %   g  - Struct, contains the geometry and Hamiltonian of the material
    %   n1 - Integer, number of cells along the first lattice vector
    %   n2 - Integer, number of cells along the second lattice vector
    %
    % Outputs:
    %   gs - Struct, the supercell Hamiltonian and related properties

    % Generate supercell Hamiltonian
    gs = MTB.ham.get_supercell_wannier(g, n1, n2);
    % Store the initial Hamiltonian for reference
    gs.iniham = gs.ham+0;
end


function [result_matrices, pairsU0, pairsU, pairsV] = findNeighbors(gs)
    % findNeighbors Finds neighboring atoms for interactions
    %
    % Inputs:
    %   gs - Struct, contains the supercell Hamiltonian and Wannier positions
    %
    % Outputs:
    %   result_matrices - Cell array, contains neighbor data for various distances
    %   pairsU          - Cell array, contains pairs for on-site and off-site interactions (U)
    %   pairsV          - Cell array, contains pairs for on-site and off-site interactions (V)

    result_matrices = cell(1, 7);
    % Loop through different neighbor distances
    for i = 1:7
        result_matrices{i} = find_neighbor_data(gs.wpos(1:2:end, :), gs.a, i, 2);
    end

    % Filter and group pairs
    pairs_onsite = result_matrices{1}(result_matrices{1}(:, 1) == result_matrices{1}(:, 2), :);
    pairs_onsite_nn = result_matrices{1}(result_matrices{1}(:, 1) ~= result_matrices{1}(:, 2), :);

    pairs_onsite_1 = pairs_onsite(mod(pairs_onsite(:,1),4)==1,:);
    pairs_onsite_2 = pairs_onsite(mod(pairs_onsite(:,1),4)==2,:);
    pairs_onsite_3 = pairs_onsite(mod(pairs_onsite(:,1),4)==3,:);
    pairs_onsite_4 = pairs_onsite(mod(pairs_onsite(:,1),4)==0,:);
    pairs_onsite_nn_13 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==1,:);
    pairs_onsite_nn_24 = pairs_onsite_nn(mod(pairs_onsite_nn(:,1),2)==0,:);


    pairs_offsite_nn = result_matrices{2};
    pairs_offsite_nn_13=result_matrices{2}(mod(result_matrices{2}(:, 1),2)~=0 & mod(result_matrices{2}(:, 2),2)~=0, :);
    pairs_offsite_nn_24=result_matrices{2}(mod(result_matrices{2}(:, 1),2)==0 & mod(result_matrices{2}(:, 2),2)==0, :);


    pairs_offsite_nnn = result_matrices{3};
    pairs_offsite_nnn_12=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==2,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_14=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 1 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==1,:)];
    pairs_offsite_nnn_34=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==0,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 0 & mod(result_matrices{3}(:,2),4)==3,:)];
    pairs_offsite_nnn_32=[result_matrices{3}(mod(result_matrices{3}(:,1),4)== 2 & mod(result_matrices{3}(:,2),4)==3,:);...
                          result_matrices{3}(mod(result_matrices{3}(:,1),4)== 3 & mod(result_matrices{3}(:,2),4)==2,:)];

    pairs_offsite_nnnn = result_matrices{4};
    pairs_offsite_nnnn_12=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==2,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_14=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 1 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==1,:)];
    pairs_offsite_nnnn_34=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==0,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 0 & mod(result_matrices{4}(:,2),4)==3,:)];
    pairs_offsite_nnnn_32=[result_matrices{4}(mod(result_matrices{4}(:,1),4)== 2 & mod(result_matrices{4}(:,2),4)==3,:);...
                          result_matrices{4}(mod(result_matrices{4}(:,1),4)== 3 & mod(result_matrices{4}(:,2),4)==2,:)];

    pairs_offsite_nnnnn = result_matrices{5};
    pairs_offsite_nnnnn_13=result_matrices{5}(mod(result_matrices{5}(:, 1),2)~=0 & mod(result_matrices{5}(:, 2),2)~=0, :);
    pairs_offsite_nnnnn_24=result_matrices{5}(mod(result_matrices{5}(:, 1),2)==0 & mod(result_matrices{5}(:, 2),2)==0, :);


    % Construct pairs for U and V
    pairsU0 = {pairs_onsite_1, pairs_onsite_2, pairs_onsite_3, pairs_onsite_4};
    pairsU = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34,...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
    pairsV = {pairs_onsite_nn_13, pairs_onsite_nn_24,...
              pairs_offsite_nn_13,pairs_offsite_nn_24,...
              pairs_offsite_nnn_12, pairs_offsite_nnn_14,...
              pairs_offsite_nnn_32, pairs_offsite_nnn_34...
              pairs_offsite_nnnn_12,pairs_offsite_nnnn_14,...
              pairs_offsite_nnnn_32,pairs_offsite_nnnn_34,...
              pairs_offsite_nnnnn_13,pairs_offsite_nnnnn_24,...
              result_matrices{6},result_matrices{7}};
end

function [xinitial, ni, si, efermi] = runHartreeFock(gs, xinitial, U, V, knum, stepmax, critial)
    % runHartreeFock Runs the Hartree-Fock self-consistent calculation
    %
    % Inputs:
    %   gs       - Struct, the supercell Hamiltonian and related properties
    %   xinitial - Cell array, initial guesses for Hartree-Fock states
    %   U        - Array, on-site Coulomb interaction values
    %   V        - Array, off-site Coulomb interaction values
    %   knum     - Integer, number of k-points along each direction
    %   stepmax  - Integer, maximum number of Hartree-Fock steps
    %   critial  - Float, convergence criterion for self-consistency
    %
    % Outputs:
    %   xinitial - Cell array, converged Hartree-Fock states
    %   ni       - Array, converged occupation numbers
    %   si       - Array, spin polarization values
    %   efermi   - Float, Fermi energy of the system

    % Generate 2D k-mesh
    [Kx, Ky, Kz] = gs.get_Bulk2Dkmesh([0, 1], [0, 1], knum);
    kpoints = [Kx(:), Ky(:), Kz(:)];

    % Set electric field to zero (default)
    Electric_field_in_evpA = 0;

    % Run the Hartree-Fock solver
    [xinitial, ni, si, efermi] = runhartreev8(gs, knum, Kx, Ky, Kz, kpoints, ...
        Electric_field_in_evpA, xinitial, stepmax, critial, U, V, 1e-10, 4);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Function to find the n-th NN neighbor pairing          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result_matrix = find_neighbor_data(coords_cartesian, lattice_vectors, n, dimensionality)
    % 动态处理二维或三维晶格，计算第 n 近邻的原子对信息
    %
    % 输入:
    % coords_cartesian: 原胞内的笛卡尔坐标 (n_atoms x 3)
    % lattice_vectors: 晶格基矢量 (3 x 3)
    % n: 第 n 近邻
    % dimensionality: 2 表示二维晶格，3 表示三维晶格
    %
    % 输出:
    % result_matrix: (num_pairs x 6) 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 1. 将笛卡尔坐标转换为分数坐标
    inv_lattice = inv(lattice_vectors);
    coords_fractional = coords_cartesian * inv_lattice;

    % 2. 构造超胞
    if dimensionality == 2
        % 仅扩展 x 和 y
        [super_coords, super_indices] = construct_supercell_2d(coords_fractional, lattice_vectors, 1);
    elseif dimensionality == 3
        % 扩展 x, y, z
        [super_coords, super_indices] = construct_supercell(coords_fractional, lattice_vectors, 1);
    else
        error('Dimensionality must be 2 or 3');
    end

    % 3. 在 Non-PBC 条件下计算超胞的距离矩阵
    if dimensionality == 2
        % 仅考虑 x-y 平面距离
        dist_matrix_non_pbc = compute_distance_matrix_2d(super_coords, lattice_vectors, false);
    elseif dimensionality == 3
        % 考虑完整三维距离
        dist_matrix_non_pbc = compute_distance_matrix(super_coords, lattice_vectors, false);
    end

    % 4. 提取 unique 距离（非零，考虑浮点误差）
    raw_distances = triu(dist_matrix_non_pbc);
    % raw_distances = raw_distances(raw_distances > 0); % 去掉 0 距离
    rounded_distances = round(raw_distances, 2); % 保留两位小数，分组距离
    unique_distances = unique(rounded_distances, 'sorted');

    % 5. 找到第 n 近邻的距离
    if n > length(unique_distances)
        error('第 %d 近邻超出最大可能的距离范围', n);
    end
    nth_distance = unique_distances(n);

    % 6. 找到满足第 n 近邻距离的原子对
    [pair_i, pair_j] = find(abs(dist_matrix_non_pbc - nth_distance) < 1e-2); % 容忍浮点误差

    % 7. 直接生成结果矩阵
    result_matrix=find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors);
end

function unique_pairs = find_unique_nth_neighbors(pair_i, pair_j, super_coords, super_indices, nth_distance, lattice_vectors)
    % 找到所有第 n 近邻的原子对并去除重复
    %
    % 输入:
    % pair_i, pair_j: 满足第 n 近邻条件的原子对索引
    % super_coords: 超胞中原子的分数坐标
    % super_indices: 超胞中原子的索引和周期性偏移
    % nth_distance: 第 n 近邻的距离
    % lattice_vectors: 晶格基矢量
    %
    % 输出:
    % unique_pairs: 矩阵，包含 [i, j, frax, fray, fraz, dis]

    % 初始化结果存储
    num_pairs = length(pair_i);
    all_pairs = zeros(num_pairs, 6); % [i, j, frax, fray, fraz, dis]

    % 遍历所有原子对
    for k = 1:num_pairs
        i = pair_i(k); % 原子 i
        j = pair_j(k); % 原子 j

        % 原胞内的原子编号
        atom_i = super_indices(i, 1);
        atom_j = super_indices(j, 1);

        % 计算周期性偏移（分数坐标差）
        delta_r = super_coords(i, :) - super_coords(j, :);

        % 将分数坐标差转换为笛卡尔坐标，用于计算实际距离
        delta_cartesian = delta_r * lattice_vectors;
        distance = sqrt(sum(delta_cartesian.^2));

        % 检查是否接近第 n 近邻距离
        if abs(distance - nth_distance) < 1e-1
            % 记录原子对信息
            all_pairs(k, :) = [atom_i, atom_j, super_indices(j,2:end)-super_indices(i,2:end), nth_distance];
        end
    end

    % 去除重复的原子对（如 i->j 和 j->i）
    % 按原子对的编号排序，并提取唯一值
    [~, unique_rows] = unique(all_pairs, 'rows');
    unique_pairs = all_pairs(unique_rows, :); % 提取唯一行
end

function dist_matrix = compute_distance_matrix(coords, lattice_vectors, pbc)
    % 计算三维距离矩阵，可选是否使用周期性边界条件 (PBC)
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            delta_cartesian = delta_r * lattice_vectors; % 转为笛卡尔坐标
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function dist_matrix = compute_distance_matrix_2d(coords, lattice_vectors, pbc)
    % 计算二维距离矩阵（仅考虑 x 和 y 方向），可选 PBC
    n = size(coords, 1); % 原子数
    dist_matrix = zeros(n, n); % 初始化距离矩阵

    for i = 1:n
        for j = i+1:n
            delta_r = coords(i, :) - coords(j, :);
            if pbc
                delta_r = delta_r - round(delta_r); % 最近镜像
            end
            % 仅保留 x 和 y 方向的笛卡尔坐标
            delta_cartesian = delta_r * lattice_vectors; 
            delta_cartesian = delta_cartesian(:, 1:2); % x 和 y 方向
            dist_matrix(i, j) = sqrt(sum(delta_cartesian.^2));
            dist_matrix(j, i) = dist_matrix(i, j); % 对称性
        end
    end
end

function [super_coords, super_indices] = construct_supercell(coords, lattice_vectors, scale)
    % 构造三维超胞
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            for k = -scale:scale
                offset = [i, j, k];
                super_coords = [super_coords; coords + offset];
                for a = 1:n_atoms
                    super_indices = [super_indices; a, i, j, k];
                end
            end
        end
    end
end

function [super_coords, super_indices] = construct_supercell_2d(coords, lattice_vectors, scale)
    % 构造二维超胞（仅扩展 x 和 y 方向）
    n_atoms = size(coords, 1);
    super_coords = [];
    super_indices = [];

    for i = -scale:scale
        for j = -scale:scale
            offset = [i, j, 0]; % z 方向固定为 0
            super_coords = [super_coords; coords + offset];
            for a = 1:n_atoms
                super_indices = [super_indices; a, i, j, 0];
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   Function to run hartree                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xinitial_U, xinitial_V] = initializeStates(U, V, pairsU, pairsV, nbands,xinitial_0)

    % initializeStates initializes the xinitial_U and xinitial_V arrays
    % Inputs:
    %   U       - Cell array or vector for U
    %   V       - Cell array or vector for V
    %   pairsU  - Cell array of pair matrices corresponding to U
    %   pairsV  - Cell array of pair matrices corresponding to V
    %   nbands  - Number of bands
    % Outputs:
    %   xinitial_U - Cell array containing initialized 3D matrices for U
    %   xinitial_V - Cell array containing initialized 3D matrices for V

    % Initialize empty cell arrays
    xinitial_U = {};
    xinitial_V = {};
    
    % Populate xinitial_U
    if ~isempty(U)
        for i = 1:length(U)
            pair = pairsU{i}; % Get the (i-1)th pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_idx = 1:num_pair
                % 提取原子对索引和分数坐标
                atom_i = pair(pair_idx, 1); % 原子 i
                atom_j = pair(pair_idx, 2);
                x=U(i);
                nj=xinitial_0(2*atom_j-1,2*atom_j-1)+xinitial_0(2*atom_j,2*atom_j);
                correlation(pair_idx,2*atom_i-1, 2*atom_i-1) = x*nj;%
                correlation(pair_idx,2*atom_i, 2*atom_i) = x*nj;%
            end
            xinitial_U = [xinitial_U, correlation];
        end
    end
    
    % Populate xinitial_V
    if ~isempty(V)
        for i = 1:length(V)
            pair = pairsV{i}; % Get the ith pair
            num_pair=size(pair,1);
            % 初始化临时矩阵，维度 (nbands, nbands)
            correlation = zeros(num_pair,nbands, nbands);
            for pair_index = 1:num_pair
                % 提取原子对索引  
                atom_i = pair(pair_index, 1); % 原子 i
                atom_j = pair(pair_index, 2); % 原子 j
                row_up_i = 2 * atom_i - 1;
                row_dn_i = 2 * atom_i;
                col_up_j = 2 * atom_j - 1;
                col_dn_j = 2 * atom_j;
                x=V(i)*0;
                correlation(pair_index,row_up_i, col_up_j) = -x; % 上自旋部分
                correlation(pair_index,row_dn_i, col_dn_j) = -x; % 下自旋部分
                correlation(pair_index,row_dn_i, col_up_j) = 0;
                correlation(pair_index,row_up_i, col_dn_j) = 0;
                correlation(pair_index,col_up_j, row_up_i) = -x; % 上自旋部分
                correlation(pair_index,col_dn_j, row_dn_i) = -x; % 下自旋部分
                correlation(pair_index,col_dn_j, row_up_i) = 0;
                correlation(pair_index,col_up_j, row_dn_i) = 0;

            end
            xinitial_V = [xinitial_V, correlation];
        end
    end
end

function V = dual_gate_potential(r, d, epsilon, n_max)
    % Calculate dual-gate Coulomb potential using image charge method
    % r: radial distance
    % d: distance to the gates
    % q: charge
    % epsilon: dielectric constant
    % n_max: number of image charges to consider
    e=1.602176634*10^(-19);
    epsilon_0=8.854187817*10^(-12);
    k_e=1/(4*pi*epsilon_0);
    V = 0; % Initialize potential
    for n = -n_max:n_max
        V = V + ((-1)^n) / sqrt(r^2 + (2*n*d)^2);
    end
    % V = e^2/(4*pi*epsilon_0)/e/epsilon * V; % Final potential meV*m
    V=e^2*k_e*1/epsilon/e * V; %to eV
end

function [xinitial,ni,si,efermi]=runhartreev8(gs,knum,Kx,Ky,Kz,kpoints,Electric_field_in_evpA,xinitial,stepmax,critial, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)

        [xinitial, metaData] = flattenNestedCell(xinitial);        % 展平操作
        objective = @(xinitial) one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV);
        xinitial = quasi_newton(objective, xinitial, 7, critial, stepmax);
        xinitial = restoreNestedCell(xinitial, metaData);  % 复原操作
        modifyHam(gs, xinitial, U, V, pairsU, pairsV)  %修改Ham
        nbands=size(gs.ham,1);
        [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(gs,Electric_field_in_evpA,Kx,Ky,Kz);
        Enk=reshape(Enk,[knum^2,nbands]);
        Unk=reshape(Unk,[nbands,nbands,knum^2]);
        [~,kindex,bandindex,efermi]=Total_energy(Enk,u2);
        [ni,si]=calonsite(Unk,kindex,bandindex);
end

function xnew = one_step_hf_v6(gs, knum, Kx, Ky, Kz, kpoints, Electric_field_in_evpA,xinitial,metaData, U0, U, V, u1, u2, pairsU0, pairsU, pairsV)
    % 1. 基本初始化
    gs.ham=gs.iniham+0; 
    nbands = size(gs.ham, 1);
    % 恢复操作
    xinitial = restoreNestedCell(xinitial, metaData);
    xinitial_0=xinitial{1};
    if ~isempty(U)
         xinitial_U=xinitial{2};
    end 
    if ~isempty(V)
        xinitial_V=xinitial{3};
    end 

    gs.onsite_modify(xinitial_0);

    if ~isempty(U)
        for i = 1:length(U)
            pairs = pairsU{i};
            xinitial_1 = xinitial_U{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_1(j,:,:)));
            end
        end
    end

    if ~isempty(V)
        for i = 1:length(V)
            pairs = pairsV{i};
            xinitial_2 = xinitial_V{i};
            for j = 1:size(pairs,1)
                gs.offsite_modify(pairs(j,3:5),squeeze(xinitial_2(j,:,:)));
            end
        end
    end


    
    % 2. 计算波函数和能量 (含电场影响)
    [Unk, Enk] = MTB.ham.get_bulk_plane_bands_add_electric(gs, Electric_field_in_evpA, Kx, Ky, Kz);
    Unk = reshape(Unk, [nbands, nbands, knum^2]);
    Enk = reshape(Enk, [knum^2, nbands]);
    C_k_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2);
    C_k = calculate_C_k_u1u2(Unk, Enk, u1, u2);

    % 3. 计算 U 和 V 的关联函数
    onsite_correlation_U =  zeros(size(xinitial_0));
    offsite_correlation_U = {};
    offsite_correlation_V = {};
    total_correlation = {};
   
    % 3.1 计算 U[0] 的关联函数 (onsite)
    if ~isempty(U0)
        for i = 1:length(U0)
            onsite_U = calculate_onsite_matrix(C_k_avg,pairsU0{i}, U0(i));
            % onsite_U = onsite_U*0.8+xinitial_0*0.2;
            onsite_correlation_U = onsite_correlation_U + onsite_U;
        end
        % To keep the TR symmetry
        % % diag_terms = diag(onsite_correlation_U)+diag(C_k_avg);
        % % onsite_correlation_U(1:nbands+1:end)=diag_terms/2;
        % % onsite_correlation_U=(onsite_correlation_U+onsite_correlation_U')/2;
    end
    
    % 3.2 计算 U[1:] 的关联函数 (offsite)
    if ~isempty(U)
        for i = 1:length(U)           
            offsite_U = calculate_offsite_U(C_k, kpoints, pairsU{i}, gs.a, U(i));
            % offsite_U = offsite_U*0.8 + xinitial{2}{i}*0.2;
            offsite_correlation_U = [offsite_correlation_U,offsite_U];
        end
    end
    
    % 3.3 计算 V 的关联函数 (offsite_V)
    if ~isempty(V)
        for i = 1:length(V)
            offsite_V = calculate_offsite_V(C_k, kpoints, pairsV{i}, gs.a, V(i));
            % offsite_V = offsite_V*0.8 + xinitial{3}{i}*0.2;         
            offsite_correlation_V = [offsite_correlation_V, offsite_V];
        end
    end
    total_correlation = {onsite_correlation_U,offsite_correlation_U,offsite_correlation_V};

    [xnew, ~] = flattenNestedCell(total_correlation);% 展平成列向量
end

function modifyHam(gs, xinitial, U, V, pairsU, pairsV)
    % modifyStates modifies onsite and offsite states using gs object
    % Inputs:
    %   gs        - Object containing methods `onsite_modify` and `offsite_modify`
    %   xinitial  - Cell array containing initial state matrices
    %   U         - Cell array or vector for U
    %   V         - Cell array or vector for V
    %   pairsU    - Cell array of pair matrices corresponding to U
    %   pairsV    - Cell array of pair matrices corresponding to V

    % Handle onsite modification for xinitial{1}
    xinitial_0 = xinitial{1};
    gs.onsite_modify(xinitial_0);

    % Handle offsite modification for U
    if ~isempty(U)
        xinitial_U = xinitial{2}; % Extract second cell for U states
        for i = 1:length(U)
            pairs = pairsU{i};  % Get pairs for the (i-1)th entry
            xinitial_1 = xinitial_U{i}; % Get initial state for U
            for j = 1:size(pairs, 1)
                % Modify offsite states for U
                gs.offsite_modify(pairs(j, 3:5), squeeze(xinitial_1(j, :, :)));
            end
        end
    end

    % Handle offsite modification for V
    if ~isempty(V)
        xinitial_V = xinitial{3}; % Extract third cell for V states
        for i = 1:length(V)
            pairs = pairsV{i};    % Get pairs for the ith entry
            xinitial_2 = xinitial_V{i}; % Get initial state for V
            for j = 1:size(pairs, 1)
                % Modify offsite states for V
                gs.offsite_modify(pairs(j, 3:5), squeeze(xinitial_2(j, :, :)));
            end
        end
    end
end

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
    % % a=a(ceil(size(Enk(:),1)*0.5)+1:end);
    % % b=b(ceil(size(Enk(:),1)*0.5)+1:end);

    % end
    %find index in Enk
    row=mod(b,knum);row(row==0)=knum;
    col=ceil(b./knum);
    Etot=sum(a,'all')/knum;
    kindex=row;
    bandindex=col;
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end

function C_total_avg = calculate_correlation_with_range_avg(Enk, Unk, u1, u2)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)，每行是一个 k 点的能带
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    % u1: 填充比例下限 (0 <= u1 <= 1)
    % u2: 填充比例上限 (0 <= u2 <= 1)
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)

    % 获取维度信息
    [knum2, nbands] = size(Enk);  % Enk 的维度
    C_total = zeros(nbands, nbands); % 初始化总关联函数矩阵

    % 1. 计算 u1 和 u2 对应的费米能级
    E_low = calculate_ef(Enk, u1); % u1 对应的费米能级
    E_high = calculate_ef(Enk, u2); % u2 对应的费米能级

    %fprintf('能级范围: %.6f eV 到 %.6f eV\n', E_low, E_high);

    % 2. 遍历每个 k 点，累加所有 C^(m)
    for m = 1:knum2
        % 提取第 m 个 k 点的本征值和本征矢量
        E_k = Enk(m, :);          % 第 m 个 k 点的能带本征值，维度 (1, nbands)
        U_k = Unk(:, :, m);       % 第 m 个 k 点的本征矢量矩阵，维度 (nbands, nbands)

        % 筛选能级：处于 [E_low, E_high] 区间的态置为 1，其余置为 0
        W_k = diag(double(E_k >= E_low & E_k <= E_high)); % 维度 (nbands, nbands)

        % 计算关联函数矩阵 C^(m)
        C_m = conj(U_k) * W_k * transpose(U_k);

        % 累加当前 k 点的关联函数矩阵
        C_total = C_total + C_m;
    end

    % 3. 求平均：除以总的 k 点数目 knum^2
    C_total_avg = C_total / knum2;

    %fprintf('关联函数计算完成，结果已对 k 点求平均。\n');
end

function C_k = calculate_C_k_u1u2(U_k, E_k, u1, u2)
    % 计算关联函数矩阵 C_k，使用 u1 和 u2 提取占据能级范围
    %
    % 输入:
    % U_k: 本征矢量矩阵，维度 (nbands, nbands, nk)
    % E_k: 本征值矩阵，维度 (nk, nbands)
    % u1, u2: 占据比例范围
    %
    % 输出:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)

    % 获取维度
    [nbands, ~, nk] = size(U_k);
    C_k = zeros(nbands, nbands, nk);

    % 确定占据能级范围 [E_low, E_high]
    total_states = numel(E_k);
    E_sorted = sort(E_k(:)); % 将所有能级排序
    num_u1 = ceil(total_states * u1);
    num_u2 = ceil(total_states * u2);
    E_low = E_sorted(num_u1);
    E_high = E_sorted(num_u2);

    % 遍历每个 k 点
    for k = 1:nk
        % 提取当前 k 点的本征值和本征矢量
        E_k_point = E_k(k, :);
        U_k_point = U_k(:, :, k);

        % 构建占据函数 f(E_k) 范围
        f_k = double(E_k_point >= E_low & E_k_point <= E_high);

        % 计算关联函数矩阵 C_k
        C_k(:, :, k) = conj(U_k_point) * diag(f_k) * transpose(U_k_point);
    end
end

function onsite_U = calculate_onsite_matrix(correlation, pairs, U)

    % 计算 Onsite 矩阵，仅保留对角线上的 2x2 block 并进行变换
    %
    % 输入:
    % correlation: 关联函数矩阵，维度 (2*nbands, 2*nbands)
    % U: Hubbard U 参数
    %
    % 输出:
    % onsite_matrix: Onsite 矩阵，维度 (2*nbands, 2*nbands)

    % 获取轨道数目
    total_bands = size(correlation, 1);
    if mod(total_bands, 2) ~= 0
        error('关联函数矩阵的行数必须是偶数，每个轨道对应上下自旋');
    end
    nbands = total_bands / 2;

    % 定义变换矩阵 T
    T = [0, -1; 1, 0];

    % 初始化输出矩阵
    onsite_matrix = zeros(total_bands, total_bands);
    num_pairs = size(pairs,1);

    % 遍历每个对角线的 2x2 block
       % 遍历每个原子对
    for pair_idx = 1:num_pairs
        i = pairs(pair_idx,1);
        % 提取对角线上的 2x2 block
        block = correlation(2*i-1:2*i, 2*i-1:2*i);

        % 应用变换矩阵 T
        transformed_block = T * block * T';

        % 保存到输出矩阵的对应位置
        onsite_matrix(2*i-1:2*i, 2*i-1:2*i) = U * transformed_block;
        % To keep the TR symmetry
    end
    onsite_U=onsite_matrix;
end

function offsite_U = calculate_offsite_U(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_U = zeros(num_pairs, nbands, nbands);
    % offsite_U_total = zeros(nbands,nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义行/列索引
        index_i = 2 * atom_i - 1; 
        index_j = 2 * atom_j - 1; 


        % 遍历所有 k 点，计算相位修正的关联值
        ni_up_sum = 0;
        ni_dn_sum = 0;
        nj_up_sum = 0;
        nj_dn_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

            % 提取 C_k 中的元素
            % ni_up = C_k(index_i, index_i, k) * phase_factor; % 上自旋
            % ni_dn = C_k(index_i+1, index_i+1, k) * phase_factor; % 下自旋
            nj_up = C_k(index_j, index_j, k) * phase_factor; % 上自旋
            nj_dn = C_k(index_j+1, index_j+1, k) * phase_factor; % 下自旋
            % 累加相位修正的贡献
            % ni_up_sum = ni_up_sum + ni_up;
            % ni_dn_sum = ni_dn_sum + ni_dn;
            nj_up_sum = nj_up_sum + nj_up;
            nj_dn_sum = nj_dn_sum + nj_dn;
        end

        % 平均化，并乘以 U
        % ni_up_sum_avg = real(ni_up_sum/nk) * U;
        % ni_dn_sum_avg = real(ni_dn_sum/nk) * U;
        nj_up_sum_avg = real(nj_up_sum/nk) * U;
        nj_dn_sum_avg = real(nj_dn_sum/nk) * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(index_i, index_i) = nj_dn_sum_avg+nj_up_sum_avg; %
        correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg+nj_dn_sum_avg; %
        % correlation_temp(index_j, index_j) = ni_dn_sum_avg+ni_up_sum_avg; %
        % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg+ni_dn_sum_avg; %
        % To keep the PT symmetry
        % % % correlation_temp(index_i, index_i) = nj_dn_sum_avg/2+nj_up_sum_avg/2; %
        % % % correlation_temp(index_i+1, index_i+1) = nj_up_sum_avg/2+nj_dn_sum_avg/2; %
        % % % correlation_temp(index_j, index_j) = ni_dn_sum_avg/2+ni_up_sum_avg/2; %
        % % % correlation_temp(index_j+1, index_j+1) = ni_up_sum_avg/2+ni_dn_sum_avg/2; %       
        % 将当前原子对的矩阵存储到输出列表中
        offsite_U(pair_idx, :, :) = correlation_temp;
        % offsite_U_total = offsite_U_total+correlation_temp;
    end
end

function offsite_V = calculate_offsite_V(C_k, kpoints, result_matrix, lattice_vectors, U)
    % 计算 Offsite Hubbard U 矩阵，输出维度为 (num_pairs, nbands, nbands)
    %
    % 输入:
    % C_k: 关联函数矩阵，维度 (nbands, nbands, nk)
    % kpoints: k 点坐标，维度 (nk, 3)
    % result_matrix: 最近邻原子对信息，包含 [i, j, R_fractional]
    % lattice_vectors: 晶格基矢量，维度 (3, 3)
    % U: Hubbard U 参数
    %
    % 输出:
    % offsite_U_list: 每个原子对的 Offsite U 矩阵，维度 (num_pairs, nbands, nbands)

    % 获取输入维度
    [nbands, ~, nk] = size(C_k);
    num_pairs = size(result_matrix, 1);

    % 初始化输出矩阵，维度为 (num_pairs, nbands, nbands)
    offsite_V = zeros(num_pairs, nbands, nbands);

    % 遍历每个原子对
    for pair_idx = 1:num_pairs
        % 提取原子对索引和分数坐标
        atom_i = result_matrix(pair_idx, 1); % 原子 i
        atom_j = result_matrix(pair_idx, 2); % 原子 j
        R_fractional = result_matrix(pair_idx, 3:5); % 原子对的分数坐标

        % 将分数坐标转换为笛卡尔坐标
        R_cartesian = R_fractional * lattice_vectors;

        % 初始化临时矩阵，维度 (nbands, nbands)
        correlation_temp = zeros(nbands, nbands);

        % 定义上下自旋的行/列索引
        col_dn_i = 2 * atom_i;       % 下自旋 (dn) 行
        col_up_i = 2 * atom_i - 1;   % 上自旋 (up) 行
        row_dn_j = 2 * atom_j;       % 下自旋 (dn) 列
        row_up_j = 2 * atom_j - 1;   % 上自旋 (up) 列

        % 遍历所有 k 点，计算相位修正的关联值
        dndn_sum = 0;
        upup_sum = 0;
        updn_sum = 0;
        dnup_sum = 0;
        for k = 1:nk
            % 计算相位因子 e^{i*k*R}
            phase_factor = exp(1i * dot(kpoints(k, :), R_cartesian));

            % 提取 C_k 中的元素
            dndn = C_k(row_dn_j, col_dn_i, k) * phase_factor; % 下自旋-下自旋
            upup = C_k(row_up_j, col_up_i, k) * phase_factor; % 上自旋-上自旋
            updn = C_k(row_up_j, col_dn_i, k) * phase_factor;
            dnup = C_k(row_dn_j, col_up_i, k) * phase_factor;

            % 累加相位修正的贡献
            dndn_sum = dndn_sum + dndn;
            upup_sum = upup_sum + upup;
            updn_sum = updn_sum + updn;
            dnup_sum = dnup_sum + dnup;
        end

        % 平均化，并乘以 U
        dndn_avg = dndn_sum / nk * U;
        upup_avg = upup_sum / nk * U;
        updn_avg = updn_sum / nk * U;
        dnup_avg = dnup_sum / nk * U;
        % 将 dndn 和 upup 填充到临时矩阵
        correlation_temp(col_up_i, row_up_j) = -upup_avg; % 上自旋部分
        correlation_temp(col_dn_i, row_dn_j) = -dndn_avg; % 下自旋部分
        correlation_temp(col_dn_i, row_up_j) = -updn_avg; % 上自旋部分
        correlation_temp(col_up_i, row_dn_j) = -dnup_avg; % 下自旋部分

        % 将当前原子对的矩阵存储到输出列表中
        offsite_V(pair_idx, :, :) = correlation_temp;
    end
end

function projector = get_projector(gs, kpoints,Unk, whichbands)
    % 计算费米能级附近填充比例 u1 和 u2 对应区间的关联函数矩阵，并对 knum^2 求平均
    %
    % 输入:
    % gs: geometry
    % kpoints: 面内mesh 维度(knum^2,3)
    % whichbands: 选择投影的band 维度(1,nbands)
    % Unk: 本征矢量矩阵，维度 (nbands, nbands, knum^2)，每层对应一个 k 点的本征矢量
    %
    %
    % 输出:
    % C_total_avg: 平均关联函数矩阵，维度 (nbands, nbands)
    [nbands,~,knum2]=size(Unk);
    projectors_k = zeros(nbands,nbands,knum2);
    % 构造 projectors_k
    for i = 1:knum2
        projectors_k(:, :, i) = diag(which_bands); % 按对角阵填充
    end

    % 应用投影算符
    vs_conj = conj(Unk);
    vs_transpose = permute(Unk, [1, 2, 1]); % 转置第一和第二维
    projectors_k = pagemtimes(pagemtimes(vs_conj, projectors_k), vs_transpose);


    rvectors=gs.hopr;
    expr=exp(1j*(kpoints*(rvectors*gs.a).'));
    % 获取维度信息
    [knum2, nbands] = size(Enk);  % Enk 的维度

    
    C_total = zeros(nbands, nbands); % 初始化总关联函数矩阵

    % 1. 计算 u1 和 u2 对应的费米能级
    E_low = calculate_ef(Enk, u1); % u1 对应的费米能级
    E_high = calculate_ef(Enk, u2); % u2 对应的费米能级

    %fprintf('能级范围: %.6f eV 到 %.6f eV\n', E_low, E_high);

    % 2. 遍历每个 k 点，累加所有 C^(m)
    for m = 1:knum2
        % 提取第 m 个 k 点的本征值和本征矢量
        E_k = Enk(m, :);          % 第 m 个 k 点的能带本征值，维度 (1, nbands)
        U_k = Unk(:, :, m);       % 第 m 个 k 点的本征矢量矩阵，维度 (nbands, nbands)

        % 筛选能级：处于 [E_low, E_high] 区间的态置为 1，其余置为 0
        W_k = diag(double(E_k >= E_low & E_k <= E_high)); % 维度 (nbands, nbands)

        % 计算关联函数矩阵 C^(m)
        C_m = conj(U_k) * W_k * transpose(U_k);

        % 累加当前 k 点的关联函数矩阵
        C_total = C_total + C_m;
    end

    % 3. 求平均：除以总的 k 点数目 knum^2
    C_total_avg = C_total / knum2;

    %fprintf('关联函数计算完成，结果已对 k 点求平均。\n');
end


function [band_gap, VBM, CBM, VBM_k, CBM_k] = compute_band_gap_1d(En, valence_band_index, conduction_band_index)
% compute_band_gap_1d: 计算指定价带和导带之间的能隙及对应的 k 点位置
%
% 输入参数:
%   En - 能量矩阵，维度为 (knum, nbands)，包含能量数据
%   valence_band_index - 价带索引 (整数)
%   conduction_band_index - 导带索引 (整数)
%
% 输出参数:
%   band_gap - 带隙 (CBM - VBM)
%   VBM - 价带最高点能量值
%   CBM - 导带最低点能量值
%   VBM_k - 价带最高点对应的 k 点位置 (整数)
%   CBM_k - 导带最低点对应的 k 点位置 (整数)

    % 获取输入矩阵的维度
    [knum, nbands] = size(En);
    assert(valence_band_index <= nbands && conduction_band_index <= nbands, ...
        '价带或导带索引超出能量矩阵的范围');

    % Step 1: 计算 VBM 和 CBM 及其索引
    [VBM, VBM_k] = max(En(:, valence_band_index)); % 价带最高点
    [CBM, CBM_k] = min(En(:, conduction_band_index)); % 导带最低点

    % Step 2: 计算带隙
    band_gap = CBM - VBM;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Function for scf convergency                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xn = quasi_newton(func, x0, q, tol, maxsteps, history)
    % 默认参数设置
    if nargin < 3, q = 10; end
    if nargin < 4, tol = 1e-7; end
    if nargin < 5, maxsteps = 100; end
    if nargin < 6, history = {}; end

    % 初始化变量
    f0 = func(x0).'; % 初始梯度，确保为列向量
    n = length(f0); % 确保 x0 和 func 输出长度一致

    alpha = 1e-3; % 初始步长
    U = zeros(n, q); % 初始化 U 矩阵
    V = zeros(n, q); % 初始化 V 矩阵

    xn = x0(:); % 确保 x0 是列向量
    fn = f0; % 初始化梯度为列向量

    disp('Starting quasi-Newton...');

    % === 预热阶段：初始化 U 和 V ===
    for i = 1:q
        fn = func(xn).'; % 确保梯度为列向量
        xn_1 = xn - alpha * fn; % 初步更新变量

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 更新 U 和 V 矩阵
        U(:, i) = fn - xn; % 梯度差
        V(:, i) = func(fn).' - fn; % 梯度映射差
        xn = xn_1; % 更新变量
    end

    % % % % % % disp('Warmed up.');

    % === 主迭代循环 ===
    for i = 1:maxsteps
        fn = func(xn).'; % 确保梯度为列向量

        % 检查维度一致性
        xn = xn(:); % 确保 xn 是列向量
        fn = fn(:); % 确保 fn 是列向量
        if size(U, 1) ~= length(xn)
            error('Dimension mismatch: U should have %d rows, but got %d.', length(xn), size(U, 1));
        end

        % 构造矩阵 C 和向量 b
        C = transpose(U) * U - transpose(U) * V; % 计算 C 矩阵
        b = transpose(U) * (xn - fn); % 计算 b 向量

        % 判断 C 是否病态，必要时正则化
        if rcond(C) < 1e-10
            C = C + 1e-8 * eye(size(C));
        end

        % 求解更新方向 delta
        delta = V * (C \ b);

        % 更新变量
        xn_1 = fn - delta;

        % 对变量进行非负修正
        xn_1 = max(0.0, xn_1);

        % 存储历史记录
        history{end+1} = xn_1;

        % 更新 U 和 V 矩阵
        fn_1 = func(xn_1).'; % 确保梯度为列向量
        U = [U(:, 2:end), fn_1 - xn_1]; % 滑动窗口更新 U
        V = [V(:, 2:end), func(fn_1).' - fn_1]; % 滑动窗口更新 V

        % 检查收敛条件
        dx = max(abs(xn_1 - xn)); % 计算变量更新幅度
        fprintf('Iteration %d, dx = %e\n', i, dx);

        if dx < tol
            fprintf('Converged at iteration %d\n', i);
            xn = xn_1.'; % 最终结果转置为行向量
            return;
        end

        xn = xn_1; % 更新当前变量
    end

    % 未收敛提示
    disp('Did not converge within the maximum steps.');
    xn = xn.'; % 如果未收敛，也转置为行向量
end

function [flattenedVector, metaData] = flattenNestedCell(C)
    flattenedVector = []; % 用于存储展平后的行向量
    metaData = {};        % 用于存储路径和尺寸信息

    % 递归处理单元数组
    function processElement(element, path)
        if iscell(element)
            % 如果是单元数组，递归处理每个子元素
            for i = 1:numel(element)
                processElement(element{i}, [path, i]);
            end
        else
            % 如果是普通数组，展平并记录信息
            flattenedVector = [flattenedVector, element(:)']; % 展平并拼接
            metaData{end+1} = struct('Path', {path}, 'Size', size(element)); % 记录路径和尺寸
        end
    end

    % 开始递归处理
    processElement(C, []);
end

function restoredCell = restoreNestedCell(flattenedVector, metaData)
    restoredCell = {};    % 初始化恢复后的单元数组
    currentIndex = 1;     % 当前展平向量的索引

    % 递归赋值函数
    function target = assignElement(target, path, sizeInfo)
        if isempty(path)
            % 如果路径为空，说明已经到达最终节点
            numElements = prod(sizeInfo); % 当前数组的总元素数
            reshapedArray = reshape(flattenedVector(currentIndex:currentIndex + numElements - 1), sizeInfo); % 恢复形状
            currentIndex = currentIndex + numElements; % 更新索引位置
            target = reshapedArray; % 直接返回恢复后的数组
        else
            % 处理嵌套单元数组
            idx = path(1);
            if numel(target) < idx || ~iscell(target{idx})
                target{idx} = {}; % 确保路径上的每一级是单元数组
            end
            target{idx} = assignElement(target{idx}, path(2:end), sizeInfo); % 递归处理下一层路径
        end
    end

    % 遍历元数据并恢复每个元素
    for i = 1:numel(metaData)
        restoredCell = assignElement(restoredCell, metaData{i}.Path, metaData{i}.Size);
    end
end


function strain_band_plot(material, data_path, s_values, Electric_field_in_evpA, nk, efermi, labels)
    % material: 材料名称 (e.g., "TaIrTe4")
    % data_path: 紧束缚模型数据路径 (e.g., "data/TaIrTe4_2d_tb/1019")
    % s_values: 应变比例的数组 (e.g., 0.85:0.01:1.05)
    % Electric_field_in_evpA: 施加的电场 (单位 eV/Å)
    % nk: k 点数量
    % efermi: 费米能级
    % labels: k 路径标签

    % 初始化几何结构
    g2 = initializeGeometry(material, data_path);
    
    % 基底变换
    T = getBasisTransformMatrix();
    g2.ham = transformBasis(g2.ham, T);

    % 记录原始的晶格常数
    a_original = g2.a;

    % 遍历应变参数 s
    for s = s_values
        % 施加应变
        strain = [1/s, 0, 0; 0, s, 0; 0, 0, 1];
        g2.a = strain * a_original;  % 修改晶格常数
        g2.wpos = setWannierPosition(g2);

        % 计算邻近矩阵
        [result_matrices2, pairsU0, pairsU, pairsV] = findNeighbors(g2);

        % 计算修正后的跳跃积分
        for i = 2:size(result_matrices2,2)-2
            delta_ij = result_matrices2{i}(:,6) - result_matrices2{i}(:,6); % 确保 result_matrices2{i} 正确
            t_ij = exp(-3 * delta_ij ./ result_matrices2{i}(:,6)); % 这里 beta = 3, 可调

            % 添加修正因子到矩阵
            result_matrices2{i} = [result_matrices2{i}, t_ij];

            % 遍历所有邻接项
            for j = 1:length(result_matrices2{i})
                bandindex_i = result_matrices2{i}(j,1);
                bandindex_j = result_matrices2{i}(j,2);
                scale = result_matrices2{i}(j,7);
                
                % 找到匹配的跳跃项
                index = find(ismember(g2.hopr, result_matrices2{i}(j,3:5), "rows"));
                if ~isempty(index)
                    g2.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) = ...
                        g2.ham(2*bandindex_i-1:2*bandindex_i, 2*bandindex_j-1:2*bandindex_j, index) * scale;
                end
            end
        end

        % 记录初始哈密顿量
        g2.iniham = g2.ham + 0;

        % 计算能带
        [nbands,~,nrpts] = size(g2.ham);
        hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
        [Energy, kpath, kindex] = MTB.ham.get_bulk_bands_add_electric(...
            g2.ham, g2.hopr, g2.wpos, Electric_field_in_evpA, nbands, nrpts, hkpoints, nk, g2.a, g2.b);
        
        % 绘制能带
        hold on;
        MTB.plot.plot_bands(Energy, nbands, efermi, kpath, labels, kindex, ...
            strcat("TaIrTe4-2d (s=", num2str(s), ")"), Electric_field_in_evpA * 10000);
    end

    hold off;
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

%% Enk(nk,nk,nband) Kx(nk,nk) Ky(nk,nk)
function writeEnk(Enk,Kx,Ky,Occ,filename)
    file=fopen(filename, 'w');
    fprintf(file, '%s\n', '# kx     ky     kz     Ev5     Ev4     Ev3     Ev2     Ev1     Ec1     Ec2     Ec3     Ec4     Ec5');
    for i=1:size(Enk,1)
        for j=1:size(Enk,2)
            fprintf(file, [repmat('%12.6f',1,13) '\n'],...
                Kx(i,j), Ky(i,j), 0.0, ...
                Enk(i,j,Occ-4), Enk(i,j,Occ-3), Enk(i,j,Occ-2), Enk(i,j,Occ-1), Enk(i,j,Occ), ...
                Enk(i,j,Occ+1), Enk(i,j,Occ+2), Enk(i,j,Occ+3), Enk(i,j,Occ+4), Enk(i,j,Occ+5));
        end
        fprintf(file,'\n');
    end
    fclose(file);
end