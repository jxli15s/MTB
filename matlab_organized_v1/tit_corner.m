%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% %[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

% g = MTB.read_poscar(g,"data/TaIrTe4_2d/obs/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/obs/wannier90_hr_p1.dat','data/TaIrTe4_2d/obs/wannier90_hr_p2.dat');
dataPath='data/TaIrTe4_2d/obs/';
    g = MTB.read_poscar(g, fullfile(dataPath, "POSCAR"));
    % Read Wannier Hamiltonian data
    [g.ham, g.hopr] = MTB.wannier.read_hr(...
        fullfile(dataPath, "wannier90_hr_p1.dat"), ...
        fullfile(dataPath, "wannier90_hr_p2.dat"));
    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Construct Hamiltonian and Basis Transform              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
% g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
% %[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

% g = MTB.read_poscar(g,"data/TaIrTe4_2d/obs/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d/obs/wannier90_hr_p1.dat','data/TaIrTe4_2d/obs/wannier90_hr_p2.dat');
dataPath='data/tit_hf/15x1/10nmd/U-epsilon-f/9-15x1/wtool/r11/';
    g = MTB.read_poscar(g, fullfile(dataPath, "POSCAR"));
    % Read Wannier Hamiltonian data
    [g.ham, g.hopr] = MTB.wannier.read_hr(...
        fullfile(dataPath, "wannier90_hr_p1.dat"), ...
        fullfile(dataPath, "wannier90_hr_p2.dat"));
%%
clc;
clear;
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4_2d_tb/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/wannier90_hr_p2.dat');
%[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/1019/wannier90_hr_p1.dat','data/TaIrTe4_2d_tb/1019/wannier90_hr_p2.dat');

%%
[nbands,~,nrpts]=size(g.ham);
g.wpos=g.atoms*g.a+0;
labels={'Y','\Gamma','X','R'}; % labels for k
hkpoints={[0.0,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points

efermi=5.8681;
% efermi=0.0;
Electric_field_in_evpA=0.00*0.529177; nk=31;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);
hold on;
% plot(kpath,Energy(61,:)-efermi,'Color','red','LineWidth',2);
% plot(kpath,Energy(62,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')
%%
knum=101;
band1=1;
band2=60+2;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);
%%

%%
g.wpos=g.atoms*g.a+0;

% Step 4: Construct Supercell
n1 = 5; n2 =20; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

gs.iniham = gs.ham + 0; % Store the initial Ham

 plot(gs.wpos(:,1),gs.wpos(:,2),'ro')
%%
% vals=sort(eig(gs.ham(:,:,5)))
[V,D]=eig(gs.ham(:,:,5)+gs.ham(:,:,5)');
[E,ind]=sort(diag(D));
Psik=V(:,ind);
%%
figure()
plot(E,'ro')
figure()
plot(gs.wpos(1:120:end,1),gs.wpos(1:120:end,2),'ro')
%%
T=kron(diag(ones(1,n1*n2)),ones(1,120))
%%
x=gs.wpos(1:120:end,1);
y=gs.wpos(1:120:end,2);
% p=zeros(size(Psik));
% num=n1*n2*62+1;
% num=1521
% figure()
% scatter(x(:), y(:), 100, abs(Psik(:,num)).^2, 'filled');
% p(1:2:end,num)=Psik(1:2:end,num)+Psik(2:2:end,num);
% p(2:2:end,num)=Psik(1:2:end,num)+Psik(2:2:end,num);
pp=abs(Psik).^2;
p=T*pp;
%%
for i=5900:6300
    fig = figure('Visible', 'off');
    % p(1:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i);
    % scatter(x(:), y(:), 100, abs(p(:,i)).^2, 'filled');
    scatter(x(:), y(:), 100, p(:,i), 'filled');
    wxname="data/TaIrTe4_2d_tb/states/"+int2str(i)+'.png';
    saveas(fig, wxname)
end
%%
% scatter3(x, y, z, sizes, colors, 'filled');
figure()

% scatter3(x, y, p(:,6153), 100, p(:,6153), 'filled');
% caxis([0.05, 0.08]);
% view(45, 30);

scatter3(x, y, p(:,6225), 100, p(:,6225), 'filled');
caxis([0.031, 0.032]);
ylim([38,222])
view(0, 90);
%%

%%
%%
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
hkpoints={[0.5,0.5,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.5,0.0]};% hkpoints-high symmetry k points
efermi=-0.4423;
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

Electric_field_in_evpA=0.00;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"TaIrTe4-2d",Electric_field_in_evpA*10000);

hold on;
plot(kpath,Energy(88,:)-efermi,'Color','red','LineWidth',2);
plot(kpath,Energy(89,:)-efermi,"Color",'blue','LineWidth',2);
% print(bandname,'-dpng','-r600')
%%
% Step 4: Construct Supercell
n1 = 10; n2 =10; % Supercell dimensions
gs = constructSupercell(g, n1, n2);

gs.iniham = gs.ham + 0; % Store the initial Ham

 plot(gs.wpos(:,1),gs.wpos(:,2),'ro')
 %%
 %%
% vals=sort(eig(gs.ham(:,:,5)))
[V,D]=eig(gs.ham(:,:,5)+gs.ham(:,:,5)');
[E,ind]=sort(diag(D));
Psik=V(:,ind);
%%
figure();
plot(E,'ro');
figure();
plot(gs.wpos(1:124:end,1),gs.wpos(1:124:end,2),'ro');
%%
T=kron(diag(ones(1,n1*n2)),ones(1,124));
%%
x=gs.wpos(1:124:end,1);
y=gs.wpos(1:124:end,2);
% p=zeros(size(Psik));
% num=n1*n2*62+1;
% num=1521
% figure()
% scatter(x(:), y(:), 100, abs(Psik(:,num)).^2, 'filled');
% p(1:2:end,num)=Psik(1:2:end,num)+Psik(2:2:end,num);
% p(2:2:end,num)=Psik(1:2:end,num)+Psik(2:2:end,num);
pp=abs(Psik).^2;
p=T*pp;
%%
for i=8700:8900
    fig = figure('Visible', 'off');
    % p(1:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i)+Psik(2:2:end,i);
    % p(2:2:end,i)=Psik(1:2:end,i);
    % scatter(x(:), y(:), 100, abs(p(:,i)).^2, 'filled');
    scatter(x(:), y(:), 100, p(:,i), 'filled');
    wxname="data/TaIrTe4_2d_tb/states/"+int2str(i)+'.png';
    saveas(fig, wxname)
end
%%
figure()
% scatter3(x, y, sum(p(:, 8700:8850), 2), 100, sum(p(:, 8700:8850), 2), 'filled');
% caxis([2.38, 2.4]);
scatter3(x, y, sum(p(:, 8795:8798), 2), 100, sum(p(:, 8795:8798), 2), 'filled');
caxis([0.20, 0.25]);
view(45, 30);
%%
knum=51;
band1=85;
band2=88;
[wx,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);


%%
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