clear;
clear all;
% parpool('local',4)

g = MTB.geometry("WS2");
g = MTB.read_poscar(g,"data/WS2/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/WS2/wannier90_hr_p1.dat","data/WS2/wannier90_hr_p2.dat");
%%
MillerIndices=[1,-1,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
nk=31;
nslab=100;
delta=0.03;% 0.03~0.05
numEigs=30;
mus=linspace(6.4464,7.0464,51);%mu for E_f-mu to E_f+mu of 100 points
kpoint=[0.0,0.0];% Gamma Point
tic;
Energy=MTB.ham.get_slab_mu_E_sparse_BdG(g.ham,g.hopr2,nslab,nbands,numEigs,nrpts,kpoint,g.a2,mus,delta);
save("Energy-mus.mat","Energy");
toc
%%

%%
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
plot(x',Energy','o')
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


%%
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

