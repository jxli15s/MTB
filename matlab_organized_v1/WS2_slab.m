clear;
clear all;
g = MTB.geometry("WS2");
g = MTB.read_poscar(g,"data/WS2/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr("data/WS2/wannier90_hr_p1.dat","data/WS2/wannier90_hr_p2.dat");
%%
MillerIndices=[1,0,0];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%%
labels={'C','\Gamma','C'}; % labels for k
hkpoints={[0.5,0.0000000000],...
          [0.0000000000,0.0000000000],...
          [0.5,0.0000000000]};% hkpoints-high symmetry k points
nk=51;
nslab=50;
%%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=6.7464; %% set Fermi Level
%%
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")
%%
g.wpos=[]
g.wpos=g.atoms*g.a
orbital_num=[10,10,6,6,6,6]
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
        repmat(g.wpos(3,:),[orbital_num(3),1]);...
        repmat(g.wpos(4,:),[orbital_num(4),1]);...
        repmat(g.wpos(5,:),[orbital_num(5),1]);...
        repmat(g.wpos(6,:),[orbital_num(6),1]);...
    ]

%%
wpos=[]
for i=1:nslab;
    wpos=[wpos;g.wpos+g.a(3,:)*(i-1)]
end




 
%% plot bulk bands
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
hkpoints={[-0.2685588623,0.7314411377,0.0000000000],...
          [0.0000000000,0.0000000000,0.0000000000],...
          [-0.2685588623,0.7314411377,0.0000000000]};% hkpoints-high symmetry k points
nk=101;%number of k points
%Energy(nbands,nrpts) kpath((nk-1)*(length(hkpoints)-1)+1) kindex(length(hkpoints))
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
efermi=6.7464; %% set Fermi Level
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex)
%%

%% plot slab bands
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'C','\Gamma','C'}; % labels for k
hkpoints={[-0.2685588623,0.7314411377,0.0000000000],...
          [0.0000000000,0.0000000000,0.0000000000],...
          [-0.2685588623,0.7314411377,0.0000000000]};% hkpoints-high symmetry k points
nk=101;%number of k points
%Energy(nbands,nrpts) kpath((nk-1)*(length(hkpoints)-1)+1) kindex(length(hkpoints))
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
efermi=6.7464; %% set Fermi Level
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex)
%%


%%
Np=1;
omegamin=-0.3;
omegamax=0.1;
omeganum=501;
omegas=linspace(omegamin,omegamax,omeganum);
nk=301;
[dos_l,dos_r,dos_bulk,kindex,kpath]=MTB.ham.get_surfstates(g.ham,g.hopr2,nbands,nrpts,hkpoints,nk,Np,g.a2,g.b2,omegamax,omegamin,omeganum);

%%
figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas);
pcolor(Kx,Ky,dos_l);
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_r)
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)

figure('Color','White')
[Kx,Ky]=meshgrid(kpath,omegas)
pcolor(Kx,Ky,dos_bulk)
colormap(slanCM('vik'))
shading interp
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xticks(kindex)
xticklabels(labels)
ylabel('E-E_f (eV)','FontSize',20)
