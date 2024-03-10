clear;
clear all;
g = MTB.geometry("MnBiTe");
g = MTB.read_poscar(g,"data/MnBiTe-xu/POSCAR-MnBiTe");
[g.ham,g.hopr] = MTB.read_hr('data/MnBiTe-xu/wannier90_hr_MnBiTe.dat');
%% Set K-path
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
labels={'K','\Gamma','K'}; % labels for k
hkpoints={[0.333333,0.333333,0.000000],...
          [0.0000000000,0.0000000000,0.000000],...
          [-0.333333,-0.333333,0.000000]};% hkpoints-high symmetry k points
nk=101;

%% Calculate bulk bands
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
efermi=0.0258; %% set Fermi Level 0.0258
%%
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"MnBiTe-bulk",0)
%% Calculate bulk bands add electric field
g.readwpos("data/MnBiTe-xu/wpos_MnBiTe");
g.wpos=g.wpos*27
Electric_field_in_evpA=0.0024; %eV/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b)
E_gate=0.0024
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig10",E_gate*10000)
%%
E_gate=0.0024
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig",E_gate*10000)
%% Calculate bulk bands add electric fields
g.readwpos("data/MnBiTe-xu/wpos_MnBiTe")
g.wpos=g.wpos*27
Electric_field_list=linspace(-0.003,0.003,11);
nk=101
for i=1:11
    Electric_field_in_evpA=Electric_field_list(i)
    [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b)
    MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Fig"+num2str(i),Electric_field_in_evpA*10000)
%    Egap(i)=(min(Energy(172,:))-max(Energy(171,:)))*1000
    Egap(i)=(min(Energy(175,:))-max(Energy(174,:)))*1000
    if Egap(i)<1
        Egap(i)=0
    end
end
%%
figure('Color','white')
hold on;
plot(Electric_field_list*10,Egap,'--o','LineWidth',2)
ylabel('E_{gap} (meV)','FontSize',24)
xlabel('Electric field(V/nm)','FontSize',24)
ax=gca;
ax.YAxis.FontSize=20;
ax.XAxis.FontSize=20;
ylim([0.0,70])
box on
print("EvsGap",'-dpng','-r600')

%% DOS cals
nk=101;
Electric_field_in_evpA=0;
Energy=MTB.ham.get_bulk_plane_bands_add_electric(g.ham, g.hopr, g.wpos, Electric_field_in_evpA, nbands,nrpts,nk,g.a,g.b);

%%
nk=101;
Enum=101;
Emin=-0.2;
Emax=0.2;
eps=0.0025;
Nband=276;
plottap=0;
Electric_field_list=linspace(-0.003,0.003,11);
Dos_all=zeros(length(Electric_field_list),Enum);
TDos_all=zeros(length(Electric_field_list),Enum);
for i=1:11
    Electric_field_in_evpA=Electric_field_list(i);
    Energy=MTB.ham.get_bulk_plane_bands_add_electric(g.ham, g.hopr, g.wpos, Electric_field_in_evpA, nbands,nrpts,nk,g.a,g.b);
    [Eaxis,Dos,TDos]=MTB.ham.get_dos(Energy-efermi,eps,Enum,Emin,Emax,Nband,nk,plottap);
    Dos_all(i,:)=Dos;
    TDos_all(i,:)=TDos;
end
%%
figure('Color','white')
%plot(Eaxis,Dos,'k-')
plot(Eaxis,Dos_all(6,:),'Linestyle','-','Color','#4DA1D7','LineWidth',2)
hold on;
%plot(Eaxis,TDos_all(6,:),'Linestyle','-','LineWidth',2)
xlabel('E(eV)')
ylabel('Dos')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
%% Plot Dos
figure('Color','white')
hold on;
TDos_all_new=TDos_all/16.5012*10^16
for i=1:1:11
    plot(Eaxis,TDos_all_new(i,:),'Linestyle','-','LineWidth',2,'Color',[(i-1)/10,100/255 255/255])
end
labels=Electric_field_list
%Eaxis_all=repmat(Eaxis,11,1)
%plot(Eaxis_all',TDos_all','LineStyle','-')
box on;
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
ylabel('n (cm^{-2})','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
legend('-3.0meV','-2.4meV','-1.8meV','-1.2meV','-0.6meV','0meV','0.6meV','1.2meV','1.8meV','2.4meV','3.0meV','Location','best','NumColumns',2)
legend('boxon')
legend('FontSize',13)
xlim([-0.1,0.1])
ylim([0.0,3*10^13])
print('MnBiTe-Density','-dpng','-r300')

%%

%%
figure('Color','white')
hold on;
for i=1:11
    plot(Eaxis,Dos_all(i,:)+i*0.00,'Linestyle','-','LineWidth',2,'Color',[i/11,30/255 255/255])
end
xline(0,'LineStyle','--','LineWidth',2)
legend('-3.0meV','-2.4meV','-1.8meV','-1.2meV','-0.6meV','0meV','0.6meV','1.2meV','1.8meV','2.4meV','3.0meV','Location','northeast','NumColumns',2)
% legend('boxoff')
legend('FontSize',13)
ax=gca;
ax.LineWidth=1;
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
xlim([-0.1,0.1])
%ylim([0.0,0.08])
yticklabels('')
ylabel('DOS','FontSize',24)
xlabel('E-E_f (eV)','FontSize',24)
box on
print('MnBiTe-DOS','-dpng','-r300')