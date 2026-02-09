clear;
clear all;
load("data/SrTiIrO/inplane/inplane.mat")
%%
efermi=6.5146;
energy=energy-6.5146;

%%
ef=-0.25;
emin=ef-0.003;
emax=ef+0.003;
energy_ef=energy(energy>emin & energy<emax);
kx_ef=kx(energy>emin & energy<emax);
ky_ef=ky(energy>emin & energy<emax);
kz_ef=kz(energy>emin & energy<emax);
sx_ef=sx(energy>emin & energy<emax);
sy_ef=sy(energy>emin & energy<emax);
sz_ef=sz(energy>emin & energy<emax);

kx_mtx=[kx_ef;kx_ef+1.0;kx_ef;kx_ef+1.0];
ky_mtx=[ky_ef;ky_ef;ky_ef+1.0;ky_ef+1.0];
kz_mtx=repmat(kz_ef,4,1);
sx_mtx=repmat(sx_ef,4,1);
sy_mtx=repmat(sy_ef,4,1);
sz_mtx=repmat(sz_ef,4,1);
%%
 %quiver3(kx_ef,ky_ef,kz_ef,sx_ef./abs(sx_ef),sy_ef./abs(sy_ef),sz_ef./abs(sz_ef))
% quiver3(kx_mtx,ky_mtx,kz_mtx,sx_mtx,sy_mtx,sz_mtx)
quiver(kx_mtx,ky_mtx,sx_mtx,sy_mtx,0.25,'r');
q.ShowArrowHead = 'off';
ylim([0.3,0.7])
xlim([0.3,0.7])


%% Load from source readmatrix

data=readmatrix("data/SrTiIrO/inplane/band-structure-all-101-101.dat");
%data=readmatrix("data/SrTiIrO/inplane/band-structure-all.dat");
%%
efermi=6.5146;
data2=data(:,[3,4,5,6,7,41,58,75]);
data2(:,5)=data2(:,5)-efermi;
data2=data2(data2(:,4)==520,:);
% kx ky kz band_index energy sx sy sz
ef=-0.25;
emin=ef-0.0025;
emax=ef+0.0025;
data_ef=data2(data2(:,5)>emin & data2(:,5)<emax,:);
kx_mtx=[data_ef(:,1);data_ef(:,1)+1.0;data_ef(:,1);data_ef(:,1)+1.0];
ky_mtx=[data_ef(:,2);data_ef(:,2);data_ef(:,2)+1.0;data_ef(:,2)+1.0];
kz_mtx=repmat(data_ef(:,3),4,1);
sx_mtx=repmat(data_ef(:,6),4,1);
sy_mtx=repmat(data_ef(:,7),4,1);
sz_mtx=repmat(data_ef(:,8),4,1);
ef_mtx=repmat(data_ef(:,5),4,1);
% kx_mtx=kx_mtx(kx_mtx>0.5);
% kx_mtx=[kx_mtx;-kx_mtx];
% ky_mtx=kx_mtx(kx_mtx>0.5);
quiver(kx_mtx,ky_mtx,sx_mtx,sy_mtx,0.5,'r');
ylim([-0.5,0.5])
xlim([-0.5,0.5])
% print("SrTiIrO_spintexture.png",'-dpng','-r600')
hold on
%%
print("SrTiIrO_spintexture.eps",'-depsc','-vector')
% print(figure_current,'-depsc','-painters',file_figure_eps_full_location);
%% plot FermiSurface
efermi=6.5146;
data2=data(:,[3,4,5,6,7,41,58,75]);
data2(:,5)=data2(:,5)-efermi;
ef=-0.25;
emin=ef-0.0025;
emax=ef+0.0025;
data_ef=data2(data2(:,5)>emin & data2(:,5)<emax,:);
kx_mtx2=[data_ef(:,1);data_ef(:,1)+1.0;data_ef(:,1);data_ef(:,1)+1.0];
ky_mtx2=[data_ef(:,2);data_ef(:,2);data_ef(:,2)+1.0;data_ef(:,2)+1.0];
ef_mtx2=repmat(data_ef(:,5),4,1);
ef_band_index=repmat(data_ef(:,4),4,1);
data_new=[kx_mtx2,ky_mtx2,ef_band_index,ef_mtx2];
% figure('Color','White')
data_new=data_new(kx_mtx2>0.3 & kx_mtx2<0.7 & ky_mtx2>0.3 & ky_mtx2<0.7 & ef_band_index==517,:,:,:,:)
data_new=sortrows(data_new,1)
data_new_1=data_new(data_new(:,1)>0.5,:)
data_new_1=sortrows(data_new_1,2)
data_new_2=data_new(data_new(:,1)<0.5,:)
data_new_2=sortrows(data_new_2,-2)
data_new=[data_new_1;data_new_2;data_new_1(1,:)]
% data_new=sortrows(data_new,2)
plot(data_new(:,1),data_new(:,2),'*-')
hold on
%plot(kx_mtx2,ky_mtx2,'o','MarkerSize',1,'MarkerFaceColor','red','MarkerEdgeColor','red')

% quiver(data_ef(:,1),data_ef(:,2),data_ef(:,6),data_ef(:,7),0.25,'r');

%%

kx=reshape(data2(:,1),101,101,624);
ky=reshape(data2(:,2),101,101,624);
energy_fs=reshape(data2(:,5),101,101,624);
energy_fs=kron(ones(2,2),energy_fs(:,:,519))
contour(energy_fs(:,:))
% contour(kx(:,:,1),ky(:,:,1),energy_fs(:,:,520))
