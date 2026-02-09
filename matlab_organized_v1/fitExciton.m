clc;
clear;
% a=[[1/2,sqrt(3)/2];[-1/2,sqrt(3)/2]]
a=[[1,0];[0.5,sqrt(3)/2]];
b=inv(a')*2*pi;
%%
figure()
hold on;
quiver(0, 0, a(1,1), a(1,2), 'AutoScale', 'on',Color='red')
quiver(0, 0, a(2,1), a(2,2), 'AutoScale', 'on',Color='red')

quiver(0, 0, b(1,1), b(1,2), 'AutoScale', 'on',Color='blue')
quiver(0, 0, b(2,1), b(2,2), 'AutoScale', 'on',Color='blue')
xlim([-8,8]);
ylim([-8,8]);
%%
x=0:0.01:1;
y=0:0.01:1;
[x,y]=meshgrid(x,y);
X=x*a(1,1)+y*a(2,1);
Y=x*a(1,2)+y*a(2,2);
%%

g1=b(2,:);
g2=-b(1,:);
g3=g2-g1;
% g1=b(1,:)+b(2,:);
% g2=b(1,:);
% g3=g2-g1;
g4=-g1;
g5=-g2;
g6=-g3;
d1=X*g1(1)+Y*g1(2);
d2=X*g2(1)+Y*g2(2);
d3=X*g3(1)+Y*g3(2);
d4=X*g4(1)+Y*g4(2);
d5=X*g5(1)+Y*g5(2);
d6=X*g6(1)+Y*g6(2);

% V1=2.3*exp(1j*deg2rad(30.8))
% V1=1.116*exp(1j*deg2rad(29.5)) %MoSe2 AA norelax 
% V1=1.498*exp(1j*deg2rad(64.72)) %MoSe2 AA relax
% V1=-0.7606*exp(1j*deg2rad(-57.03))% MoSe2 AB relax
% V1=0.3431*exp(1j*deg2rad(30.42)) %MoS2 AA norelax
% V1=1.05*exp(1j*deg2rad(67.11)) %MoS2 AA relax
% V1=0.18*exp(1j*deg2rad(0)) %MoS2 AB norelax
 V1=0.8832*exp(1j*deg2rad(98.51)) %MoS2 AB relax
V2=conj(V1)
% V1=1
% V2=1
Z=V1*exp(1j*d1)+V2*exp(1j*d4)+ ...
  V2*exp(1j*d2)+V1*exp(1j*d5)+ ...
  V1*exp(1j*d3)+V2*exp(1j*d6);

figure('Color','white')
% contour(X, Y, Z, 50)
surf(X, Y, Z)
%pcolor(X,Y,Z)
xlim([0,1.6]);
ylim([0,1.6]);
shading interp
colormap(slanCM('RdBu'))
colorbar;
axis off;
grid off;
view(0,90)
axis off;
% print('exciton_MoSe2','-dpng','-r600')
% print('exciton_MoSe2_relax','-dpng','-r600')
print('exciton_MoS2_relax_AB','-dpng','-r600')
% caxis([-10,10])
%% For moire
theta=deg2rad(10);
% theta=1;
am =[[1,0];[0.5,sqrt(3)/2]]/theta;
tmp=cross([g1,0],[0,0,1])
gm1=tmp(1:2)*theta
tmp=cross([g2,0],[0,0,1])
gm2=tmp(1:2)*theta
tmp=cross([g3,0],[0,0,1])
gm3=tmp(1:2)*theta
tmp=cross([g4,0],[0,0,1])
gm4=tmp(1:2)*theta
tmp=cross([g5,0],[0,0,1])
gm5=tmp(1:2)*theta
tmp=cross([g6,0],[0,0,1])
gm6=tmp(1:2)*theta
%% Fit
x=0:0.01:1;
y=0;
X=x*am(1,1)+y*am(2,1);
Y=x*am(1,2)+y*am(2,2);
d1=X*gm1(1)+Y*gm1(2);
d2=X*gm2(1)+Y*gm2(2);
d3=X*gm3(1)+Y*g3(2);
d4=X*gm4(1)+Y*gm4(2);
d5=X*gm5(1)+Y*gm5(2);
d6=X*gm6(1)+Y*gm6(2);
Z=V1*exp(1j*d1)+V2*exp(1j*d4)+ ...
  V2*exp(1j*d2)+V1*exp(1j*d5)+ ...
  V1*exp(1j*d3)+V2*exp(1j*d6);

figure('Color','white')
plot(X,Z)
%%
x=0:0.01:5;
y=0:0.01:5;
[x,y]=meshgrid(x,y);
X=x*am(1,1)+y*am(2,1);
Y=x*am(1,2)+y*am(2,2);

d1=X*gm1(1)+Y*gm1(2);
d2=X*gm2(1)+Y*gm2(2);
d3=X*gm3(1)+Y*gm3(2);
d4=X*gm4(1)+Y*gm4(2);
d5=X*gm5(1)+Y*gm5(2);
d6=X*gm6(1)+Y*gm6(2);
Z=V1*exp(1j*d1)+V2*exp(1j*d4)+ ...
  V2*exp(1j*d2)+V1*exp(1j*d5)+ ...
  V1*exp(1j*d3)+V2*exp(1j*d6);

figure('Color','white')
% contour(X, Y, Z, 50)
surf(X, Y, Z)
%pcolor(X,Y,Z)
% xlim([13.2,26.5]);
% ylim([4.2,18.5]);
xlim([13.5,29]);
ylim([6.5,22]);
% % xlim([16,30]);
% % ylim([4,19]);
% xlim([0,1.6]);
% ylim([0,1.6]);
shading interp
colormap(slanCM('RdBu'))
colorbar;
shading interp;
view(0,90)
axis off;
% print('exciton_MoSe2_moire_relax','-dpng','-r600')
% print('exciton_MoSe2_moire','-dpng','-r600')
print('exciton_MoS2_moire_relax_AB','-dpng','-r600')
%% Fit
a=[[1,0];[0.5,sqrt(3)/2]];
b=inv(a')*2*pi;
x=[0,0.333333]
y=[0,0.333333]
X=x*a(1,1)+y*a(2,1);
Y=x*a(1,2)+y*a(2,2);
Z=[10,-10]
g1=b(2,:);
g2=-b(1,:);
g3=g2-g1;
g4=-g1;
g5=-g2;
g6=-g3;

d1=X*g1(1)+Y*g1(2);
d2=X*g2(1)+Y*g2(2);
d3=X*g3(1)+Y*g3(2);
d4=X*g4(1)+Y*g4(2);
d5=X*g5(1)+Y*g5(2);
d6=X*g6(1)+Y*g6(2);
V1=1.116*exp(1j*deg2rad(30))
V2=conj(V1)

myfittype=fittype(@(V,phi,g11,g12,g21,g22,g31,g32,g41,g42,g51,g52,g61,g62,X,Y)  V*exp(1j*deg2rad(phi)*(X*g11+Y*g12))+ ...
                                V*exp(1j*deg2rad(phi)*(X*g31+Y*g32))+ ...
                                V*exp(1j*deg2rad(phi)*(X*g51+Y*g52))+...
                               conj(V)*exp(1j*deg2rad(phi)*(X*g21+Y*g22))+...
                               conj(V)*exp(1j*deg2rad(phi)*(X*g41+Y*g42))+...
                               conj(V)*exp(1j*deg2rad(phi)*(X*g61+Y*g62)), ...
                               'problem', {'g11','g12','g21','g22','g31','g32', ...
                               'g41','g42','g51','g52','g61','g62'},...
                               'independent', {'X','Y'}, ...
                               'dependent', 'Z')

f2 = fit( [X',Y'] ,Z, myfittype, 'problem', {g1(1),g1(2),g2(1),g2(2),g3(1),g3(2),g4(1),g4(2),g5(1),g5(2),g6(1),g6(2)}, 'start', [1, 2] )



