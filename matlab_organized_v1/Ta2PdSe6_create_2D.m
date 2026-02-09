clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Ta2PdSe6");
g = MTB.read_poscar(g,"data/Ta2PdSe6/POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a
%
g.wpos=g.wpos-g.wpos(5,:)
%
theta=-angle(8.8011+5.5846j)
R=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]
R=[1,0,0;0,cos(theta),-sin(theta);0,-sin(theta),sin(theta)]
g.a(2,:)=(R*g.a(2,:)')';
g.a(3,:)=cross(g.a(1,:),g.a(2,:))/norm( cross(g.a(1,:),g.a(2,:)))*20
g.atoms=g.wpos*inv(g.a);
g.atoms(:,3)=g.atoms(:,3)+0.5;

newwpos=R*g.wpos'
newga=[g.a(1,:);g.a(2,:);0,0,20]