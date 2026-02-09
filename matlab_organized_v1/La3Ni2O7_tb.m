clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get band structure of TB La3Ni2O7               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

efermi=0
hkpoints=[0.0,0.0;0.5,0.0;0.5,0.5;0.0,0.0]*2*pi;
knum=20;
kpoints=[];
for i=1:size(hkpoints,1)-1
    k1=hkpoints(i,:);
    k2=hkpoints(i+1,:);
    kpoints=[kpoints;linspace(k1(1),k2(1),knum)',linspace(k1(2),k2(2),knum)'];
end
Energy=zeros(4,size(kpoints,1))
parfor i=1:length(kpoints)
    kx=kpoints(i,1)
    ky=kpoints(i,2)
    hk=get_hk(kx,ky)
    vals=sort(eig(hk));
    Energy(:,i)=vals;
end

kpath=linspace(0,1,size(kpoints,1))
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'Color','black','LineWidth',2);
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Get Parity of La3Ni2O7               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tx2 tz2 bx2 bz2
hk=get_hk(0.0*2*pi,0.00*2*pi);
[V,D]=eig(hk);
[Energy,ind]=sort(diag(D));
Psik=V(:,ind);
s1=[0,1;1,0]
P=kron(s1,eye(2))
b=Psik'*P*Psik

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get band structure of TB La3Ni2O7               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hkpoints=[0.0,0.0;0.5,0.0;0.5,0.5;0.0,0.0]*2*pi;
knum=20;
kpoints=[];
for i=1:size(hkpoints,1)-1
    k1=hkpoints(i,:);
    k2=hkpoints(i+1,:);
    kpoints=[kpoints;linspace(k1(1),k2(1),knum)',linspace(k1(2),k2(2),knum)'];
end
Energy=zeros(3,size(kpoints,1))
U=[1,0,0,0;0,1/sqrt(2),0,1/sqrt(2);0,0,1,0;0,1/sqrt(2),0,-1/sqrt(2)];
parfor i=1:length(kpoints)
    kx=kpoints(i,1)
    ky=kpoints(i,2)
    hk=get_hk(kx,ky)
    hk=U*hk*inv(U)
    hk=hk(1:3,1:3);
    vals=sort(eig(hk));
    Energy(:,i)=vals;
end

kpath=linspace(0,1,size(kpoints,1))
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'Color','black','LineWidth',2);
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Get Parity of La3Ni2O7               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tx2 tz2 bx2 bz2
hk=get_hk(0.0*2*pi,0.00*2*pi);
hk=U*hk*inv(U);
hk=hk(1:3,1:3);
[V,D]=eig(hk);
[Energy,ind]=sort(diag(D));
Psik=V(:,ind);
P=[0,0,1;0,1,0;1,0,0];
b=Psik'*P*Psik



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get band structure of TB La3Ni2O7               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hkpoints=[0.0,0.0,0.0;0.5,0.0,0.0;0.5,0.5,0.0;0.0,0.0,0.0;0.0,0.0,0.5]*2*pi;
knum=20;
kpoints=[];
for i=1:size(hkpoints,1)-1
    k1=hkpoints(i,:);
    k2=hkpoints(i+1,:);
    kpoints=[kpoints;linspace(k1(1),k2(1),knum)',linspace(k1(2),k2(2),knum)',linspace(k1(3),k2(3),knum)'];
end
Energy=zeros(8,size(kpoints,1));
hk=[]
parfor i=1:length(kpoints)
    kx=kpoints(i,1)
    ky=kpoints(i,2)
    kz=kpoints(i,3)
    hk=get_hk_3d(kx,ky,kz);
    vals=sort(eig(hk));
    Energy(:,i)=vals;
end

kpath=linspace(0,1,size(kpoints,1))
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'Color','black','LineWidth',2);
    hold on
end

%%


function hk=get_hk(kx,ky)
ex2=10.929;
t11x=-0.505;
t11xy=0.060;
t11xx=-0.047;
ez2=10.474;
t22x=-0.126;
t22xy=-0.003;
t22xx=-0.021;
s110=-0.049;
s11x=0.001;
s11xy=0.035;
s11xx=-0.022;
s220=-0.628;
s22x=0.011;
s22xy=-0.032;
s22xx=0.024;
t12x=0.253;
t12xx=0.037;
s12x=-0.038;
s12xx=-0.009;

h11=ex2+2*t11x*(cos(kx)+cos(ky))+4*t11xy*cos(kx)*cos(ky)+2*t11xx*(cos(2*kx)+cos(2*ky));
h33=h11;
h22=ez2+2*t22x*(cos(kx)+cos(ky))+4*t22xy*cos(kx)*cos(ky)+2*t22xx*(cos(2*kx)+cos(2*ky));
h44=h22;
h12=2*t12x*(cos(kx)-cos(ky))+2*t12xx*(cos(2*kx)-cos(2*ky));
h13=s110+2*s11x*(cos(kx)+cos(ky))+4*s11xy*cos(kx)*cos(ky)+2*s11xx*(cos(2*kx)+cos(2*ky));
h14=2*s12x*(cos(kx)-cos(ky))+2*s12xx*(cos(2*kx)-cos(2*ky));
h24=s220+2*s22x*(cos(kx)+cos(ky))+4*s22xy*cos(kx)*cos(ky)+2*s22xx*(cos(2*kx)+cos(2*ky));
h23=h14;
h34=h12;
hk=[h11,h12,h13,h14;h12',h22,h23,h24;h13',h23',h33,h34;h14',h24',h34',h44];
end


function hk_3d=get_hk_3d(kx,ky,kz)



ex2=10.920;
t11x=-0.512;
t11xy=0.062;
t11xx=-0.053;
ez2=10.501;
t22x=-0.123;
t22xy=-0.005;
t22xx=-0.022;
s110=-0.030;
s11x=0.001;
s11xy=0.023;
s11xx=-0.015;
s220=-0.601;
s22x=0.025;
s22xy=-0.011;
s22xx=0.011;
t12x=0.250;
t12xx=0.035;
s12x=-0.038;
s12xx=-0.006;
s11z=0.001;
s13z=-0.0004;
s22z=-0.006;
s24z=-0.031;

h11=ex2+2*t11x*(cos(kx)+cos(ky))+4*t11xy*cos(kx)*cos(ky)+2*t11xx*(cos(2*kx)+cos(2*ky));
h33=h11;
h22=ez2+2*t22x*(cos(kx)+cos(ky))+4*t22xy*cos(kx)*cos(ky)+2*t22xx*(cos(2*kx)+cos(2*ky));
h44=h22;
h12=2*t12x*(cos(kx)-cos(ky))+2*t12xx*(cos(2*kx)-cos(2*ky));
h13=s110+2*s11x*(cos(kx)+cos(ky))+4*s11xy*cos(kx)*cos(ky)+2*s11xx*(cos(2*kx)+cos(2*ky));
h14=2*s12x*(cos(kx)-cos(ky))+2*s12xx*(cos(2*kx)-cos(2*ky));
h24=s220+2*s22x*(cos(kx)+cos(ky))+4*s22xy*cos(kx)*cos(ky)+2*s22xx*(cos(2*kx)+cos(2*ky));
h34=h12;

h13=h13*exp(1j*0.3051*kz);
h14=h14*exp(1j*0.3051*kz);
h24=h24*exp(1j*0.3051*kz);
h23=h14;
hk_2d=[h11,h12,h13,h14;h12',h22,h23,h24;h13',h23',h33,h34;h14',h24',h34',h44];
hz11=8*s11z*cos(kz/2)*cos(kx/2)*cos(ky/2);
hz22=8*s22z*cos(kz/2)*cos(kx/2)*cos(ky/2);
hz13=4*s13z*exp(-1j*(1/2-0.3051)*kz)*cos(kx/2)*cos(ky/2);
hz24=4*s24z*exp(-1j*(1/2-0.3051)*kz)*cos(kx/2)*cos(ky/2);
hz=[hz11,0,hz13,0;...
    0,hz22,0,hz24;...
    hz13',0,hz11,0;...
    0,hz24',0,hz22];
hk_3d=[hk_2d,hz;hz',hk_2d];
end

