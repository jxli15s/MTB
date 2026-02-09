clc;
clear;

a=1
n=6;
L=n*a;
r=1:1:L;
k=2*pi/L*2;
result=exp(1j*k*r);
result_sum_r=sum(exp(1j*k*r))/length(r);
R=0;
k=2*pi/L*r;
result_k=exp(1j*k*R);
result_sum_k=sum(exp(1j*k*R))/length(k)

figure()
plot(real(result),imag(result),'o')

% plot(real(result_k),imag(result_k),'o')

%%
a=1;
a1=[1/2,sqrt(3)/2]*a;
a2=[-1/2,sqrt(3)/2]*a;
b=2*pi*inv([a1;a2])';
b1=b(1,:);
b2=b(2,:);
n=81;
n1=1:n;
n2=1:n;
[N1,N2]=meshgrid(n1,n2);
nk1=linspace(-1,1,n);
nk2=linspace(-1,1,n);
[Nk1,Nk2]=meshgrid(nk1,nk2);
kx=Nk1*b1(1,1)+Nk2*b2(1,1);
ky=Nk1*b1(1,2)+Nk2*b2(1,2);
rx=N1*a1(1,1)+N2*a2(1,1);
ry=N1*a1(1,2)+N2*a2(1,2);
S=zeros(size(kx));
for i=1:size(kx,1)
    for j=1:size(ky,2)
        k=[kx(i,j),ky(i,j)];
        S(i,j)=get_S_at_single_k(k,rx,ry);
    end
end
%%
k=2*b1+5*b2;
c=get_S_at_single_k(k,rx,ry);
%%
% 绘制相图
figure();
imagesc(nk1,nk2,real(S*S')); % 使用 imagesc 显示网格化数据
% 设置自定义柔和颜色（橙色、灰色、蓝色）
% colormap([1 0.6 0.2; 0.7 0.7 0.7; 0.2 0.4 0.8]); % 橙色、灰色、蓝色
colormap(slanCM('inferno'))
% colormap(flipud(colormap));
shading interp
% caxis([-10, 1000]); % 设置颜色范围（单位: eV）
%%
function sk=get_S_at_single_k(k,rx,ry)
    sk=0;
    for i=1:size(rx,1)
        for j=1:size(ry,2)
            r=[rx(i,j),ry(i,j)];
            sk=sk+exp(1j*k*r');
        end
    end
end