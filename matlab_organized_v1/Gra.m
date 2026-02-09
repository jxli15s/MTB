clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/Graphene/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.read_hr('data/Graphene/Graphene_hr.dat');


%% Calculate bulk bands
[nbands,~,nrpts]=size(g.ham);
labels={'M','K','\Gamma','M'}; % labels for k
hkpoints={[0.5,0,0.0],...
          [1/3,1/3,0.0],...
          [0,0,0.0],...
          [0.5,0.0,0.0]};% hkpoints-high symmetry k points
efermi=-1.2533;
nk=101;
g.wpos=[];
g.wpos=g.atoms*g.a
orbital_num=[1,1]
g.wpos=[repmat(g.wpos(1,:),[orbital_num(1),1]);...
        repmat(g.wpos(2,:),[orbital_num(2),1]);...
    ]
Electric_field_in_evpA=0.00;
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"Gra-2d",Electric_field_in_evpA*10000);
%%
nk=81;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
Electric_field_in_evpA=0;
nsband=1:1;
% [Unk,Enk]=MTB.ham.get_bulk_plane_bands_add_electric(g,Electric_field_in_evpA,Kx,Ky,Kz);
% [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
%%
[Omega_k,KX,KY,KZ]=MTB.ham.get_Berrycurvature_cop(g,Kx,Ky,Kz,Enk,Unk,nsband,0.01,1);

%% Calculate Berry Curvature by LOOP method
nk=100;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,nk);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
plottap=1;
bandindex=1;
[Omega_k,KX,KY] = MTB.ham.get_Berry_curvature(bandindex,Unk,Kx,Ky,plottap);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Check Time Reversal Symmetry            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T^{-1}conj(H(k))T=H(Tk)=H(-k) T=i*sigma_y*k Here no spin only orbital so
% T=UK=IK
% TH(k)T=H(-k) UH(k)*U^\dagger=H(-k)
clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/Graphene/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/wannier90_hr_p1.dat','data/Graphene/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;
[nbands,~,nrpts]=size(g.ham);

kpoint=[0.3,0.1,0.0];
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=-kpoint;
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
has_spin=false;
T = get_time_reversal_operator(nbands, has_spin);
%%
H_Trev = T.apply(hk1);
h=H_Trev-hk2;
max(h,[],'all')
% h1=T*conj(hk1)*inv(T)-hk2; %% TH=HT T=i*sigma_y*K  THT^-1=H 
% max(h1,[],'all')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check Inversion Symmetry for Gra     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PH(k)P=H(-k) Pz parity=-1
[nbands,~,nrpts]=size(g.ham);
s1=[0  -1
    -1  0];
kpoint=[0.5,0.2,0.0];
[~,~,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
kpoint=[-0.5,-0.2,0.0];
[~,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
P=s1;
c=P*hk1*inv(P)-hk2; %% PH(k)P^{-1}=H(Pk)
max(c,[],'all')
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Check Rotation Symmetry for Gra      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CH(k)C^{}
theta_deg=120;
sigma_z=[1,0;0,1];
theta = theta_deg * pi / 180;
R = rotation_matrix(theta);
kpoint1=[0.333333,0.333333,0.0];
C = construct_rotation_operator_periodic(g.wpos, g.a, g.b,R, kpoint1);
% kpoint=[0.5,0.2,0.0];
% [~,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);%[Energy,Psik,hk]
[~,Psik,hk1]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint1,g.a,g.b);
kpoint2=(R*(kpoint1*g.b)')'*inv(g.b);
% [~,~,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
[~,~,hk2]=MTB.ham.get_parity_singleK_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,kpoint2,g.a,g.b);
C=expm(1i*2*pi/3*sigma_z)
h=C*hk1*C'-hk2
%%
cn=eig(Psik(:,1:2)'*C*Psik(:,1:2))
Jz=angle(cn)/theta

%%
function T = get_time_reversal_operator(Norb, has_spin)
% Norb: 总轨道数（含自旋）
% has_spin: 是否包含自旋
if ~has_spin
    T.U = eye(Norb);
    T.K = @(x) conj(x);
    T.apply = @(Hk) conj(Hk);  % 只做复共轭
else
    % 假设轨道排序为: [up1, dn1, ..., upx, dnx, ...]
    N = Norb / 2;
    sigma_y = [0, -1i; 1i, 0];
    U = kron(eye(N), -1i*sigma_y);  % 2N x 2N
    T.U = U;
    T.K = @(x) conj(x);
    T.apply = @(Hk) U' * conj(Hk) * U;
end
end

function P = get_inversion_operator(gs)
% 构造轨道空间的反演算符 P

Norb = size(gs.wpos, 1);
P = zeros(Norb);

% 原子坐标（假设在笛卡尔空间）
wpos = gs.wpos(:, 1:3);

% 对每个轨道 i 找到其反演后的坐标 -r_i
for i = 1:Norb
    ri = wpos(i, :);
    r_inv = -ri;

    % 找与 -ri 最接近的轨道 j
    dists = vecnorm(wpos - r_inv, 2, 2);
    [~, j] = min(dists);

    % 设为 +1（偶宇称）或 -1（奇宇称），可以根据轨道类型定义
    parity = -1;  % 对于 p_z
    P(j, i) = parity;
end
end


function C = construct_rotation_operator_periodic(wpos, a, b, R, k)
% wpos: Norb×3，轨道位置（笛卡尔坐标）
% a: 3×3 实空间晶格基矢
% R: 3×3 旋转矩阵
% k: 1×3 动量点
Norb = size(wpos, 1);
C = zeros(Norb);

% 平移搜索范围：-1, 0, 1
range = -2:2;
[X, Y, Z] = ndgrid(range, range, 0);  % 平移方向：xy方向±1
shifts = [X(:), Y(:), Z(:)];  % 共 3×3×1 = 9 个平移
Tlist = shifts * a;  % 将分数坐标转为笛卡尔坐标

wpos_rot = (R * wpos')';  % 所有轨道旋转后的位置
k=k*b;
for i = 1:Norb
    ri_rot = wpos_rot(i, 1:3);  % 旋转后的轨道位置

    min_dist = inf;
    j_best = -1;
    T_best = [0 0 0];

    for j = 1:Norb
        for t = 1:size(Tlist,1)
            T = Tlist(t,:);
            diff = wpos(j,:) - (ri_rot + T);
            d = norm(diff);
            if d < min_dist
                min_dist = d;
                j_best = j;
                T_best = T;
            end
        end
    end

    if min_dist > 1e-3
        error('No matching site found for orbital %d under rotation.', i);
    end


    % R=round((ri_rot-wpos(j_best,:))*inv(a));
    % phase = exp(1j * dot(k, R));
    % C(j_best, i) = phase;
    phase = exp(1j * dot(k, wpos(j_best,:) - ri_rot));
    C(j_best, i) = phase;

end
end



function R = rotation_matrix(theta)
% 绕 z 轴旋转 theta 弧度
R = eye(3);
R(1:2,1:2) = [cos(theta), -sin(theta); sin(theta), cos(theta)];
end
