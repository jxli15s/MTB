%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Construct BHZ kp Hamiltonian                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%=== Pauli matrices and gamma matrices ===
s0 = eye(2);
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];

% Gamma matrices for 4×4 BHZ model
Gamma1 = kron(sx, sz);  % σ_x ⊗ s_z
Gamma2 = kron(sy, s0);  % σ_y ⊗ I
Gamma5 = kron(sz, s0);  % σ_z ⊗ I

%=== Define parameters and symbolic variables ===
syms A B M real
syms kx ky kz real

%=== Define H(k) with lattice regularization ===
% H(k) = A*(sin(kx)Γ1 + sin(ky)Γ2) + [M + 2B(2 - cos(kx) - cos(ky))] * Γ5
Hk = A * sin(kx) * Gamma1 ...
   + A * sin(ky) * Gamma2 ...
   + (M + 2*B*(2 - cos(kx) - cos(ky))) * Gamma5;

%=== Define hopping directions (max: NN) ===
deltas = [
     0,  0;
     1,  0;
    -1,  0;
     0,  1;
     0, -1;
];
fprintf("Computing tight-binding hopping terms from H(k)...\n");

% === Compute hopping terms by inverse Fourier transform ===
hoppings = {};
for i = 1:size(deltas,1)
    dx = deltas(i,1);
    dy = deltas(i,2);

    % Phase factor: e^{-i(k·delta)}
    phase = exp(-1i * (kx*dx + ky*dy));

    % Inverse FT: T(delta) = ∫∫ H(k) e^{-i k·δ} dkx dky / (2π)^2
    Tdelta = int(int(Hk * phase, kx, -pi, pi), ky, -pi, pi) / (2*pi)^2;

    % Simplify and convert to numeric coefficients but keep A, B, M symbolic
    Tdelta = vpa(simplify(Tdelta, 'Steps', 50), 6);

    % Store
    hoppings{end+1,1} = [dx, dy];
    hoppings{end,2} = Tdelta;
end

% === Display results ===
fprintf("\nTight-binding hopping terms (numerical coefficients, symbolic parameters):\n");
for i = 1:length(hoppings)
    delta = hoppings{i,1};
    T = hoppings{i,2};
    fprintf("δ = [%2d, %2d]:\n", delta(1), delta(2));
    disp(T);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Construct BHZ kp and Plot bands                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
% === Create BHZ 4 band kp bands ===
syms A B M real
syms kx ky kz real
H_sym = bhz4band_tb_symbolic(kx, ky, kz, M, A, B);
% === Define values of variables ===
m=0.1;a=0.2;b=-0.2;
params = {
    M, m;
    A, a;
    B, b;
};
% === Transform the symbolic to functions ===
H_fixed = subs(H_sym, params(:,1), params(:,2));
H_func = matlabFunction(H_fixed,'Vars',{kx,ky,kz});

% === Define high-symmetry path Γ–M–K ===
G = [0,  0,   0];
M = [pi, 0, 0];
K = [1/3*pi,1/3*pi,0];
kpts = [M; G; K];
labels = {'M', '\Gamma', 'K'};

% === Calculate the eigenvalues along the HSP ===
Nk = 100;  % points per segment
[Energy, kdist, k_tick] = get_bands_kp(H_func, kpts, labels, Nk);
efermi=0.0;
nbands=size(H_func(0,0,0),1);
MTB.plot.plot_bands_Electric(Energy',nbands,efermi,kdist,labels,k_tick,"kp",0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              Construct FeSe 6 bands and Plot bands                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% === Pauli matrices and gamma matrices ===
s0 = eye(2);
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];
Gamma1 = kron(sx, sz);  % sigma_x ⊗ s_z
Gamma2 = kron(sy, s0);  % sigma_y ⊗ I
Gamma3 = kron(sx, sx);  % sigma_x ⊗ s_x
Gamma4 = kron(sx, sy);  % sigma_x ⊗ s_y
Gamma5 = kron(sz, s0);  % sigma_z ⊗ I

% === Define parameters and symbolic variables ===
syms M10 M11 M12 M20 M21 M22 A1 A2 B1 B2 dso real
syms kx ky kz real
M1=M10+M11*(kx^2+ky^2)+M12*kz^2
M2=M20+M21*(kx^2+ky^2)+M22*kz^2
D=B1*(kx^2-ky^2)-1i*B2*kx*ky;
kp=kx+1i*ky
km=kx-1i*ky
% === Define H_sym with lattice regularization ===
H_sym = [    M1,       0,  A2*kz,   -A1*km,   A1*kp,        0;...
             0,      M1,  A1*kp,    A2*kz,       0,   -A1*km;...
         A2*kz,   A1*km,     M2,       0,        0,  conj(D);...
        -A1*kp,   A2*kz,      0,       M2,       D,        0;...
         A1*km,       0,      0,  conj(D),  M2+dso,        0;...
             0,  -A1*kp,      D,        0,       0,   M2+dso;
    ];

% === Define values of symbolic variables ===
m10=1;m11=m10;m12=-0.5*m10;
m20=-m10;m21=-m11;m22=-m12;
a1=0.5;a2=0.1;b1=0;b2=0;dsoc=0.5;
k1=0.1;k2=0.2;k3=0;

% Substitute parameters into H_sym % H_func = matlabFunction(H_sym, 'Vars', {kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso});
H_fixed = subs(H_sym, {M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso}, {m10, m11, m12, m20, m21, m22, a1, a2, b1, b2, dsoc});
H_func = matlabFunction(H_fixed,'Vars',{kx,ky,kz});

% === Define high-symmetry path M-Γ–Z ===
G = [0,  0,   0];
Z = [0, 0,  pi];
M = [pi, pi, 0];
kpts = [M; G; Z];
labels = {'M', '\Gamma', 'Z'};

% Interpolation
Nk = 100;  % points per segment
klist = [];
k_tick = 0;
for i = 1:size(kpts,1)-1
    seg_kx = linspace(kpts(i,1), kpts(i+1,1), Nk)';
    seg_ky = linspace(kpts(i,2), kpts(i+1,2), Nk)';
    seg_kz = linspace(kpts(i,3), kpts(i+1,3), Nk)';
    klist = [klist; [seg_kx, seg_ky, seg_kz]];
end

% Compute path length (k distance)
kdist = [0; cumsum(sqrt(sum(diff(klist).^2, 2)))];

% Compute x-ticks (positions of high-symmetry points)
k_tick = [0];
for i = 2:size(kpts,1)
    d = norm(kpts(i,:) - kpts(i-1,:));
    k_tick(end+1) = k_tick(end) + d;
end

% === Calculate band energies ===
nBands = size(H_func(0,0,0),1);
Ebands = zeros(length(klist), nBands);
for i = 1:length(klist)
    kx_val = klist(i,1);
    ky_val = klist(i,2);
    kz_val = klist(i,3);
    H = H_func(kx_val, ky_val, kz_val);
    Ebands(i,:) = sort(real(eig(H)))';
end

% === Plot band structure ===
figure;
plot(kdist, Ebands, 'LineWidth', 1.5);
xlabel('k-path'); ylabel('Energy');
title('Band structure of 4×4 lattice model');
xticks(k_tick);
xticklabels(labels);
xlim([kdist(1), kdist(end)]);
grid on;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Construct FeSe 6 bands tb and get real space hoppings in 3D     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% === Define parameters and symbolic variables ===
syms M10 M11 M12 M20 M21 M22 A1 A2 B1 B2 dso real
syms kx ky kz real
M1=M10+M11*(2*(1-cos(kx))+2*(1-cos(ky)))+M12*2*(1-cos(kz));
M2=M20+M21*(2*(1-cos(kx))+2*(1-cos(ky)))+M22*2*(1-cos(kz));
D=B1*(2*(1-cos(kx))-2*(1-cos(ky)))-1i*B2*sin(kx)*sin(ky);
kp=sin(kx)+1i*sin(ky);
km=sin(kx)-1i*sin(ky);
% === Define H_sym with lattice regularization ===
H_sym = [    M1,       0,  A2*sin(kz),   -A1*km,   A1*kp,        0;...
             0,      M1,  A1*kp,    A2*sin(kz),       0,   -A1*km;...
         A2*sin(kz),   A1*km,     M2,       0,        0,  conj(D);...
        -A1*kp,   A2*sin(kz),      0,       M2,       D,        0;...
         A1*km,       0,      0,  conj(D),  M2+dso,        0;...
             0,  -A1*kp,      D,        0,       0,   M2+dso;
    ];
% === Define hopping directions (max: NN) ===
deltas = [
    0 0 0;   % on-site
    1 0 0;  
   -1 0 0;
    0 1 0;   
   0 -1 0;
    0 0 1;  
   0 0 -1
];

fprintf("Computing tight-binding hopping terms from H_sym(kx,ky,kz)...\n");

hoppings = {};
for i = 1:size(deltas,1)
    dx = deltas(i,1);
    dy = deltas(i,2);
    dz = deltas(i,3);
    % Phase factor: e^{-i(k·delta)}
    phase = exp(-1i * (kx*dx + ky*dy + kz*dz));

    % Inverse Fourier transform
    Tdelta = int(int(int(H_sym * phase, kx, -pi, pi), ky, -pi, pi), kz, -pi, pi) / (2*pi)^3;

    % Simplify & store
    Tdelta = vpa(simplify(Tdelta, 'Steps', 50), 6);
    hoppings{end+1,1} = [dx, dy, dz];
    hoppings{end,2} = Tdelta;
end

% === Display results ===
fprintf("\nTight-binding hopping terms (numerical coefficients, symbolic parameters):\n");
for i = 1:length(hoppings)
    delta = hoppings{i,1};
    T = hoppings{i,2};
    fprintf("δ = [%2d, %2d, %2d]:\n", delta(1), delta(2),delta(3));
    disp(T);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       Construct FeSe 6 bands tb and get bands by functions        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
syms M10 M11 M12 M20 M21 M22 A1 A2 B1 B2 dso real
syms kx ky kz real
% get the symbolic functions
H_sym = FeSe_6bandkp_symbolic(kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso);
%  Define values
m10=1;m11=m10;m12=-0.5*m10;
m20=-m10;m21=-m11;m22=-m12;
a1=0.5;a2=0.1;b1=0;b2=0;dsoc=0.5;
k1=0.1;k2=0.2;k3=0;
% Substitute the symbolic parameters
params = {
    M10, m10;
    M11, m11;
    M12, m12;
    M20, m20;
    M21, m21;
    M22, m22;
    A1, a1;
    A2, a2;
    B1, b1;
    B2, b2;
    dso, dsoc;
};
H_fixed = subs(H_sym, params(:,1), params(:,2));
H_func = matlabFunction(H_fixed,'Vars',{kx,ky,kz});

nbands=size(H_func(0,0,0),1);
efermi=0.0;
% === Define high-symmetry path Γ–X–M–Γ ===
G = [0,  0,   0];
Z = [0, 0,  pi];
M = [pi, pi, 0];
kpts = [M; G; Z];
labels = {'M', '\Gamma', 'Z'};
Nk = 100;  % points per segment

% Calculate the band eigenvalues
[Energy, kdist, k_tick] = get_bands_kp(H_func, kpts, labels, Nk);
plot_band_kp(Energy, kdist, k_tick, kpts, labels)
MTB.plot.plot_bands_Electric(Energy',nbands,efermi,kdist,labels,k_tick,"kp",0)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       Construct FeSe 6 bands kp and get real space hoppings       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
syms M10 M11 M12 M20 M21 M22 A1 A2 B1 B2 dso real
syms kx ky kz real
% get the symbolic functions
% H_sym = FeSe_6bandkp_symbolic(kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso);
H_sym = FeSe_6bandtb_symbolic(kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso);
%  Define values
% % m10=1;m11=m10;m12=-0.5*m10;
% % m20=-m10;m21=-m11;m22=-m12;
% % a1=0.5;a2=0.1;b1=0;b2=0;dsoc=0.5;

%  Define values
m10=2;m11=1;m12=-0.5*m10;
m20=-m10;m21=-m11;m22=-m12;
a1=0.5;a2=0.1;b1=0;b2=0;dsoc=0.5;
%para in original code

% Substitute the symbolic parameters
params = {
    M10, m10;
    M11, m11;
    M12, m12;
    M20, m20;
    M21, m21;
    M22, m22;
    A1, a1;
    A2, a2;
    B1, b1;
    B2, b2;
    dso, dsoc;
};
H_fixed = subs(H_sym, params(:,1), params(:,2));
% H_func = matlabFunction(H_fixed,'Vars',{kx,ky,kz});

% === Define hopping directions (max: NN) ===
deltas = [
    0 0 0;   % on-site
    1 0 0;  
   -1 0 0;
    0 1 0;   
   0 -1 0;
    0 0 1;  
   0 0 -1
];

fprintf("Computing tight-binding hopping terms from H_sym(kx,ky,kz)...\n");

hoppings=get_hopping_inreal(H_fixed,deltas);

%%
% === Display results ===
fprintf("\nTight-binding hopping terms (numerical coefficients, symbolic parameters):\n");
for i = 1:length(hoppings)
    delta = hoppings{i,1};
    T = hoppings{i,2};
    fprintf("δ = [%2d, %2d, %2d]:\n", delta(1), delta(2),delta(3));
    disp(T);
end
%%
g = MTB.geometry("Fe_6band_tb");
g.a=[1,0,0;...
    0,1,0;...
    0,0,1];
g.b=inv(g.a)*2*pi;
g.atoms=[0,0,0;... %atom1
         0,0,0;... %atom2
         0,0,0;... %atom3
         0,0,0;...
         0,0,0;...
         0,0,0];    %atom4
g.wpos=g.atoms*g.a;
nhoprs=size(hoppings,1);
nbands=size(hoppings{1,2},1);
%
g.ham=zeros(nbands,nbands,nhoprs);
g.hopr=zeros(nhoprs,3);
for idx=1:nhoprs
    g.ham(:,:,idx) = double(hoppings{idx,2});
    g.hopr(idx,:) = hoppings{idx,1};
end

% Calculate bulk bands
% g=creat_wannier_from_kp(hoppings);
%%
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','Z'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.0,0.0,0.5]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kindex,"FeSe",Electric_field_in_evpA*10000);

%%
filename="./FeSe_6band_wannier90_par2.dat"
MTB.write_hr(g,filename)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Construct Kai Chern 3bands and get real space hoppings in 3D    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% === Define parameters and symbolic variables ===
syms tdd tpd tpp tpp_prime delta onsite_d real
syms kx ky kz real

H_sym = Kai_3bandChern_tb_symbolic(kx, ky, kz, tdd, tpd, tpp, tpp_prime,delta, onsite_d);

%  Define values
tdd_ = 1;
tpd_ = 1;
tpp_ = 1;
delta_=2.8;
tpp_prime_=tpp_*delta_/(4*tpp_+delta_);
onsite_d_=-4*tdd_+2*tpp_+delta_-2*tpp_*delta_/(4*tpp_+delta_);

%para in original code
% Substitute the symbolic parameters
params = {
    tdd, tdd_;
    tpd, tpd_;
    tpp, tpp_;
    tpp_prime, tpp_prime_;
    delta, delta_;
    onsite_d, onsite_d_;
};
H_fixed = subs(H_sym, params(:,1), params(:,2));

% === Define hopping directions (max: NN) ===
deltas = [
    0 0 0;   % on-site
    1 0 0;  
   -1 0 0;
    0 1 0;   
   0 -1 0;
];


fprintf("Computing tight-binding hopping terms from H_sym(kx,ky,kz)...\n");

hoppings=get_hopping_inreal(H_fixed,deltas);
%
g = MTB.geometry("Kai_3band_chern");
g.a=[1,0,0;...
    0,1,0;...
    0,0,1];
g.b=inv(g.a)*2*pi;
g.atoms=[0,0,0;... %atom1
         0,0,0;... %atom2
         0,0,0];   %atom3

g.wpos=g.atoms*g.a;
nhoprs=size(hoppings,1);
nbands=size(hoppings{1,2},1);
%
g.ham=zeros(nbands,nbands,nhoprs);
g.hopr=zeros(nhoprs,3);
for idx=1:nhoprs
    g.ham(:,:,idx) = double(hoppings{idx,2});
    g.hopr(idx,:) = hoppings{idx,1};
end

%%
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kindex,"FeSe",Electric_field_in_evpA*10000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Calculate the BC in xy plane            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=51;
kxline=[0,1];
kyline=[0,1];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[Unk,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
% Calculate Berry Curvature by LOOP method
plottap=2;
bandindex=2;
[Omega_k,KX,KY] = MTB.ham.get_Berry_curvature(bandindex,Unk,Kx,Ky,plottap);
Chern=sum(Omega_k,'all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Get the Wilson Loop              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=51;
kx=linspace(0,1,knum);
band1=2;
band2=2;
[wx1,unk]=MTB.ham.get_wilsonloop(g,knum,band1,band2);
[wx2,unk]=MTB.ham.get_wilsonloop_ky(g,knum,band1,band2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Construct Kai Chern 2bands and get real space hoppings in 3D    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
% === Define parameters and symbolic variables ===
syms t t1p t2p tpp phi real
syms kx ky kz real

% H_sym = Kai_3bandChern_tb_symbolic(kx, ky, kz, tdd, tpd, tpp, tpp_prime,delta, onsite_d);


H_sym = Kai_2bandChern_tb_symbolic(kx, ky, kz,t,t1p,t2p,tpp,phi)

%
%  Define values
t_ = -1;
t1p_ = t_/(2+sqrt(2));
t2p_ = -t1p_;
tpp_=t_/(2+2*sqrt(2));
phi_=pi/4;


%para in original code
% Substitute the symbolic parameters
params = {
    t, t_;
    t1p, t1p_;
    t2p, t2p_;
    tpp, tpp_;
    phi, phi_;
};

H_fixed = subs(H_sym, params(:,1), params(:,2));

% === Define hopping directions (max: NN) ===
deltas = [
    0 0 0;   % on-site
    1 0 0;  
   -1 0 0;
    0 1 0;   
   0 -1 0;
   1 1 0;
   -1 1 0;
   -1 -1 0;
   1 -1 0
];


fprintf("Computing tight-binding hopping terms from H_sym(kx,ky,kz)...\n");

hoppings=get_hopping_inreal(H_fixed,deltas);
 % hoppings=get_hopping_inreal(H_sym,deltas);
%
%
%
g = MTB.geometry("Kai_2band_chern");
g.a=[1,0,0;...
    0,1,0;...
    0,0,1];
g.b=inv(g.a)*2*pi;
g.atoms=[0,0,0;... %atom1
         0.5,0.5,0];   %atom2

g.wpos=g.atoms*g.a;
nhoprs=size(hoppings,1);
nbands=size(hoppings{1,2},1);
%
g.ham=zeros(nbands,nbands,nhoprs);
g.hopr=zeros(nhoprs,3);
for idx=1:nhoprs
    g.ham(:,:,idx) = double(hoppings{idx,2});
    g.hopr(idx,:) = hoppings{idx,1};
end

%
[nbands,~,nrpts]=size(g.ham);
labels={'X','\Gamma','M'}; % labels for k
hkpoints={[0.5,0.0,0.0],...
          [0.0,0.0,0.0],...
          [0.5,0.5,0.0]};% hkpoints-high symmetry k points
efermi=0.0;
nk=101;
Electric_field_in_evpA=0.0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kindex,"FeSe",Electric_field_in_evpA*10000);


%%
filename="./Kai_2band_Chern.dat"
MTB.write_hr(g,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Get slab bands for edges           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs=g;
MillerIndices=[0,1,0];
Umatrix=gs.MillerIndicestoumatrix(MillerIndices);
Urot=gs.surfab;
nslab=51;
%
[nbands,~,nrpts]=size(gs.ham);
labels={'X','\Gamma','X'};
hkpoints={[-0.5,0.0],...
          [0.0,0.0],...
          [0.5,0.0]};% hkpoints-high symmetry k points
nk=31;
efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(gs.ham,gs.hopr2,nslab,nbands,nrpts,hkpoints,nk,gs.a2,gs.b2);
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"TaIrTe4-slab",0);



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================= FUNCTIONS =============================%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       H_sym for FeSe 6band kp        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = FeSe_6bandkp_symbolic(kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso)
% 构造 6×6 符号 H(k) 哈密顿量，包含动量和材料参数
% Pauli matrices
s0 = eye(2);
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];

% Gamma matrices
Gamma1 = kron(sx, sz);
Gamma2 = kron(sy, s0);
Gamma3 = kron(sx, sx);
Gamma4 = kron(sx, sy);
Gamma5 = kron(sz, s0);

% Expressions
M1 = M10 + M11*(kx^2 + ky^2) + M12*kz^2;
M2 = M20 + M21*(kx^2 + ky^2) + M22*kz^2;
D  = B1*(kx^2 - ky^2) - 1i*B2*kx*ky;
kp = kx + 1i*ky;
km = kx - 1i*ky;

% Hamiltonian
H = [ M1,    0,  A2*kz,  -A1*km,  A1*kp,    0;
       0,   M1,  A1*kp,   A2*kz,     0, -A1*km;
    A2*kz, A1*km,   M2,      0,      0,  conj(D);
   -A1*kp, A2*kz,    0,     M2,      D,      0;
    A1*km,    0,     0,  conj(D), M2+dso,  0;
        0, -A1*kp,   D,     0,     0,   M2+dso];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       H_sym for FeSe 6band tb        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = FeSe_6bandtb_symbolic(kx, ky, kz, M10, M11, M12, M20, M21, M22, A1, A2, B1, B2, dso)
% 构造 6×6 符号 H(k) 哈密顿量，包含动量和材料参数
% Pauli matrices
s0 = eye(2);
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];

% Gamma matrices
Gamma1 = kron(sx, sz);
Gamma2 = kron(sy, s0);
Gamma3 = kron(sx, sx);
Gamma4 = kron(sx, sy);
Gamma5 = kron(sz, s0);

% Expressions
M1=M10+M11*(2*(1-cos(kx))+2*(1-cos(ky)))+M12*2*(1-cos(kz));
M2=M20+M21*(2*(1-cos(kx))+2*(1-cos(ky)))+M22*2*(1-cos(kz));
D=B1*(2*(1-cos(kx))-2*(1-cos(ky)))-1i*B2*sin(kx)*sin(ky);
kp=sin(kx)+1i*sin(ky);
km=sin(kx)-1i*sin(ky);

% Hamiltonian
H = [ M1,    0,  A2*kz,  -A1*km,  A1*kp,    0;
       0,   M1,  A1*kp,   A2*kz,     0, -A1*km;
    A2*kz, A1*km,   M2,      0,      0,  conj(D);
   -A1*kp, A2*kz,    0,     M2,      D,      0;
    A1*km,    0,     0,  conj(D), M2+dso,  0;
        0, -A1*kp,   D,     0,     0,   M2+dso];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       H_sym for BHZ 4band kp         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = bhz4band_tb_symbolic(kx, ky, kz, M, A, B)
% 构造 4×4 符号 H(k) 哈密顿量，包含动量和材料参数
s0 = eye(2);
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];

% Gamma matrices
Gamma1 = kron(sx, sz);
Gamma2 = kron(sy, s0);
Gamma3 = kron(sx, sx);
Gamma4 = kron(sx, sy);
Gamma5 = kron(sz, s0);

% H(k) = A*(sin(kx)Γ1 + sin(ky)Γ2) + [M + 2B(2 - cos(kx) - cos(ky))] * Γ5
H = A * sin(kx) * Gamma1 ...
   + A * sin(ky) * Gamma2 ...
   + (M + 2*B*(2 - cos(kx) - cos(ky))) * Gamma5 ...
   + kz*0;
% Pauli matrices
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       H_sym for 3band Chern tb        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = Kai_3bandChern_tb_symbolic(kx, ky, kz, tdd, tpd, tpp, tpp_prime,delta, onsite_d)

M11=-2*tdd*(cos(kx)+cos(ky))+onsite_d;
M22=2*tpp*cos(kx)-2*tpp_prime*cos(ky);
M33=2*tpp*cos(ky)-2*tpp_prime*cos(kx);

% === Define H_sym with lattice regularization ===
H = [         M11,      2*1j*tpd*sin(kx), 2*1j*tpd*sin(ky);...
         -2*1j*tpd*sin(kx),      M22,            1j*delta;     ...
         -2*1j*tpd*sin(ky),      -1j*delta,          M33];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       H_sym for Kai 2band Chern tb   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = Kai_2bandChern_tb_symbolic(kx, ky, kz,t,t1p,t2p,tpp,phi)

% phi=pi/4;
M1=(t1p+t2p)*(cos(kx)+cos(ky))+4*tpp*cos(kx)*cos(ky);
M12=4*t*cos(phi)*(cos(kx/2)*cos(ky/2))-1j*4*t*sin(phi)*(sin(kx/2)*sin(ky/2))
M21=4*t*cos(phi)*(cos(kx/2)*cos(ky/2))+1j*4*t*sin(phi)*(sin(kx/2)*sin(ky/2))

U=diag([1,exp(1j*(kx+ky)/2)])
M2=(t1p-t2p)*(cos(kx)-cos(ky))

H=[M1+M2,M12;...
   M21,M1-M2]
H=U*H*U'
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       Calculate bands for kp H_func(kx,ky,kz)  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ebands, kdist, k_tick] = get_bands_kp(H_func, kpts, labels, Nk)
% 计算并绘制3D高对称路径下的能带图

% === Interpolate k-path ===
klist = [];
for i = 1:size(kpts,1)-1
    seg_kx = linspace(kpts(i,1), kpts(i+1,1), Nk)';
    seg_ky = linspace(kpts(i,2), kpts(i+1,2), Nk)';
    seg_kz = linspace(kpts(i,3), kpts(i+1,3), Nk)';
    klist = [klist; [seg_kx, seg_ky, seg_kz]];
end

% === k distance ===
kdist = [0; cumsum(sqrt(sum(diff(klist).^2, 2)))];

% Compute x-ticks (positions of high-symmetry points)
k_tick = [0];
for i = 2:size(kpts,1)
    d = norm(kpts(i,:) - kpts(i-1,:));
    k_tick(end+1) = k_tick(end) + d-10^-6;
end

% === Compute band structure ===
nBands = size(H_func(0,0,0),1);
Ebands = zeros(length(klist), nBands);
for i = 1:length(klist)
    kx = klist(i,1); ky = klist(i,2); kz = klist(i,3);
    H = H_func(kx, ky, kz);
    Ebands(i,:) = sort(real(eig(H)))';
end
end

function plot_band_kp(Ebands,kdist,k_tick,kpts,labels)
figure;
plot(kdist, Ebands, 'LineWidth', 1.2);

set(gca, 'XTick', k_tick);
set(gca, 'XTickLabel', labels, 'FontSize', 12);
xline(kdist(1), '--k');
for i = 2:size(kpts,1)-1
    xline(k_tick(i), '--k');
end
xlim([kdist(1), kdist(end)]);
grid on;
title('Band structure along high-symmetry path');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              FFT to real space hoppings        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hoppings=get_hopping_inreal(H_sym,deltas)
syms kx ky kz real
hoppings = {};
for i = 1:size(deltas,1)
    dx = deltas(i,1);
    dy = deltas(i,2);
    dz = deltas(i,3);
    % Phase factor: e^{-i(k·delta)}
    phase = exp(-1i * (kx*dx + ky*dy + kz*dz));

    % Inverse Fourier transform
    Tdelta = int(int(int(H_sym * phase, kx, -pi, pi), ky, -pi, pi), kz, -pi, pi) / (2*pi)^3;

    % Simplify & store
    Tdelta = vpa(simplify(Tdelta, 'Steps', 50), 6);
    hoppings{end+1,1} = [dx, dy, dz];
    hoppings{end,2} = Tdelta;
end

% === Display results ===
fprintf("\nTight-binding hopping terms (numerical coefficients, symbolic parameters):\n");
for i = 1:length(hoppings)
    delta = hoppings{i,1};
    T = hoppings{i,2};
    fprintf("δ = [%2d, %2d, %2d]:\n", delta(1), delta(2),delta(3));
    disp(T);
end
end

%%
function g=creat_wannier_from_kp(hoppings)
g = MTB.geometry("BBH");
g.a=[1,0,0;...
     0,1,0;...
     0,0,1];
g.b=inv(g.a)*2*pi;
nbands=size(hoppings{1,2},1);
nrpts=size(hoppings,1);
g.atoms=kron(ones(nbands,1),[0,0,0]);
g.wpos=g.atoms*g.a;
g.ham=zeros(nbands,nbands,nrpts);
g.hopr=zeros(nrpts,3);
for idx=1:nrpts
    g.ham(:,:,idx) = hoppings{idx,2};
    g.hopr(idx,:) =  hoppings{idx,1};
end
end