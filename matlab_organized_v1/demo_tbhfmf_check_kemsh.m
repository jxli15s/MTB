%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Construct the g.ham                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g = MTB.geometry("Rgra_5s");
g = MTB.read_poscar(g,"data/Graphene/5s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/5s/wannier90_hr_p1.dat','data/Graphene/5s/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;

%%
mesh=test_build_kmesh_patch_with_g(g);

%%
function mesh=test_build_kmesh_patch_with_g(g)
%TEST_BUILD_KMESH_PATCH_WITH_G  Use your real g to validate build_kmesh_patch.

% -------- user choices --------
Nk      = [129, 129];         % [Nkx, Nky]  推荐偶数便于 FFT
L_frac  = [0.2, 0.2];       % patch size in reduced coords (fraction of b1/b2)
K0_frac = [1/3, 2/3];         % valley center in reduced coords (e.g. graphene K)
% 如果你要 Γ patch：K0_frac=[0,0];

% -------- build mesh --------
mesh = tbHFMF.build_kmesh_patch(g, K0_frac, L_frac, Nk);

% -------- run your check --------
check_kmesh(mesh);

% -------- extra sanity checks --------

% (A) Check that center k equals K0 in Cartesian too
Nkx = mesh.Nkx; Nky = mesh.Nky;
icx = floor(Nkx/2)+1;
icy = floor(Nky/2)+1;

k0_cart = K0_frac(1)*g.b(1,:) + K0_frac(2)*g.b(2,:);
k_center = [mesh.Kx(icy,icx), mesh.Ky(icy,icx), mesh.Kz(icy,icx)];

fprintf('\n[Extra A] ||k_center - k0_cart|| = %.3e (Å^-1)\n', norm(k_center - k0_cart));

% (B) Check that qabs is really |p| (since q grid is p grid on patch)
p_center = [mesh.Px(icy,icx), mesh.Py(icy,icx), mesh.Pz(icy,icx)];
fprintf('[Extra B] ||p_center|| = %.3e, qabs(center)=%.3e\n', norm(p_center), mesh.qabs(icy,icx));

% (C) Check step size in reduced coords (should be L/N)
dx = mesh.PXf(icy,icx+1) - mesh.PXf(icy,icx);
dy = mesh.PYf(icy+1,icx) - mesh.PYf(icy,icx);
fprintf('[Extra C] dPXf=%.6e (target %.6e), dPYf=%.6e (target %.6e)\n', ...
  dx, L_frac(1)/Nk(1), dy, L_frac(2)/Nk(2));

% (D) Visual check: qabs minimum should be at center
[~,idxmin] = min(mesh.qabs(:));
[iymin,ixmin] = ind2sub(size(mesh.qabs), idxmin);
fprintf('[Extra D] qabs min at (%d,%d), center=(%d,%d)\n', iymin, ixmin, icy, icx);

end


function check_kmesh(mesh)
Nkx = mesh.Nkx; Nky = mesh.Nky;
icx = floor(Nkx/2)+1;
icy = floor(Nky/2)+1;

fprintf('Nkx=%d Nky=%d, center index = (%d,%d)\n', Nkx, Nky, icy, icx);

% (1) p=0 at center
fprintf('PXf(center)=%.3e, PYf(center)=%.3e\n', mesh.PXf(icy,icx), mesh.PYf(icy,icx));

% (2) k=K0 at center (reduced coord)
fprintf('KXf(center)=%.12f (target %.12f)\n', mesh.KXf(icy,icx), mesh.K0_frac(1));
fprintf('KYf(center)=%.12f (target %.12f)\n', mesh.KYf(icy,icx), mesh.K0_frac(2));

% (3) qabs(center)=0
fprintf('qabs(center)=%.3e (should be ~0)\n', mesh.qabs(icy,icx));

% (4) shift sanity: marker at center should move to (1,1)
A = zeros(Nky,Nkx); A(icy,icx)=1;
A_to = mesh.to_fft(A);
[i0,j0] = find(A_to==1);
fprintf('to_fft moves center marker -> (%d,%d) (should be (1,1))\n', i0, j0);

% (5) BZ area identity check: ABZ ?= (2pi)^2/Acell
lhs = mesh.ABZ;
rhs = (2*pi)^2 / mesh.Acell;
fprintf('ABZ = %.12e, (2pi)^2/Acell = %.12e, rel.err=%.3e\n', lhs, rhs, abs(lhs-rhs)/rhs);

% (6) pref_patch check
fprintf('area_frac=%.6f, Acell=%.6f Å^2, pref_patch=%.6e Å^-2\n', ...
  mesh.area_frac, mesh.Acell, mesh.pref_patch);
end

