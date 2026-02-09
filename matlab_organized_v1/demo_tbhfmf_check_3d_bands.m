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
g.Rcart = g.hopr * g.a;
%%
% ---- choose patch ----
K0_frac = [1/3, 2/3];
L_frac  = [0.04, 0.04];
Nk      = [128, 128];         % [Nkx, Nky]

mesh = tbHFMF.build_kmesh_patch(g, K0_frac, L_frac, Nk);

% ---- build hk on 2D patch ----
hk = tbHFMF.build_hk_atomgauge(g, mesh);

% ---- check hk builder consistency ----
% tbHFMF.check_build_hk_atomgauge(g, mesh, hk, 12);

% ---- diagonalize: choose active bands ----
% 如果你想先看"完整面"，就把 act_idx = 1:size(g.ham,1)
act_idx = 9:12;
act = tbHFMF.diag_hk_active(hk, act_idx);

band_ids=1:length(act_idx);
% ---- plot one band surface ----
% plot_band_surface_patch(mesh, act, band_ids, 'imagesc');  % lowest band
plot_band_surface_patch(mesh, act, band_ids, 'surf');    % second band
%%


check_build_hk_atomgauge(g, mesh, hk, 4)
%%
function plot_band_surface_patch(mesh, act, band_ids, mode)
%PLOT_BAND_SURFACE_PATCH  Plot E_n(kx,ky) on the patch for one band.
%
% Inputs:
%   mesh    : from build_kmesh_patch
%   act     : from diag_hk_active
%   band_id : 1..Nact (index within active set)
%   mode    : 'imagesc' or 'surf' (default 'imagesc')

if nargin < 4, mode = 'imagesc'; end

Nky = mesh.Nky; Nkx = mesh.Nkx;
figure();
hold on;

for i =1:length(band_ids)
band_id=band_ids(i);
E = reshape(act.eps_flat(band_id,:), Nky, Nkx);   % centered order
Px = mesh.Px;  Py = mesh.Py;                      % Å^-1 (local momentum)

switch lower(mode)
  case 'imagesc'
    % figure;
    imagesc(Px(1,:), Py(:,1), E);   % x uses row of Px, y uses col of Py
    axis xy; axis image;
    xlabel('p_x (Å^{-1})'); ylabel('p_y (Å^{-1})');
    title(sprintf('Active band %d: E(p) on patch', band_id));
    colorbar;

  case 'surf'
    % figure;
    surf(Px, Py, E, 'EdgeColor','none');
    xlabel('p_x (Å^{-1})'); ylabel('p_y (Å^{-1})'); zlabel('E (eV)');
    title(sprintf('Active band %d: E(p) on patch', band_id));
    view(35,35); box on;

  otherwise
    error('mode must be imagesc or surf');
end
axis tight
end
end

function rep = check_build_hk_atomgauge(g, mesh, hk, ntest)
%CHECK_BUILD_HK_ATOMGAUGE  Validate build_hk_atomgauge against single-k builder.
%
% Checks:
%   (i) max Hermiticity error of hk over random k points
%  (ii) max matrix difference between hk(:,:,:, :) and get_hk_atomgauge_single at same k

if nargin < 4, ntest = 10; end

[Norb,~,Nky,Nkx] = size(hk);
Nk = Nky*Nkx;

% flatten index mapping: (iy,ix) <-> ik
% ik = (ix-1)*Nky + iy  if we reshape with (Nky,Nkx) order
% but easiest: just pick (iy,ix) directly.
rng(1);
ixs = randi(Nkx, [ntest,1]);
iys = randi(Nky, [ntest,1]);

maxHerm = 0;
maxDiff = 0;

for t = 1:ntest
  ix = ixs(t); iy = iys(t);

  % hk from bulk builder
  H1 = hk(:,:,iy,ix);
  herm = norm(H1 - H1','fro') / max(1, norm(H1,'fro'));
  maxHerm = max(maxHerm, herm);

  % same k (reduced) for single-k builder
  k_frac = [mesh.KXf(iy,ix), mesh.KYf(iy,ix), 0]; % 2D patch -> kz=0
  H2 = tbHFMF.get_hk_atomgauge_single(g, k_frac);

  diff = norm(H1 - H2,'fro') / max(1, norm(H2,'fro'));
  maxDiff = max(maxDiff, diff);
end

rep = struct();
rep.maxHerm_rel = maxHerm;
rep.maxDiff_rel = maxDiff;

fprintf('[check_build_hk_atomgauge] max rel herm err  = %.3e\n', rep.maxHerm_rel);
fprintf('[check_build_hk_atomgauge] max rel H diff     = %.3e\n', rep.maxDiff_rel);

end

