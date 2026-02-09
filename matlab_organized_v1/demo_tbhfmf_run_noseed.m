%%
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
g.wpos(:,3)=g.wpos(:,3)-mean(g.wpos(:,3));
% g.wpos=g.wpos*0
g.Rcart = g.hopr * g.a;
%%
Delta_seed=0.005
tol = 1e-3; % Set tolerance for layer identification
Hseed = seed_LAF_from_wpos_z_spinblock(g, Delta_seed)
%
kfrac=[1/3,2/3,0];
Hk = tbHFMF.get_hk_atomgauge_single(g, kfrac);
Hk_seed=Hk+Hseed;
[e,u]=tbHFMF.solve_hk(Hk_seed)

%%
% Options
opts = struct();
opts.Nk      = [101, 101];        % (Nkx,Nky)
opts.K0_frac = [1/3, 2/3];      % patch center in reduced coords (Γ). Try [1/3,2/3] if you like.
opts.L_frac  = [0.2, 0.2];    % patch size in reduced coords
opts.kT      = 5e-3;            % eV
opts.max_iter = 200;
opts.tol      = 1e-8;
opts.mix      = 0.35;
opts.verbose  = true;

% Active-band projection: keep only the lowest band
opts.act_idx = 9:12;

mesh = tbHFMF.build_kmesh_patch(g, opts.K0_frac, opts.L_frac, opts.Nk);

% Filling constraint on the patch (per unit cell contribution from this patch):
% For Nact=1: n_cell = area_frac * <f>_k. Choose target <f>=0.50.
Nact = numel(opts.act_idx);
opts.n_target_cell = mesh.area_frac*(0.50 * Nact);

%
% Hartree (optional)
opts.use_hartree = true;
opts.hartree_subtract_mean = true;

% Double-gate Coulomb parameters (toy)
eps_r = 2;
opts.Vbuilder = @tbHFMF.build_Vq_doublegate_layered;
opts.Vpars = struct();
opts.Vpars.d_gate = 40.0;                 % Å
opts.Vpars.e2_over_eps = 14.3996/eps_r;   % eV·Å
opts.Vpars.q_small = 1e-6;                % Å^-1
opts.Vpars.include_q2 = false;
opts.Vpars.z0_mode = 'center';

% Run
solver = tbHFMF.HFMF_BandProj(g, opts);
out = solver.run();

%%
% ---- build hk on 2D patch ----
hk = tbHFMF.build_hk_atomgauge(g, mesh);
act = tbHFMF.diag_hk_active(hk, opts.act_idx);
mu=tbHFMF.solve_mu_fill_flat(act.eps_flat, mesh, opts);
%%
ek=reshape(out.ek_flat,20,51,51)
figure()
hold on;
surf(squeeze(ek(10,:,:)))
surf(squeeze(ek(11,:,:)))
%%

figure()
for i=1:4
    hold on;
plot(out.ek_flat(i,:))
end


%%
function Hseed = seed_LAF_from_wpos_z_spinblock(g, Delta)
%LAF seed: Hseed = Delta * diag( sz .* Lz )
% - g.wpos: Norb x 3, already contains spin DOF (duplicated positions)
% - spin block order assumed: first half = up, second half = down
% - Lz is taken directly from centered z and normalized to [-1,1]

Norb = size(g.wpos,1);
assert(mod(Norb,2)==0, 'Norb must be even for spin-block assumption.');

% layer operator from z (top->bottom decreasing)
z  = g.wpos(:,3);
z  = z - mean(z);                 % safety (you already did this)
zmax = max(abs(z));
assert(zmax > 1e-12, 'z variation too small.');

Lz = z / zmax;                    % in [-1,1]

% spin operator sz for spin-block basis
sz = [ones(Norb/2,1); -ones(Norb/2,1)];

% (optional) quick sanity check: up/down should share same z distribution
dz = norm(z(1:Norb/2) - z(Norb/2+1:end)) / max(1,norm(z(1:Norb/2)));
fprintf('[seed] up/down z mismatch ratio = %.3e (should be ~0 for perfect spin-block)\n', dz);

Hseed = Delta * diag(sz .* Lz);
Hseed = (Hseed + Hseed')/2;
end



