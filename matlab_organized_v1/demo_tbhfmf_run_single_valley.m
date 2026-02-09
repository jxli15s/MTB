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

% g = MTB.geometry("Rgra_15s");
% g = MTB.read_poscar(g,"data/Graphene/15s/fplo/encut_25/wannier90_formula//POSCAR");
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/15s/fplo/encut_25/wannier90_formula/wannier90_hr_p1.dat','data/Graphene/15s/fplo/encut_25/wannier90_formula/wannier90_hr_p2.dat');
%%
g.wpos=g.atoms*g.a;
g.wpos(:,3)=g.wpos(:,3)-mean(g.wpos(:,3));
% g.wpos=g.wpos*0
g.Rcart = g.hopr * g.a;
%%
[out,opts] = main_bandproj_LAF_seeded(g);
%%

%%
ek=reshape(out.enk,20,301,301);
figure()
hold on;

for ik=151:151
for i =9:12
    plot(squeeze(ek(i,ik,:)))
end
end

figure()
hold on;
for i=9:12
    plot(squeeze(out.enk(i,:)))
end
%%
figure()
hold on;
for i =9:12
    surf(squeeze(ek(i,:,:)),'EdgeColor','none')
end
%%
unk=reshape(out.unk,20,20,101,101);
%%
unk_10=unk(:,10,51,51);
unk10=unk_10.*conj(unk_10)
unk_9=unk(:,9,51,51);
unk9=unk_9.*conj(unk_9)

%%
function [out,opts] = main_bandproj_LAF_seeded(g)
%MAIN_BANDPROJ_LAF_SEEDED  Band-projected HF with decaying LAF seed.
%
% Requirements (already in your +tbHFMF, not printed here):
%   tbHFMF.HFMF_BandProj
%   tbHFMF.build_Vq_doublegate_layered
%   tbHFMF.fock_fft_atomgauge_flat
%   tbHFMF.hartree_orbital_flat
%   tbHFMF.precompute_dtau_phases
%   tbHFMF.build_kmesh_patch, build_hk_atomgauge, diag_hk_active, etc.

%% ---------- user knobs ----------
opts = struct();

% Patch mesh (reduced coords)
opts.K0_frac = [1/3, 2/3];        % valley center in reduced coords
opts.L_frac  = [0.15, 0.15];      % patch size in reduced coords
opts.Nk      = [301, 301];          % [Nkx, Nky]

% Active bands (indices AFTER sorting eigenvalues ascending)
opts.act_idx  = 1:20;

% Temperature + SCF
opts.kT       = 1e-3;             % eV
opts.max_iter = 100;
opts.mix      = 0.7;
opts.tol      = 1e-6;
opts.verbose  = true;

% Hartree toggle
opts.use_hartree = true;
opts.hartree_subtract_mean = true;

% Coulomb kernel builder (double gate)
opts.Vbuilder = @tbHFMF.build_Vq_doublegate_layered;
eps_r = 10;                       % relative dielectric (example)
opts.Vpars = struct();
opts.Vpars.d_gate      = 400;      % Å  (example)
opts.Vpars.e2_over_eps = 14.3996/eps_r;  % eV·Å
opts.Vpars.q_small     = 1e-10;     % Å^-1
opts.Vpars.include_q2  = false;
opts.Vpars.z0_mode     = 'center';

% ----- Seed (LAF) -----
opts.Delta_seed = 1e-3;            % eV (2 meV)
opts.lambda0    = 1.0;
opts.seed_eta   = 0.80;
opts.lambda_min = 1e-3;

% Filling mode (simplest): target average occupancy per active band
%   < sum_n f_n(k) >_k = nu * Nact
opts.nu = 0.50;                    % 0.5 filling within active subspace

%% ---------- ensure wpos ready ----------
% you said you already do this; keep it here as safety
g.wpos = g.atoms * g.a;
g.wpos(:,3) = g.wpos(:,3) - mean(g.wpos(:,3));

%% ---------- build band-proj object (uses your existing tbHFMF code) ----------
obj = tbHFMF.HFMF_BandProj(g, opts);

% We will override obj.rho_band_flat using seeded initialization:
Uflat   = obj.act.U_flat;      % Norb x Nact x Nk
epsflat = obj.act.eps_flat;    % Nact x Nk

%% ---------- build seed in orbital basis, then project to band basis ----------
Hseed_orb = seed_LAF_orb_from_z_spinblock(g, opts.Delta_seed);

seedBand = project_seedBand_from_Uflat(Uflat, Hseed_orb); % Nact x Nact x Nk

%% ---------- initialize rho_band from seeded band Hamiltonian ----------
rho0 = init_rho_band_from_seedBand(epsflat, seedBand, opts);

obj.rho_band_flat = rho0;

%% ---------- run SCF with decaying seed ----------
out = run_bandproj_decaying_seed(obj, seedBand, opts);
out.mesh=obj.mesh;

end

function Hseed = seed_LAF_orb_from_z_spinblock(g, Delta)
%SEED_LAF_ORB_FROM_Z_SPINBLOCK  Orbital-basis LAF seed:
%   Hseed = Delta * diag( sz .* (z/max|z|) )
% Assumes spin block ordering: [up block; down block].

Norb = size(g.wpos,1);
assert(mod(Norb,2)==0, 'Spin-block assumption requires even Norb.');

z = g.wpos(:,3);
z = z - mean(z);

zmax = max(abs(z));
if zmax < 1e-12
  error('z variation too small; cannot form layer seed.');
end
Lz = z / zmax;                               % [-1,1], top->bottom

sz = [ones(Norb/2,1); -ones(Norb/2,1)];      % spin block
op = sz .* Lz;
% sz = kron(ones(15,1),[1,1]');
% op = sz .* Lz;

op = op - mean(op);                           % remove uniform offset (recommended)

Hseed = Delta * diag(op);
Hseed = (Hseed + Hseed')/2;
end

function seedBand = project_seedBand_from_Uflat(Uflat, Hseed_orb)
%PROJECT_SEEDBAND_FROM_UFLAT  seedBand(k)=U(k)'*Hseed_orb*U(k)
% Uflat : Norb x Nact x Nk
% output: Nact x Nact x Nk

[Norb,Nact,Nk] = size(Uflat);
seedBand = zeros(Nact,Nact,Nk,'like',1+1i);

Hseed_orb = (Hseed_orb + Hseed_orb')/2;

parfor ik = 1:Nk
  U = Uflat(:,:,ik);
  S = U' * Hseed_orb * U;
  seedBand(:,:,ik) = (S + S')/2;
end
end


function rho_band = init_rho_band_from_seedBand(epsflat, seedBand, opts)
%INIT_RHO_BAND_FROM_SEEDBAND  Initialize rho_band from
%   Heff0(k)=diag(epsflat(:,k)) + lambda0*seedBand(:,:,k)
% and solve mu by target nu filling in active subspace.

[Nact,Nk] = size(epsflat);
lam0 = opts.lambda0;

Ek = zeros(Nact,Nk);
Wk = zeros(Nact,Nact,Nk,'like',1+1i);

parfor ik = 1:Nk
  H = diag(epsflat(:,ik)) + lam0*seedBand(:,:,ik);
  H = (H+H')/2;
  [W,D] = eig(H,'vector');
  [e,ord] = sort(real(D),'ascend');
  Wk(:,:,ik) = W(:,ord);
  Ek(:,ik)   = e(ord);
end

mu = solve_mu_nu_simple(Ek, opts.kT, opts.nu);

rho_band = zeros(Nact,Nact,Nk,'like',1+1i);
parfor ik = 1:Nk
  W = Wk(:,:,ik);
  f = fermi_vec(Ek(:,ik), mu, opts.kT);
  R = W * (f .* (W'));                % W*diag(f)*W'
  rho_band(:,:,ik) = (R + R')/2;
end
end

function out = run_bandproj_decaying_seed(obj, seedBand, opts)
%RUN_BANDPROJ_DECAYING_SEED  Same logic as HFMF_BandProj.run, but with:
%   Heff(k) = diag(eps) + Sigma_band(k) + lambda(it)*seedBand(k)
% and mu solved by simple nu filling (no area/target complications).

mesh = obj.mesh;
act  = obj.act;
V    = obj.V;
ph   = obj.ph;

Norb = size(obj.g.ham,1);
Nact = act.Nact;
Nk   = mesh.Nk;

Uflat   = act.U_flat;      % Norb x Nact x Nk
epsflat = act.eps_flat;    % Nact x Nk

rho_band_flat = obj.rho_band_flat;

hist = struct();
hist.it     = [];
hist.err    = [];
hist.mu     = [];
hist.lambda = [];

for it = 1:opts.max_iter

  % seed schedule (no IF inside parfor)
  lam = opts.lambda0 * (opts.seed_eta^(it-1));
  if lam < opts.lambda_min, lam = 0; end

  % (1) rho_band -> rho_orb
  rho_orb_flat = zeros(Norb,Norb,Nk,'like',1+1i);
  parfor ik = 1:Nk
    U  = Uflat(:,:,ik);
    rb = rho_band_flat(:,:,ik);
    rho_orb_flat(:,:,ik) = U * rb * U';
  end

  % (2) Fock
  SigmaF_orb_flat = tbHFMF.fock_fft_atomgauge_flat(mesh, V, ph, rho_orb_flat);

  % (3) Hartree
  if isfield(opts,'use_hartree') && opts.use_hartree
    SigmaH = tbHFMF.hartree_orbital_flat(V.V0, rho_orb_flat, mesh, opts);
  else
    SigmaH = zeros(Norb,Norb,'like',1+1i);
  end

  Sigma_orb_flat = SigmaF_orb_flat + SigmaH;

  % (4) Project Σ_orb -> Σ_band
  Sigma_band_flat = zeros(Nact,Nact,Nk,'like',1+1i);
  parfor ik = 1:Nk
    U = Uflat(:,:,ik);
    S = U' * Sigma_orb_flat(:,:,ik) * U;
    Sigma_band_flat(:,:,ik) = (S + S')/2;
  end

  % (5) Diagonalize Heff in active space WITH seed
  Ek_flat = zeros(Nact,Nk);
  W_flat  = zeros(Nact,Nact,Nk,'like',1+1i);
  ham_flat = zeros(Nact,Nact,Nk,'like',1+1i);
  parfor ik = 1:Nk
    Heff = diag(epsflat(:,ik)) + Sigma_band_flat(:,:,ik) + lam*seedBand(:,:,ik);
    Heff = (Heff + Heff')/2;
    ham_flat(:,:,ik) = Heff;
    [W,D] = eig(Heff,'vector');
    [E,ord] = sort(real(D),'ascend');
    W = W(:,ord);

    Ek_flat(:,ik)  = E;
    W_flat(:,:,ik) = W;
  end

  % (6) μ by nu filling (simple)
  mu = solve_mu_nu_simple(Ek_flat, opts.kT, opts.nu);

  % (7) Update rho_band
  rho_new_flat = zeros(Nact,Nact,Nk,'like',1+1i);
  parfor ik = 1:Nk
    W = W_flat(:,:,ik);
    f = fermi_vec(Ek_flat(:,ik), mu, opts.kT);
    R = W * (f .* (W'));              % W*diag(f)*W'
    rho_new_flat(:,:,ik) = (R + R')/2;
  end

  % mixing + convergence
  rho_mixed = (1-opts.mix)*rho_band_flat + opts.mix*rho_new_flat;
  err = max(abs(rho_mixed(:) - rho_band_flat(:)));

  rho_band_flat = rho_mixed;

  hist.it(end+1)     = it; %#ok<AGROW>
  hist.err(end+1)    = err;
  hist.mu(end+1)     = mu;
  hist.lambda(end+1) = lam;

  if isfield(opts,'verbose') && opts.verbose
    fprintf('it=%3d  err=%.3e  mu=%.6f  lambda=%.2e\n', it, err, mu, lam);
  end

  if err < opts.tol
    break;
  end
end

obj.rho_band_flat = rho_band_flat;
obj.mu = mu;

out = struct();
out.rho_band_flat = rho_band_flat;
out.mu = mu;
out.history = hist;
out.enk = Ek_flat;
out.unk = W_flat;
out.ham = ham_flat;
end

function mu = solve_mu_nu_simple(Ek, kT, nu)
%SOLVE_MU_NU_SIMPLE  Solve mu for target nu filling in active subspace:
%   mean_k sum_n f(E_n(k)-mu) = nu * Nact
%
% Ek : Nact x Nk

E = real(Ek);
[Nact,~] = size(E);
target = nu * Nact;

g = @(mu) mean(sum(fermi_stable(E,mu,kT),1));

emin = min(E,[],'all');
emax = max(E,[],'all');
W = 50*max(kT,1e-6) + 1;
muL = emin - W;
muR = emax + W;

gL = g(muL);
gR = g(muR);
if target <= gL, mu = muL; return; end
if target >= gR, mu = muR; return; end

for it = 1:80
  muM = 0.5*(muL+muR);
  gM  = g(muM);
  if gM > target
    muR = muM;
  else
    muL = muM;
  end
  if abs(muR-muL) < 1e-12, break; end
end
mu = 0.5*(muL+muR);
end

function F = fermi_stable(E, mu, kT)
if kT <= 0
  F = double(E < mu);
  return
end
x = (E - mu)./kT;
x = max(min(x,50),-50);
F = 1./(1+exp(x));
end

function f = fermi_vec(e, mu, kT)
e = real(e(:));
if kT <= 0
  f = double(e < mu);
  return
end
x = (e - mu)./kT;
x = max(min(x,50),-50);
f = 1./(1+exp(x));
end

