% clear;
% clear all;
g = MTB.geometry("Rgra_5s");
g = MTB.read_poscar(g,"data/Graphene/5s/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/Graphene/5s/wannier90_hr_p1.dat','data/Graphene/5s/wannier90_hr_p2.dat');
g.wpos=g.atoms*g.a;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Calculate the band structures                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.wpos=g.atoms*g.a;
[nbands,~,nrpts]=size(g.ham);
labels={'\Gamma','K','M'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [2/3,1/3,0.0],...
%           [0.5,0.5,0.0]};% hkpoints-high symmetry k points
hkpoints={[2/3,1/3,0.0]*0.9,...
          [2/3,1/3,0.0],...
          [2/3,1/3,0.0]+([0.5,0.5,0.0]-[2/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
hkpoints={[1/3,2/3,0.0]*0.9,...
          [1/3,2/3,0.0],...
          [1/3,2/3,0.0]+([0.5,0.5,0.0]-[1/3,2/3,0.0])*0.15,...
          };% hkpoints-high symmetry k points
nk=251;
% efermi=-0.117;

Electric_field_in_evpA=0.00; %0.08-0.12 V/A
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands_add_electric(g.ham,g.hopr,g.wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,g.a,g.b);
MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"bulk",0)

%%
% ===== Example driver: Non-uniform (patch) HFMF without FFT =====


% ---- 0) Load or define your Wannier model 'g' ----
% Required fields:
%   g.a     (3x3) real-space lattice vectors (Å)
%   g.b     (3x3) reciprocal lattice vectors (Å^-1), g.b * g.a = 2*pi*I
%   g.hopr  (Lx3) hopping lattice vectors in fractional coords
%   g.ham   (m x m x L) hopping matrices H(R)
%   g.atoms (Nat x 3) atom positions in fractional coords   <-- per your note
%   g.wpos  (m x 3)  orbital positions in Cartesian (Å)     <-- per your note
% load('g.mat');

% ---- 1) Build non-uniform hex k-mesh ([0,1) frac) ----
b1 = g.b(1,1:2); b2 = g.b(2,1:2);
Ncoarse = 10; Ndense = 31; R_dense = 0.08;
% [k_cart2D, k_frac2D, w_k, region_id] = make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense);
[k_cart2D, k_frac2D, w_k, region_id] = make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense);
klist_frac = [k_frac2D, zeros(size(k_frac2D,1),1)]; % 2D->3D compatible
klist_cart = klist_frac*g.b;

% % region23=~(region_id==1);
% % klist_frac_23=klist_frac(region23,:);
% % figure()
% % hold on
% % plot(klist_frac(:,1),klist_frac(:,2),'ro')
% % plot(klist_frac_23(:,1),klist_frac_23(:,2),'bo')
% region23=~(region_id==1);
% klist_frac=klist_frac(region23,:);
% w_k=w_k(region23);

%
% ---- 2) HF options ----
opts = struct;
opts.nspin = 1;
opts.kT    = 0;
opts.vmodel= '2D';              % or 'keldysh'
opts.eps   = 25.0;
opts.e2    = 14.3996;
opts.r0    = 10.0;
opts.hartree_enable = true;
opts.fock_drop_q0   = true;
opts.qeps  = 1e-6;
opts.mixing= 0.6;
opts.tol   = 1e-6;
opts.max_iter = 100;
opts.verbose = true;
opts.Efilling = size(g.ham,1)*opts.nspin/2;

% ---- 3) Run HF (direct-sum Fock on patch mesh) ----
out = hf_run_wannier_patch(g, opts, klist_frac, w_k);

% ---- 4) Quick plot ----
Eall = out.E_k(:);
figure; histogram(Eall, 60); xlabel('E (eV)'); ylabel('count'); title('HF eigenvalues (all k/bands)');

%%
figure()
for i=1:20
    hold on;
plot(out.E_k(:,i))
end
%%
region23=~(region_id==1);
klist_frac_23=klist_frac(region23,:);
figure()
hold on
plot(klist_frac(:,1),klist_frac(:,2),'ro')
plot(klist_frac_23(:,1),klist_frac_23(:,2),'bo')

%%
nk=Ndense+1;
totalk=nk^2*2;
E_k_dense=out.E_k(region23,:);
klist_frac_dense=klist_frac(region23,:);
klist_frac_mesh1x=reshape(klist_frac_dense(1:totalk/2,1),nk,nk);
klist_frac_mesh1y=reshape(klist_frac_dense(1:totalk/2,2),nk,nk);
klist_frac_mesh2x=reshape(klist_frac_dense(totalk/2+1:end,1),nk,nk);
klist_frac_mesh2y=reshape(klist_frac_dense(totalk/2+1:end,2),nk,nk);

Enk1=reshape(E_k_dense(1:totalk/2,:),nk,nk,[]);
[E_sorted1, idx] = sort(Enk1, 3);   % idx: (Nk, Nk, Nband)
Enk2=reshape(E_k_dense(totalk/2+1:end,:),nk,nk,[]);
[E_sorted2, idx] = sort(Enk2, 3);   % idx: (Nk, Nk, Nband)
figure()
hold on
plot(klist_frac(:,1),klist_frac(:,2),'ro')
plot(klist_frac_dense(1:totalk/2,1),klist_frac_dense(1:totalk/2,2),'bo')
plot(klist_frac_dense(totalk/2+1:end,1),klist_frac_dense(totalk/2+1:end,2),'bo')
plot(klist_frac_mesh1x(6,:),klist_frac_mesh1y(6,:),'ko');
plot(klist_frac_mesh2x(1,:),klist_frac_mesh2y(1,:));
%%
figure()
hold on;
for j=6:6
for i =9:12
   plot(E_sorted1(j,:,i))
   plot(E_sorted2(j,:,i))
end
end
%%
figure()
hold on;
surf(klist_frac_mesh1x,klist_frac_mesh1y,E_sorted1(:,:,9),'EdgeColor','none')
surf(klist_frac_mesh1x,klist_frac_mesh1y,E_sorted1(:,:,10),'EdgeColor','none')
surf(klist_frac_mesh1x,klist_frac_mesh1y,E_sorted1(:,:,11),'EdgeColor','none')
surf(klist_frac_mesh1x,klist_frac_mesh1y,E_sorted1(:,:,12),'EdgeColor','none')
surf(klist_frac_mesh2x,klist_frac_mesh2y,E_sorted2(:,:,9),'EdgeColor','none')
surf(klist_frac_mesh2x,klist_frac_mesh2y,E_sorted2(:,:,10),'EdgeColor','none')
surf(klist_frac_mesh2x,klist_frac_mesh2y,E_sorted2(:,:,11),'EdgeColor','none')
surf(klist_frac_mesh2x,klist_frac_mesh2y,E_sorted2(:,:,12),'EdgeColor','none')
shading interp

%%
nk=Ndense+1;
totalk=nk^2*2;
klist_frac_mesh1x=reshape(klist_frac(1:totalk/2,1),nk,nk);
klist_frac_mesh1y=reshape(klist_frac(1:totalk/2,2),nk,nk);
klist_frac_mesh2x=reshape(klist_frac(totalk/2+1:end,1),nk,nk);
klist_frac_mesh2y=reshape(klist_frac(totalk/2+1:end,2),nk,nk);

Enk1=reshape(out.E_k(1:totalk/2,:),nk,nk,[]);
[E_sorted1, idx] = sort(Enk1, 3);   % idx: (Nk, Nk, Nband)
Enk2=reshape(out.E_k(totalk/2+1:end,:),nk,nk,[]);
[E_sorted2, idx] = sort(Enk2, 3);   % idx: (Nk, Nk, Nband)
figure()
hold on
plot(klist_frac(:,1),klist_frac(:,2),'ro')
plot(klist_frac(1:totalk/2,1),klist_frac(1:totalk/2,2),'bo')
plot(klist_frac_mesh1x(1,:),klist_frac_mesh1y(1,:));
plot(klist_frac_mesh2x(1,:),klist_frac_mesh2y(1,:));
figure()
hold on;
for j=1:nk
for i =1:20
   plot(E_sorted1(j,:,i))
   plot(E_sorted2(j,:,i))
end
end

%%

b1 = g.b(1:2,1); b2 = g.b(1:2,2);
Ncoarse = 1; Ndense = 20; R_dense = 0.01;
% [k_cart2D, k_frac2D, w_k, region_id] = make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense);
[k_cart2D, k_frac2D, w_k, region_id] = make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense);
klist_frac = [k_frac2D, zeros(size(k_frac2D,1),1)]; % 2D->3D compatible
%%
klist_cart = klist_frac*g.b;

figure()
hold on;
region1=region_id==1;
region2=region_id==2;
region3=region_id==3;
plot(klist_cart(region1,1),klist_cart(region1,2),'ko');
plot(klist_cart(region2,1),klist_cart(region2,2),'bo');
plot(klist_cart(region3,1),klist_cart(region3,2),'ro');
% plot(klist_frac(region1,1),klist_frac(region1,2),'ko');
% plot(klist_frac(region2,1),klist_frac(region2,2),'bo');
% plot(klist_frac(region3,1),klist_frac(region3,2),'ro');
%%
klist_cart_1=klist_cart(region1,:);
klist_cart_2=klist_cart(region2,:);
klist_cart_3=klist_cart(region3,:);

klist_cart_2=reshape(klist_cart_2,21,21,[]);
% plot(klist_cart_2(1,:),)
plot(klist_cart_2(1,:,1),klist_cart_2(1,:,2))

%%
function out = hf_run_wannier_patch(g, opts, klist_frac, w_k)
% HF self-consistent calculation on a NON-UNIFORM k-list (patch mesh).
% Fock: direct definition (NO FFT), periodic wrap in fractional coords.

% ---------- defaults ----------
if ~isfield(opts,'nspin'),           opts.nspin = 2;    end
if ~isfield(opts,'kT'),              opts.kT    = 1e-4; end
if ~isfield(opts,'max_iter'),        opts.max_iter = 100; end
if ~isfield(opts,'tol'),             opts.tol   = 1e-6; end
if ~isfield(opts,'mixing'),          opts.mixing= 0.6;  end
if ~isfield(opts,'vmodel'),          opts.vmodel= '2D'; end
if ~isfield(opts,'eps'),             opts.eps   = 4.0;  end
if ~isfield(opts,'e2'),              opts.e2    = 14.3996; end % eV·Å
if ~isfield(opts,'r0'),              opts.r0    = 10.0; end
if ~isfield(opts,'hartree_enable'),  opts.hartree_enable = true; end
if ~isfield(opts,'fock_drop_q0'),    opts.fock_drop_q0  = true; end
if ~isfield(opts,'qeps'),            opts.qeps  = 1e-6; end
if ~isfield(opts,'verbose'),         opts.verbose = true; end

% 2D input -> 3D compatible
if size(klist_frac,2)==2
    klist_frac = [klist_frac, zeros(size(klist_frac,1),1)];
end
% wrap to [-1/2,1/2)
klist_frac = wrap_frac_pm(klist_frac);

Nk = size(klist_frac,1);
m  = size(g.ham,1);

% weights normalization
if ~isfield(opts,'Wnorm') || isempty(opts.Wnorm)
    opts.Wnorm = sum(w_k);
end
Wnorm = opts.Wnorm;

% orbital positions: directly use Cartesian Å (per your data)
tau = g.wpos;
assert(size(tau,2)==3, 'g.wpos must be (m x 3) Cartesian Å');

% ---------- H0(k) ----------
H0_k = zeros(Nk,m,m);
for i = 1:Nk
    k_cart = (g.b * klist_frac(i,:).').';
    H0_k(i,:,:) = H0k_from_wannier(g, k_cart);
end

% ---------- initial diag & rho ----------
E_k  = zeros(Nk,m);
U_k  = cell(Nk,1);
rho_k= zeros(Nk,m,m);
for i = 1:Nk
    Hk = squeeze(H0_k(i,:,:));
    [V,D] = eig((Hk+Hk')/2);
    E_k(i,:) = real(diag(D));
    U_k{i}   = V;
end

% electrons per cell (including spin)
if ~isfield(opts,'Efilling') || isempty(opts.Efilling)
    opts.Efilling = m * opts.nspin / 2; % default: half-filling
end
Efilling = opts.Efilling;

mu = find_mu_patch(E_k, w_k, Wnorm, Efilling, opts.kT, opts.nspin);
for i = 1:Nk
    Ek = squeeze(E_k(i,:)).';
    fk = fermi_dirac(Ek, mu, opts.kT);
    Uk = U_k{i};
    rho_spin = Uk * diag(fk) * Uk';
    rho_k(i,:,:) = opts.nspin * rho_spin; % total
end

% ---------- SCF ----------
history = zeros(opts.max_iter,1);
logmat  = nan(opts.max_iter,5);
mu_prev = mu; Etot_prev = NaN;

for it = 1:opts.max_iter
    % Fock: direct sum on nonuniform mesh (periodized via wrap)
    SigmaF_k = fock_direct_patch(g, klist_frac, w_k, Wnorm, rho_k, opts);

    % Hartree: uniform (k-independent), orbital-resolved, diagonal
    if opts.hartree_enable
        SigmaH_k = hartree_uniform_smallq_patch(g, tau, rho_k, w_k, Wnorm, opts);
    else
        SigmaH_k = zeros(Nk,m,m);
    end

    % Mean-field Hamiltonian
    HMF_k = H0_k + SigmaH_k + SigmaF_k;

    % Diagonalize and update rho
    E_k_new = zeros(Nk,m);
    U_k_new = cell(Nk,1);
    parfor i = 1:Nk
        Hk = squeeze(HMF_k(i,:,:));
        [V,D] = eig((Hk+Hk')/2);
        E_k_new(i,:) = real(diag(D));
        U_k_new{i}   = V;
    end
    mu = find_mu_patch(E_k_new, w_k, Wnorm, Efilling, opts.kT, opts.nspin);

    rho_new = zeros(Nk,m,m);
    parfor i = 1:Nk
        Ek = squeeze(E_k_new(i,:)).';
        fk = fermi_dirac(Ek, mu, opts.kT);
        Uk = U_k_new{i};
        rho_spin = Uk * diag(fk) * Uk';
        rho_new(i,:,:) = opts.nspin * rho_spin;
    end

    % diagnostics (before mixing)
    drho   = rho_new - rho_k;
    nr     = norm(drho(:)) / max(1, norm(rho_k(:)));
    nr_inf = max(abs(drho(:)));
    Etot   = hf_total_energy_patch(H0_k, SigmaH_k, SigmaF_k, rho_new, w_k, Wnorm);

    % if opts.verbose
    %     if isnan(Etot_prev)
    %         fprintf('SCF %3d | nr=%.3e | max|dρ|=%.3e | Δμ=% .3e eV | E=% .10f eV\n',...
    %             it, nr, nr_inf, (mu-mu_prev), Etot);
    %     else
    %         fprintf('SCF %3d | nr=%.3e | max|dρ|=%.3e | Δμ=% .3e eV | E=% .10f eV (ΔE=% .3e)\n',...
    %             it, nr, nr_inf, (mu-mu_prev), Etot, (Etot-Etot_prev));
    %     end
    % end
    fprintf('SCF %3d | nr=%.3e | max|dρ|=%.3e | Δμ=% .3e eV | E=% .10f eV\n',...
                it, nr, nr_inf, (mu-mu_prev), Etot);

    logmat(it,:) = [it, nr, nr_inf, mu, Etot];
    mu_prev   = mu;
    Etot_prev = Etot;

    % mixing & stop
    rho_k = (1-opts.mixing)*rho_k + opts.mixing*rho_new;
    E_k   = E_k_new; U_k = U_k_new;

    history(it) = nr;
    if nr < opts.tol
        history = history(1:it);
        logmat  = logmat(1:it,:);
        if opts.verbose
            fprintf('HF converged at iter %d: nr=%.3e\n', it, nr);
        end
        break;
    end
    if it==opts.max_iter && opts.verbose
        fprintf('HF stopped at max_iter %d: nr=%.3e\n', it, nr);
        logmat  = logmat(1:it,:);
    end
end

% output
out.rho_k    = rho_k;
out.E_k      = E_k;
out.U_k      = U_k;
out.SigmaF_k = SigmaF_k;
out.SigmaH_k = SigmaH_k;
out.H0_k     = H0_k;
out.HMF_k    = HMF_k;
out.mu       = mu;
out.history  = history;
out.log      = logmat;
out.klist_frac = klist_frac;
out.Wnorm      = Wnorm;
end


function SigmaF_k = fock_direct_patch(g, klist_frac, w_k, Wnorm, rho_k, opts)
% Parfor over k_i; vectorized q-computation per i.
Nk = size(klist_frac,1);
m  = size(rho_k,2);

tau = g.wpos;
dZ  = abs(tau(:,3) - tau(:,3).');

SigmaF_k = zeros(Nk,m,m);

parfor i = 1:Nk
    ki = klist_frac(i,:);

    % vectorized q for all j
    qf_all = wrap_frac_pm( ki - klist_frac );    % Nk x 3
    qv_all = (g.b * qf_all.').';                 % Nk x 3
    qmag    = sqrt(sum(qv_all.^2,2));            % Nk x 1
    drop    = (opts.fock_drop_q0) & (qmag < opts.qeps);

    S = zeros(m,m);
    for j = 1:Nk
        if drop(j), continue; end

        % Vscalar
        qm = max(qmag(j), opts.qeps);
        if qm==0, continue; end
        switch lower(opts.vmodel)
          case '2d'
            Vscalar = 2*pi*opts.e2/(opts.eps*qm);
          case 'keldysh'
            Vscalar = opts.e2 /( qm*(1+opts.r0*qm) );
          otherwise
            error('Unknown vmodel');
        end
        if Vscalar==0, continue; end

        M = exp(-qmag(j) * dZ);
        if isfield(opts,'use_tau_phase') && opts.use_tau_phase
            dq = qv_all(j,:);
            phase = exp(-1i * (dq(1)*(tau(:,1)-tau(:,1).') + ...
                               dq(2)*(tau(:,2)-tau(:,2).') + ...
                               dq(3)*(tau(:,3)-tau(:,3).')));
        else
            phase = 1.0;
        end
        Vab = Vscalar * M .* phase;

        Rj  = squeeze(rho_k(j,:,:));
        S   = S - (w_k(j)/Wnorm) * ( Vab .* Rj );
    end
    SigmaF_k(i,:,:) = 0.5*(S + S');
end
end


function SigmaH_k = hartree_uniform_smallq_patch(g, tau, rho_k, w_k, Wnorm, opts)
% k-independent, diagonal Hartree shift; layer-resolved; remove common mode.

[Nk,m,~] = size(rho_k);

% average occupation per orbital (per cell, total including spin)
nbar = zeros(m,1);
for a=1:m
    nbar(a) = (1/Wnorm) * sum( w_k(:) .* real( reshape(rho_k(:,a,a),[],1) ) );
end

% reference/background
if isfield(opts,'nref') && numel(opts.nref)==m
    nref = opts.nref(:);
else
    nref = ones(m,1) * (sum(nbar)/m);
end
dn = nbar - nref;

% approximate V_ab(q->0) by sampling a few tiny q along b1,b2,b3
qdirs = eye(3);
Vsmall = zeros(m,m,3);
z = tau(:,3);
for iq=1:3
    qf = 1e-4 * qdirs(iq,:);   % tiny step in fractional coords
    qv = (g.b * qf.').';  qmag = norm(qv); if qmag < 1e-16, continue; end
    switch lower(opts.vmodel)
        case '2d'
            Vscalar = 2*pi*opts.e2/(opts.eps*qmag);
        case 'keldysh'
            Vscalar = opts.e2 /( qmag*(1+opts.r0*qmag) );
        otherwise
            error('Unknown vmodel');
    end
    Vsmall(:,:,iq) = Vscalar * exp(-qmag*abs(z(:)-z(:).'));
end
Vab0 = mean(Vsmall,3,'omitnan');

SigH_diag = Vab0 * dn;
SigH_diag = SigH_diag - mean(SigH_diag);  % remove common mode

SigmaH_k = zeros(Nk,m,m);
for a=1:m
    SigmaH_k(:,a,a) = SigH_diag(a);
end
end


function mu = find_mu_patch(E_k, w_k, Wnorm, Nelec, kT, nspin)
emin = min(E_k(:)) - 5*max(kT,1e-4);
emax = max(E_k(:)) + 5*max(kT,1e-4);
for it=1:80
    mu = 0.5*(emin+emax);
    f  = fermi_dirac(E_k, mu, kT);
    Ne = nspin * (1/Wnorm) * sum( w_k(:) .* sum(f,2) );
    if Ne > Nelec, emax = mu; else, emin = mu; end
    if abs(Ne-Nelec) < 1e-12*max(1,Nelec), break; end
end
end

function f = fermi_dirac(E, mu, kT)
if kT<=0
    f = double(E < mu);
else
    x = (E - mu)/kT;
    x = max(min(x, 40), -40);
    f = 1 ./ (1 + exp(x));
end
end

function Etot = hf_total_energy_patch(H0_k, SigmaH_k, SigmaF_k, rho_k, w_k, Wnorm)
Nk = size(rho_k,1);
Et = 0.0;
for i=1:Nk
    H0  = squeeze(H0_k(i,:,:));  H0=(H0+H0')/2;
    SH  = squeeze(SigmaH_k(i,:,:)); SH=(SH+SH')/2;
    SF  = squeeze(SigmaF_k(i,:,:)); SF=(SF+SF')/2;
    rho = squeeze(rho_k(i,:,:));    rho=(rho+rho')/2;
    Et = Et + w_k(i) * real(trace( H0*rho + 0.5*(SH+SF)*rho ));
end
Etot = Et / Wnorm;
end

function Hk = H0k_from_wannier(g, kvec)
m = size(g.ham,1);
Hk = zeros(m,m);
L  = size(g.hopr,1);
for l = 1:L
    Rfrac = g.hopr(l,:);
    Rvec  = (g.a * Rfrac.').';   % Å  (g.hopr is fractional)
    phase = exp(1i * dot(kvec, Rvec));
    Hk = Hk + phase * g.ham(:,:,l);
end
Hk = (Hk + Hk')/2;
end

function xw = wrap_frac_pm(x)
% wrap each component to [-1/2, 1/2)
xw = x - round(x);
end



function [k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense)
%MAKE_KMESH_HEX_PATCH_01  k-mesh on hexagonal BZ (frac coords in [0,1))
%   Coarse mesh on full BZ + dense patches around K1=(1/3,2/3), K2=(2/3,1/3).
%
% INPUT:
%   b1, b2   : [2x1] reciprocal lattice vectors, k = k1*b1 + k2*b2
%   Ncoarse  : coarse grid density per direction (fractional coords in [0,1))
%   Ndense   : dense patch density per direction
%   R_dense  : half-width of dense patch (in fractional coords)
%
% OUTPUT:
%   k_cart   : [Nk x 2] cartesian k points (kx, ky)
%   k_frac   : [Nk x 2] fractional coords (k1, k2) in [0,1)
%   weight   : [Nk x 1] integration weights, sum(weight) = area_BZ
%   region_id: [Nk x 1], 1 = coarse region, 2 = K1 patch, 3 = K2 patch

    % ---- 0. 输入检查 ----
    if nargin < 5
        error('Usage: make_kmesh_hex_patch_01(b1, b2, Ncoarse, Ndense, R_dense)');
    end

    b1 = b1(:);  b2 = b2(:);
    if numel(b1) ~= 2 || numel(b2) ~= 2
        error('b1 and b2 must be 2x1 column vectors.');
    end

    B = [b1, b2];                 % 2x2
    area_BZ = abs(det(B));        % BZ 平行四边形面积

    % ---- 1. K1 / K2 的 frac 坐标（[0,1)）----
    % 按你说的约定：Gamma=(0,0), K1=(1/3,2/3), K2=(2/3,1/3)
    K1_frac = [1/3, 2/3];
    K2_frac = [2/3, 1/3];

    % ---- 2. whole-BZ coarse grid: frac in [0,1) x [0,1) ----
    d_coarse = 1 / Ncoarse;
    k1c = 0 : d_coarse : 1 - d_coarse;
    k2c = 0 : d_coarse : 1 - d_coarse;
    [K1c, K2c] = meshgrid(k1c, k2c);
    kfrac_coarse = [K1c(:), K2c(:)];   % Nc x 2, in [0,1)

    % ---- 3. 在 K1 / K2 周围构造 dense patch (fractional) ----
    d_dense = (2 * R_dense) / Ndense;
    d1 = -R_dense : d_dense : R_dense;
    d2 = -R_dense : d_dense : R_dense;
    [D1, D2] = meshgrid(d1, d2);

    kfrac_dense_K1 = [D1(:) + K1_frac(1), D2(:) + K1_frac(2)];
    kfrac_dense_K2 = [D1(:) + K2_frac(1), D2(:) + K2_frac(2)];

    % 折回 [0,1)（mod 1）
    kfrac_dense_K1 = kfrac_dense_K1 - floor(kfrac_dense_K1);
    kfrac_dense_K2 = kfrac_dense_K2 - floor(kfrac_dense_K2);

    % ---- 4. 从 coarse 网格中移除落在两个 dense patch 里的点 ----
    % patch 条件（直接用 frac 空间的矩形窗口）
    deltaK1 = abs(kfrac_coarse - K1_frac);
    inK1    = (deltaK1(:,1) <= R_dense) & (deltaK1(:,2) <= R_dense);

    deltaK2 = abs(kfrac_coarse - K2_frac);
    inK2    = (deltaK2(:,1) <= R_dense) & (deltaK2(:,2) <= R_dense);

    in_patch = inK1 | inK2;
    kfrac_coarse_out = kfrac_coarse(~in_patch, :);

    % ---- 5. 合并 & 去重 ----
    kfrac_all = [kfrac_coarse_out;
                 kfrac_dense_K1;
                 kfrac_dense_K2];

    key = round(kfrac_all * 1e8);
    [~, ia] = unique(key, 'rows', 'stable');
    k_frac = kfrac_all(ia, :);          % 所有 frac 点 (k1,k2) in [0,1)

    Nk = size(k_frac, 1);

    % ---- 6. region_id: 1 coarse, 2 K1, 3 K2 ----
    region_id = ones(Nk, 1);

    deltaK1_all = abs(k_frac - K1_frac);
    inK1_all    = (deltaK1_all(:,1) <= R_dense) & (deltaK1_all(:,2) <= R_dense);

    deltaK2_all = abs(k_frac - K2_frac);
    inK2_all    = (deltaK2_all(:,1) <= R_dense) & (deltaK2_all(:,2) <= R_dense);

    region_id(inK1_all) = 2;
    region_id(inK2_all) = 3;

    % ---- 7. 积分权重（先粗略按区域分配，再归一化）----
    area_patch_frac  = (2 * R_dense)^2;   % 每个 patch 的 frac 面积
    area_coarse_frac = 1 - 2 * area_patch_frac;

    idx_coarse = (region_id == 1);
    idx_K1     = (region_id == 2);
    idx_K2     = (region_id == 3);

    N_coarse = nnz(idx_coarse);
    N_K1     = nnz(idx_K1);
    N_K2     = nnz(idx_K2);

    weight_frac = zeros(Nk, 1);

    if N_coarse > 0 && area_coarse_frac > 0
        weight_frac(idx_coarse) = area_coarse_frac / N_coarse;
    else
        weight_frac(idx_coarse) = 1;   % 极端情况占位
    end
    if N_K1 > 0
        weight_frac(idx_K1) = area_patch_frac / N_K1;
    end
    if N_K2 > 0
        weight_frac(idx_K2) = area_patch_frac / N_K2;
    end

    % 全局归一化：sum(weight) = area_BZ
    S_frac = sum(weight_frac);
    if S_frac == 0
        weight_frac(:) = 1 / Nk;
        S_frac = 1;
    end
    weight_frac = weight_frac / S_frac;
    weight = area_BZ * weight_frac;

    % ---- 8. 变成真实 k 向量 ----
    k_cart = (B * k_frac.').';     % Nk x 2, (kx,ky)

end

function [k_cart, k_frac, weight, region_id] = make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense)
%MAKE_KMESH_HEX_PATCH  Build a k-mesh on a hexagonal lattice BZ
%   with coarse sampling everywhere and dense patches around K and K'.
%
% INPUT:
%   b1, b2   : [2x1] real-space column vectors of reciprocal lattice
%              (units: 1/length). They define k = k1*b1 + k2*b2.
%   Ncoarse  : scalar, coarse grid points per direction in fractional
%              coordinates (k1,k2) in [-0.5,0.5).
%   Ndense   : scalar, dense patch points per direction around each valley.
%   R_dense  : half-width of the dense patch (in fractional coordinates).
%              Patch in frac-space: |k1 - K1| <= R_dense, |k2 - K2| <= R_dense.
%
% OUTPUT:
%   k_cart   : [Nk x 2] array, cartesian (kx, ky) of all k-points.
%   k_frac   : [Nk x 2] array, fractional (k1, k2) such that k = k1*b1 + k2*b2.
%   weight   : [Nk x 1] integration weights (in k-space area units).
%              Sum_j weight(j) = area of the BZ parallelogram = |det([b1 b2])|.
%   region_id: [Nk x 1] integer flag:
%              1 = coarse region
%              2 = dense patch around K
%              3 = dense patch around K'
%
% NOTES:
%   - The BZ here is taken as the parallelogram in fractional coords
%     (k1, k2) in [-0.5, 0.5) x [-0.5, 0.5).
%   - K, K' are placed at (k1,k2) = (2/3, 1/3) and (1/3, 2/3),
%     and then folded into (-0.5, 0.5] via frac = frac - round(frac).
%   - R_dense should be small enough that the two patches do not overlap.
%
% EXAMPLE:
%   % Suppose real-space primitive vectors are:
%   a = 1.0;
%   a1 = a * [1; 0];
%   a2 = a * [0.5; sqrt(3)/2];
%   % Reciprocal lattice:
%   A  = [a1, a2];
%   B  = 2*pi * inv(A).';   % columns b1, b2
%   b1 = B(:,1); b2 = B(:,2);
%
%   [k_cart, k_frac, w, region_id] = make_kmesh_hex_patch(b1, b2, ...
%                                     20, 40, 0.08);
%
%   % k_cart, w 就可以直接用在 HFMF 的 k 积分里:
%   % sum_j w(j) * f(k_cart(j,:))

    % -------- 0. 基本检查 --------
    if nargin < 5
        error('Usage: make_kmesh_hex_patch(b1, b2, Ncoarse, Ndense, R_dense)');
    end

    b1 = b1(:);
    b2 = b2(:);
    if numel(b1) ~= 2 || numel(b2) ~= 2
        error('b1 and b2 must be 2x1 column vectors.');
    end

    % Reciprocal lattice matrix and its determinant (BZ area)
    B = [b1, b2];                 % 2x2
    area_BZ = abs(det(B));        % area of parallelogram BZ in k-space

    % -------- 1. 定义 K / K' 的 fractional 坐标并折回 --------
    % 原始 frac 坐标
    K_frac  = [2/3, 1/3];
    Kp_frac = [1/3, 2/3];

    % 折回到 (-0.5, 0.5] 区间（mod 1）
    K_frac  = K_frac  - round(K_frac);
    Kp_frac = Kp_frac - round(Kp_frac);

    % -------- 2. 构造 whole-BZ coarse grid (fractional coords) --------
    d_coarse = 1 / Ncoarse;  % step in frac space
    k1c = -0.5 : d_coarse : 0.5 - d_coarse;  % [-0.5, 0.5)
    k2c = -0.5 : d_coarse : 0.5 - d_coarse;
    [K1c, K2c] = meshgrid(k1c, k2c);
    kfrac_coarse = [K1c(:), K2c(:)];   % Nc x 2

    % -------- 3. 在 K / K' 周围构造 dense patch (fractional coords) --------
    % patch 是中心在 K / K'，边长 2*R_dense 的方形区域
    d_dense = (2*R_dense) / Ndense;
    d1 = -R_dense : d_dense : R_dense;
    d2 = -R_dense : d_dense : R_dense;
    [D1, D2] = meshgrid(d1, d2);

    kfrac_dense_K  = [D1(:) + K_frac(1),  D2(:) + K_frac(2)];
    kfrac_dense_Kp = [D1(:) + Kp_frac(1), D2(:) + Kp_frac(2)];

    % 再次折回到 (-0.5, 0.5] 以保证都在同一个 fundamental domain
    kfrac_dense_K  = kfrac_dense_K  - round(kfrac_dense_K);
    kfrac_dense_Kp = kfrac_dense_Kp - round(kfrac_dense_Kp);

    % -------- 4. 从 coarse 网格中移除落在两个 dense patch 里的点 --------
    % 使用 frac 空间定义 patch: |k1 - K1| <= R_dense 且 |k2 - K2| <= R_dense
    % 注意这里已经在 (-0.5,0.5] 里，因此不考虑额外的周期性
    deltaK  = abs(kfrac_coarse - K_frac);   % Nc x 2
    inK     = (deltaK(:,1) <= R_dense) & (deltaK(:,2) <= R_dense);

    deltaKp = abs(kfrac_coarse - Kp_frac);
    inKp    = (deltaKp(:,1) <= R_dense) & (deltaKp(:,2) <= R_dense);

    in_patch = inK | inKp;

    kfrac_coarse_out = kfrac_coarse(~in_patch, :);   % coarse outside patches

    % -------- 5. 合并 coarse_out + dense_K + dense_Kp 并去重 --------
    kfrac_all = [kfrac_coarse_out;
                 kfrac_dense_K;
                 kfrac_dense_Kp];

    % 为了 unique 稳定一点，先做缩放+round 再 unique
    key = round(kfrac_all * 1e8);
    [~, idx_unique] = unique(key, 'rows', 'stable');
    kfrac = kfrac_all(idx_unique, :);
    % 
    % % -------- 6. 标记 region_id: coarse(1), K_patch(2), Kp_patch(3) --------
    % % 我们再次用 frac 空间的窗口来检查
    % Nk = size(kfrac, 1);
    % region_id = zeros(Nk, 1);
    % 
    % % coarse region: 先设为 1
    % region_id(:) = 1;
    % 
    % % K patch
    % deltaK_all  = abs(kfrac - K_frac);
    % inK_all     = (deltaK_all(:,1) <= R_dense) & (deltaK_all(:,2) <= R_dense);
    % 
    % % K' patch
    % deltaKp_all = abs(kfrac - Kp_frac);
    % inKp_all    = (deltaKp_all(:,1) <= R_dense) & (deltaKp_all(:,2) <= R_dense);
    % 
    % region_id(inK_all)  = 2;
    % region_id(inKp_all) = 3;

    % -------- 6. 标记 region_id: coarse(1), K_patch(2), Kp_patch(3) --------
Nk = size(kfrac, 1);
region_id = ones(Nk, 1);   % 先全部设为 coarse=1

% 给边界留一点浮点余量
tol = 10 * eps(R_dense);   % 或者 tol = 1e-10; 根据你 R_dense 大小调整

% K patch
deltaK_all  = abs(kfrac - K_frac);
inK_all     = (deltaK_all(:,1) <= R_dense + tol) & ...
              (deltaK_all(:,2) <= R_dense + tol);

% K' patch
deltaKp_all = abs(kfrac - Kp_frac);
inKp_all    = (deltaKp_all(:,1) <= R_dense + tol) & ...
              (deltaKp_all(:,2) <= R_dense + tol);

region_id(inK_all)  = 2;
region_id(inKp_all) = 3;


    % -------- 7. 计算积分权重 --------
    % frac-space 中 BZ 面积是 1.0 (平行四边形 [0,1)x[0,1) 或 [-0.5,0.5)x[-0.5,0.5))
    % 我们有两类区域：
    %   - coarse outside patches
    %   - two dense patches around K / K'
    %
    % coarse 区域的面积（frac-space）:
    %   area_total_frac = 1
    %   area_patch_K_frac  = (2*R_dense)^2
    %   area_patch_Kp_frac = (2*R_dense)^2
    %   area_coarse_frac   = 1 - area_patch_K_frac - area_patch_Kp_frac

    area_patch_K_frac  = (2*R_dense)^2;
    area_patch_Kp_frac = (2*R_dense)^2;
    area_coarse_frac   = 1 - area_patch_K_frac - area_patch_Kp_frac;

    if area_coarse_frac < 0
        warning('R_dense too large: patches overlap or exceed BZ area in frac space.');
    end

    % 每个区域的点数
    idx_coarse = find(region_id == 1);
    idx_K      = find(region_id == 2);
    idx_Kp     = find(region_id == 3);

    N_coarse = numel(idx_coarse);
    N_K      = numel(idx_K);
    N_Kp     = numel(idx_Kp);

    weight_frac = zeros(Nk, 1);

    % 均匀分配每个区域的 frac-area
    if N_coarse > 0
        w_coarse = area_coarse_frac / N_coarse;
        weight_frac(idx_coarse) = w_coarse;
    end
    if N_K > 0
        w_K = area_patch_K_frac / N_K;
        weight_frac(idx_K) = w_K;
    end
    if N_Kp > 0
        w_Kp = area_patch_Kp_frac / N_Kp;
        weight_frac(idx_Kp) = w_Kp;
    end

    % 从 frac-space 面积转换到真实 k-space 面积（乘以 |det(B)|）
    weight = area_BZ * weight_frac;

    % -------- 8. 变成真实 k 向量 --------
    k_cart = (B * kfrac.').';   % Nk x 2, each row (kx, ky)
    k_frac=kfrac;

end