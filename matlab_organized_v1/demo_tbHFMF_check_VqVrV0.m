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
% 1) ---- choose patch ----
K0_frac = [1/3, 2/3];
L_frac  = [0.04, 0.04];
Nk      = [31, 31];       % [Nkx, Nky]

mesh = tbHFMF.build_kmesh_patch(g, K0_frac, L_frac, Nk);
%%
% 2) ---- dtau ----
ph = tbHFMF.precompute_dtau_phases(mesh, g.wpos);
check_dtau_phases(mesh,ph)

%%
pars.d_gate      = 60.0;
pars.e2_over_eps = 14.3996/6.0;
pars.q_small     = 1e-5;
pars.include_q2  = false;
pars.z0_mode     = 'center';

V = tbHFMF.build_Vq_doublegate_layered(g, mesh, pars);
check_Vq_Vr_reconstruction(g, mesh, pars, V);
%%
% Vq=reshape(V.Vq,101,101,20,20);
% Vq_K=squeeze(Vq(51,51,:,:))

% 4) convolution equality (pick any orbital pair)
demo_check_fock_conv_singlepair(g, K0_frac, L_frac, Nk, pars, 1, 3);
demo_check_fock_conv_singlepair(g, K0_frac, L_frac, Nk, pars, 1, 2);
%%
function rep = check_dtau_phases(mesh, ph)
%CHECK_DTAU_PHASES  Sanity checks for exp(±i p·(tau_a-tau_b)) pages.

Norb = ph.Norb;
Nky  = mesh.Nky; Nkx = mesh.Nkx;

icx = floor(Nkx/2)+1;
icy = floor(Nky/2)+1;

% 1) center p=0 should give phase = 1 for all (a,b)
e1 = 0;
for p = 1:ph.P
  e1 = max(e1, abs(ph.phase_plus(icy,icx,p) - 1));
  e1 = max(e1, abs(ph.phase_minus(icy,icx,p) - 1));
end

% 2) phase_minus should be conj(phase_plus)
e2 = max(abs(ph.phase_minus(:) - conj(ph.phase_plus(:))));

% 3) phase_plus * phase_minus = 1
e3 = max(abs(ph.phase_plus(:).*ph.phase_minus(:) - 1));

% 4) ab<->ba relation: phase_plus(ab) = conj(phase_plus(ba))
e4 = 0;
for a = 1:Norb
  for b = 1:Norb
    pab = a + (b-1)*Norb;
    pba = b + (a-1)*Norb;
    tmp = max(abs(ph.phase_plus(:,:,pab) - conj(ph.phase_plus(:,:,pba))), [], 'all');
    e4 = max(e4, tmp);
  end
end

rep = struct();
rep.center_err = e1;
rep.conj_err   = e2;
rep.unit_err   = e3;
rep.swap_err   = e4;

fprintf('[dtau] center phase err      = %.3e\n', rep.center_err);
fprintf('[dtau] phase_minus-conj err  = %.3e\n', rep.conj_err);
fprintf('[dtau] phase*phase^-1 err    = %.3e\n', rep.unit_err);
fprintf('[dtau] ab<->ba conj err      = %.3e\n', rep.swap_err);

end

function rep = check_Vq_Vr_reconstruction(g, mesh, pars, V)
%CHECK_VQ_VR_RECONSTRUCTION
% 1) Reconstruct Vq from Vr: Vq_rec = from_fft(ifft2(Vr))
% 2) Compare Vq_rec(center) with V0(a,b)
% 3) Compare sampled points with analytic double-gate formula

Norb = size(g.wpos,1);
Nky  = mesh.Nky; Nkx = mesh.Nkx;
icx = floor(Nkx/2)+1;
icy = floor(Nky/2)+1;

% reconstruct Vq from Vr
Vq_fft = ifft2(V.Vr);                           % FFT order
Vq_rec = fftshift(fftshift(Vq_fft,1),2);        % centered order (from_fft)

% check center vs V0
e0 = 0;
for a = 1:Norb
  for b = 1:Norb
    p = a + (b-1)*Norb;
    e0 = max(e0, abs(real(Vq_rec(icy,icx,p)) - real(V.V0(a,b))));
  end
end

% sample compare with analytic formula away from center
rng(2);
nsamp = 20;
ixs = randi(Nkx,[nsamp,1]);
iys = randi(Nky,[nsamp,1]);

% z positions (same shift rule as builder)
z = g.wpos(:,3);
if isfield(pars,'z0_mode') && ischar(pars.z0_mode) && strcmpi(pars.z0_mode,'center')
  z = z - 0.5*(max(z)+min(z));
elseif isfield(pars,'z0_mode') && isnumeric(pars.z0_mode)
  z = z - pars.z0_mode;
end

d  = pars.d_gate;
e2 = pars.e2_over_eps;
q  = mesh.qabs;

e1 = 0;
for t = 1:nsamp
  iy = iys(t); ix = ixs(t);
  if (iy==icy && ix==icx), continue; end

  qq = q(iy,ix);
  den = sinh(2*qq*d);

  % test a random orbital pair
  a = randi(Norb); b = randi(Norb);
  zgt = max(z(a),z(b));
  zlt = min(z(a),z(b));
  A = d - zgt;
  B = d + zlt;

  num = 2*sinh(qq*A)*sinh(qq*B);
  Vana = (2*pi*e2) * (num/den) / qq;    % eV·Å^2

  pidx = a + (b-1)*Norb;
  Vrec = real(Vq_rec(iy,ix,pidx));

  e1 = max(e1, abs(Vrec - Vana));
end

rep = struct();
rep.center_V0_abs_err = e0;
rep.sample_abs_err    = e1;

fprintf('[Vq/Vr] |Vq_rec(center)-V0| max abs = %.3e (eV·Å^2)\n', rep.center_V0_abs_err);
fprintf('[Vq/Vr] sample |Vq_rec - V_analytic| max abs = %.3e (eV·Å^2)\n', rep.sample_abs_err);

end



function demo_check_fock_conv_singlepair(g, K0_frac, L_frac, Nk, pars, a, b)
%DEMO_CHECK_FOCK_CONV_SINGLEPAIR  Verify FFT conv == direct cyclic sum for one (a,b).

mesh = tbHFMF.build_kmesh_patch(g, K0_frac, L_frac, Nk);

% build V and phases
V  = tbHFMF.build_Vq_doublegate_layered(g, mesh, pars);
ph = tbHFMF.precompute_dtau_phases(mesh, g.wpos);

Nky = mesh.Nky; Nkx = mesh.Nkx;
icx = floor(Nkx/2)+1;
icy = floor(Nky/2)+1;

% Reconstruct Vq pages from Vr
Vq_fft = ifft2(V.Vr);
Vq_rec = fftshift(fftshift(Vq_fft,1),2);     % centered
pidx = a + (b-1)*size(g.wpos,1);
Vab  = real(Vq_rec(:,:,pidx));               % centered V(q) for this pair

% random rho_ab(k) on centered grid
rng(3);
rho = randn(Nky,Nkx) + 1i*randn(Nky,Nkx);

% ---- FFT method (same math as our fock_fft) for this pair only ----
to_fft   = @(A) ifftshift(ifftshift(A,1),2);
from_fft = @(A) fftshift(fftshift(A,1),2);

rho_tilde = rho .* ph.phase_plus(:,:,pidx);

Sigma_tilde = from_fft( ifft2( fft2(to_fft(Vab)) .* fft2(to_fft(rho_tilde)) ) );
Sigma_fft = - Sigma_tilde .* ph.phase_minus(:,:,pidx);   % include vertex and minus sign

% ---- Direct cyclic sum (TWIST/UNTWIST version, matches FFT trick) ----
Sigma_dir = zeros(Nky,Nkx);
cX = icx; cY = icy;

for iky = 1:Nky
  for ikx = 1:Nkx
    s = 0;
    for iqy = 1:Nky
      for iqx = 1:Nkx
        % k-q with CENTERED kernel indexing: q corresponds to (iq-c)
        jky = mod( (iky - iqy + cY - 1), Nky ) + 1;
        jkx = mod( (ikx - iqx + cX - 1), Nkx ) + 1;

        % twist rho at (k-q)
        s = s + Vab(iqy,iqx) * ( rho(jky,jkx) * ph.phase_plus(jky,jkx,pidx) );
      end
    end
    % untwist at k
    Sigma_dir(iky,ikx) = - ph.phase_minus(iky,ikx,pidx) * s;
  end
end


err = max(abs(Sigma_fft(:) - Sigma_dir(:)));
fprintf('[Fock conv check] max|Sigma_fft - Sigma_dir| = %.3e\n', err);

end

