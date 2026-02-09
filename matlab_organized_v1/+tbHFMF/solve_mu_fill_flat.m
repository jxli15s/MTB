function mu = solve_mu_fill_flat(Ek_flat, mesh, opts)
% Ek_flat: Nact x Nk
E  = real(Ek_flat);
kT = opts.kT;

Nact = size(E,1);
Nk   = size(E,2);

use_cell = isfield(opts,'n_target_cell') && ~isempty(opts.n_target_cell);
if ~use_cell
  error('set opts.n_target_cell');
end
target = opts.n_target_cell;

% n(mu) = area_frac * <sum_n f_n(k)>_k
n_of_mu = @(mu) mesh.area_frac * mean( sum( fermi_stable(E, mu, kT), 1 ) );

emin = min(E,[],'all');
emax = max(E,[],'all');
W    = 50*max(kT,1e-6) + 1;

muL = emin - W;
muR = emax + W;

nL = n_of_mu(muL);
nR = n_of_mu(muR);

% 如果 target 不可达，直接 clamp（否则你会觉得"mu 很怪"）
if target <= nL, mu = muL; return; end
if target >= nR, mu = muR; return; end

for it = 1:80
  muM = 0.5*(muL+muR);
  nM  = n_of_mu(muM);
  if nM > target
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
x = (E - mu) ./ kT;
x = max(min(x, 50), -50);
F = 1 ./ (1 + exp(x));
end
