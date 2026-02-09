function out = solve_band_path_E(g, klist_frac)
%SOLVE_BAND_PATH_E  Compute eigenvalues along a path (PARFOR, no if).
%
% Inputs:
%   g          : must have g.Rcart (use prep_tb_cache first)
%   klist_frac : Nk x 3 reduced coords (2D -> third column zeros)
%
% Output:
%   out.E      : Norb x Nk
%   out.kdist  : Nk x 1 cumulative distance (Å^-1)

Nk   = size(klist_frac,1);
Norb = size(g.ham,1);

% kcart + kdist
kcart = zeros(Nk,3);
for ik = 1:Nk
  kf = klist_frac(ik,:);
  kcart(ik,:) = kf(1)*g.b(1,:) + kf(2)*g.b(2,:) + kf(3)*g.b(3,:);
end
kdist = zeros(Nk,1);
for ik = 2:Nk
  kdist(ik) = kdist(ik-1) + norm(kcart(ik,:) - kcart(ik-1,:));
end

E = zeros(Norb,Nk);

parfor ik = 1:Nk
  Hk = tbHFMF.get_hk_atomgauge_single(g, klist_frac(ik,:));
  [ev,~] = tbHFMF.solve_hk(Hk);
  E(:,ik) = ev;
end

out.E = E;
out.kdist = kdist;
end
