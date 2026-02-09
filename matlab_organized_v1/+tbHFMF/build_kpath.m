function path = build_kpath(g, k_nodes_frac, n_per_seg, labels)
%BUILD_KPATH  Build interpolated k-path with correct tick indices (nodes).
%
% Convention:
%   Each segment includes both endpoints.
%   To avoid duplicating internal nodes, we drop the FIRST point of segment s>1.
% Then tick_idx(m) always points to the exact node m in klist.

% enforce 3 columns (2D -> add kz=0)
if size(k_nodes_frac,2) == 2
  k_nodes_frac = [k_nodes_frac, zeros(size(k_nodes_frac,1),1)];
end

M = size(k_nodes_frac,1);
tick_idx = zeros(1,M);

klist = zeros(0,3);

for s = 1:(M-1)
  % node s will be the first point of this segment in the final list
  tick_idx(s) = size(klist,1) + 1;

  A = k_nodes_frac(s,:);
  B = k_nodes_frac(s+1,:);

  t  = linspace(0,1,n_per_seg).';
  ks = (1-t).*A + t.*B;          % n_per_seg x 3

  if s > 1
    ks = ks(2:end,:);            % drop duplicated node A
  end

  klist = [klist; ks]; %#ok<AGROW>
end

% last node M is the last point in the list
tick_idx(M) = size(klist,1);

Nk = size(klist,1);

% compute kcart and kdist (Å^-1)
kcart = zeros(Nk,3);
for ik = 1:Nk
  kf = klist(ik,:);
  kcart(ik,:) = kf(1)*g.b(1,:) + kf(2)*g.b(2,:) + kf(3)*g.b(3,:);
end

kdist = zeros(Nk,1);
for ik = 2:Nk
  kdist(ik) = kdist(ik-1) + norm(kcart(ik,:) - kcart(ik-1,:));
end

path = struct();
path.klist_frac = klist;
path.kcart      = kcart;
path.kdist      = kdist;
path.tick_idx   = tick_idx;

if nargin >= 4
  path.labels = labels;
else
  path.labels = {};
end

end
