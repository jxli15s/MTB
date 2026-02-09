function g = prep_tb_cache(g)
%PREP_TB_CACHE  Precompute cached quantities for fast H(k).
% Adds:
%   g.Rcart : Nh x 3, Cartesian hop vectors (Å)

g.Rcart = g.hopr * g.a;   % Nh x 3
end
