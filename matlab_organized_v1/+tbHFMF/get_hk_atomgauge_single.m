function Hk = get_hk_atomgauge_single(g, k_frac)
%GET_HK_ATOMGAUGE_SINGLE  Atom-gauge Bloch Hamiltonian at a single k.
%
% Convention:
%   H(k) = sum_R  t(R) * exp(+i k·R)
%
% Inputs:
%   g.Rcart : Nh x 3 (Å)   (must be precomputed by prep_tb_cache)
%   g.b     : 3x3 (rows are b1,b2,b3) (Å^-1)
%   g.ham   : Norb x Norb x Nh (eV)
%   k_frac  : 1x3 reduced coords [k1 k2 k3]
%
% Output:
%   Hk      : Norb x Norb, Hermitian enforced

% reduced -> cart (Å^-1)
k_cart = k_frac(1)*g.b(1,:) + k_frac(2)*g.b(2,:) + k_frac(3)*g.b(3,:);  % 1x3

% phases for each hop
phi = exp( 1i * (g.Rcart * k_cart(:)) );  % Nh x 1

% accumulate
Norb = size(g.ham,1);
Nh   = size(g.ham,3);
Hk = zeros(Norb,Norb,'like',1+1i);
for l = 1:Nh
  Hk = Hk + g.ham(:,:,l) * phi(l);
end

% enforce Hermitian
Hk = (Hk + Hk')/2;

end
