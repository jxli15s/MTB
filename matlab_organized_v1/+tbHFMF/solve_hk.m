function [E,U] = solve_hk(Hk)
%SOLVE_HK  Diagonalize Hk with fixed conventions (no if).
%
% Always:
%   Hk <- (Hk+Hk')/2
%   [U,D] = eig(Hk)
%   E sorted ascending
%
% Output:
%   E : Norb x 1 (real)
%   U : Norb x Norb, columns are eigenvectors

Hk = (Hk + Hk')/2;
[U,D] = eig(Hk,'vector');
E = real(D);
[E,ord] = sort(E,'ascend');
U = U(:,ord);

end
