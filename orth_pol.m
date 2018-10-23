function [Q,H] = orth_pol(A)
% Computes the orthogonal factor of the polar decomposition A = HQ
% Author : E. Massart
% Version : March 8, 2018


[U,S,V] = svd(A);
H = U*S*U';
Q = U*V';

end