% compute_alpha(p, A, m, d): 
% 		 Computes the intermediate alpha of a C^2-Bezier spline in R^n,
% 		 such that 
% 				alpha = b_i^m = 'sum_{k=m-d}^{m+d} A_{k,m} p_k
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	p : [cell] a line cell of points on a manifold M.
% 				(1 x M)
%
% 			A : [matrix] a matrix defining the linear combination.
% 				(M-2 x M-2)
%
%			m : [int] the position in the row
%
% 			d : [int] the size of the window of coefficients taken in A.
% 				(optional)
%
% OUTPUT: 	alpha : [cell] intermediate alphas in the computation of b.
%
% ------------------------------------------------------------
% Authors:
% 	Pierre-Yves Gousenbourger
% 	Paul Striewski
%
% Versions
% 	04/11/2015: First version.

% This file comes from the project "C1 bezier paths on surfaces"
% by Gousenbourger et al. 
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt

% ------------------------------------------------------------

function alpha = compute_alpha(p, A, m, d)
	
	% Parameters
	M 	  = length(p);
	
	% Computation
	switch m
	case 1
		alpha = p{1};
	case M
		alpha = p{M};
	otherwise
		upperBound = min(M-1,m+d);
		lowerBound = max(2,m-d);
		alpha = zeros(size(p{m}));
		for k = lowerBound:upperBound
			alpha = alpha + A(k-1,m-1)*phat(p,k); 
		end
	end
end
