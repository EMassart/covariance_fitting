% phat(p, pref, i, j): 
% 		 Compputes the k_th element of the P(p) vector.
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	P 	 : [cell] the interpolation points of the problem.
% 					(1 x n)
% 
% 			k 	 : [int] The requested component of phat.
%
% OUTPUT: 	PHAT : [matrix] The k_th element of P(p)
% ------------------------------------------------------------
% Authors:
% 	Pierre-Yves Gousenbourger
% 	Paul Striewski
%
% Versions
% 	04/11/2015: First version.
% ------------------------------------------------------------

function phat = phat(p,k)
	
	global_variables;
	
	% Useful data
	M = length(p); assert(M > 2, 'the length of p must be > 2');
	
	if M == 3
		phat = p{2} - p{1}/6 - p{3}/6;
	else
		if k == 2
			phat = 	p{2} - p{1}/6;
		elseif k == M-1
			phat =	p{k} - p{k+1}/6;
		else
			phat = 	p{k};
		end
	end
end
