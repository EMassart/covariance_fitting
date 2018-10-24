% FUNCTION EXP_PSD_QUOTIENT(X,V,T): 
% 		 Computes the exponential map from a point X on the
% 		 manifold of fixed rank PSD matrices, with initial velocity V and at time T.
% ------------------------------------------------------------ 
% INPUT: 	X : A point on the manifold of fixed rank PSD matrices (given as a Y
%              factor, and not the full PSD matrix)
% 			V : A starting velocity.
% 			T : A time instant.
%
% OUTPUT: 	Y : The point reached on the Euclidean space
% ------------------------------------------------------------
% Author: Estelle Massart
%
% This file is a extension of the files of the project "C1 bezier paths on surfaces"
% by Gousenbourger et al to the manifold of PSD matrices.
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt
%
% ------------------------------------------------------------
% Last modification: October 24, 2018
% ------------------------------------------------------------
function y = exp_psd_quotient(x,v,t)

    if nargin < 3
        t = 1;
    end
    
    y = x;
    if isa(x,'cell') && isa(v,'cell')
        for i=1:length(x)
            y{i} = x{i} + v{i}.*t;
        end
    elseif isa(x,'double') && isa(v,'double')
        y = x + v.*t;
    else
        error('x and v must be same type');
    end
end
