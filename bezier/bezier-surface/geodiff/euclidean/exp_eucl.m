% FUNCTION EXP_EUCL(X,V,T): 
% 		 Computes the exponential map from a point X on the
% 		 Euclidean space, with initial velocity V and at time T.
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces"
% 
% INPUT: 	X : A point on the Euclidean space of size n
% 			V : A starting velocity.
% 			T : A time instant.
%
% OUTPUT: 	Y : The point reached on the Euclidean space
% ------------------------------------------------------------
% Author: Pierre-Yves Gousenbourger
% ------------------------------------------------------------
% Versions
% 	19/03/2014: first version.
% 	18/06/2015: header changed.
% ------------------------------------------------------------
function y = exp_eucl(x,v,t)
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
