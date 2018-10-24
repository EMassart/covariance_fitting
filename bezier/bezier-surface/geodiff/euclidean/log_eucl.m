% FUNCTION LOG_EUCL(X,Y): 
% 		Computes the logarithmic map from on the Euclidean space
% 		from X to Y.
% ------------------------------------------------------------
% 
% INPUT: 	X : A point
% 			Y : Another point
%
% OUTPUT: 	V : The initial velocity at T_x M
% ------------------------------------------------------------
% Authors: Pierre-Yves Gousenbourger
% ------------------------------------------------------------
% Versions
% 	19/03/2014: first version.
% 	18/06/2015: header changed.

% This file comes from the project "C1 bezier paths on surfaces"
% by Gousenbourger et al. 
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt

% ------------------------------------------------------------
function v = log_eucl(x,y)
    v = x;
    if isa(x,'cell') && isa(y,'cell')
        for i=1:length(x)
            v{i} = y{i} - x{i};
        end
    elseif isa(x,'double') && isa(y,'double')
        v = y - x;
    else
        error('x and y must be same type');
    end
end
