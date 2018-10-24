% FUNCTION TP_EUCL(X,Y,V): 
% 		Parallel transport of a vector V from X to Y on the
%		Euclidean space.
% ------------------------------------------------------------
% 
% INPUT: 	X : A point.
% 			Y : Another point.
% 			V : Vector on the tangent space of X to transport.
%
% OUTPUT: 	S : The transported vector.
% ------------------------------------------------------------
% Author: Pierre-Yves Gousenbourger
% ------------------------------------------------------------
% Versions
% 	19/03/2014: first version.
% 	18/06/2015: header changed.

% This file comes from the project "C1 bezier paths on surfaces"
% by Gousenbourger et al. 
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt


% ------------------------------------------------------------

function s = tp_eucl(x,y,v)
    s = v;
end
