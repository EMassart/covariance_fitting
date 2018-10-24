% FUNCTION NORM_EUCL(X): 
% 		 Norms  a point X on the Euclidean space.
% ------------------------------------------------------------
% 
% INPUT: 	X : A point on the Euclidean space of size n
%
% OUTPUT: 	Y : X normed.
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
function y = norm_eucl(x)
    y = x;
    if isa(x,'cell')
        for i=1:length(x); y{i} = x{i}./norm(x{i}); end
    elseif isa(x,'double')
        for i=1:size(x,1); y(i,:) = x(i,:)./norm(x(i,:)); end
    end
end
