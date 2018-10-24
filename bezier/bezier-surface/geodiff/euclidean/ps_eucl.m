% FUNCTION PS_EUCL(S,T,X,Y): 
% 		Computes the scalar product on the Euclidean space.
% ------------------------------------------------------------
% 
% INPUT: 	S : A point on the tangent space in X.
% 			T : A point on the tangent space in Y.
% 			X : [not mandatory] A point.
% 			Y : [not mandatory] Another point.
%
% OUTPUT: 	P : The scalar product between S and Y.
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

function p = ps_eucl(s,t,x,y)
    if nargin==3; x = 0; y = 0; end
    assert(size(s,1) == size(t,1) && size(s,2) == size(t,2));
    
    p = zeros(length(s),1);
    n = length(s);
    
    if isa(s,'cell') && isa (t,'cell') && isa(x,'cell') && isa(y,'cell')  
        m = length(s{1});
        for i=1:n; p(i) = reshape(s{i},1,m)*reshape(t{i},m,1); end
    elseif isa(s,'double') && isa(t,'double') && isa(x,'double') && isa(y,'double')
        m = length(s(:,:,1));
        for i=1:n p(i) = reshape(s(:,:,1),1,m)*reshape(t(:,:,i),m,1); end;
    end
end
