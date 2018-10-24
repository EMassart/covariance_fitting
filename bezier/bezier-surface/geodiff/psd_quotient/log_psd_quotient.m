% FUNCTION LOG_PSD_QUOTIENT(X,Y): 
% 		Computes the logarithmic map on the manifold of fixed rank PSD matrices 
%
% ------------------------------------------------------------
% 
% INPUT: 	X : A point on the manifold of fixed rank PSD matrices (given as a Y
%              factor, and not the full PSD matrix)
% 			Y : Another point on the manifold of fixed rank PSD matrices (given as a Y
%              factor, and not the full PSD matrix)
%
% OUTPUT: 	V : The initial velocity at T_x M
% ------------------------------------------------------------
% Author: Estelle Massart
%
% This file is a extension of the files of the project "C1 bezier paths on surfaces"
% by Gousenbourger et al to the manifold of PSD matrices.
% The original project is downloadable at 
% https://perso.uclouvain.be/pygousenbourger/#nt
% ------------------------------------------------------------
% Last modification: October 24, 2018
% ------------------------------------------------------------
function v = log_psd_quotient(x,y)

    % use a polar decomposition to compute the best representatve of the
    % equivalence class
    
    v = x;

    if isa(x,'cell') && isa(y,'cell')
        for i=1:length(x)
            [U,~,V] = svd(x{i}'*y{i});
            Q = U*V';
            y{i} = y{i}*Q';
            v{i} = y{i}-x{i};
        end
    elseif isa(x,'double') && isa(y,'double')
            [U,~,V] = svd(x'*y);
            Q = U*V';
            y = y*Q';
            v = y-x;
    else
        error('x and y must be same type');
    end
end
