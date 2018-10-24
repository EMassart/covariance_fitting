% bijmn(p, Am, An, i, j, m, n, d): 
% 	Computes the control point b_ij^mn
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	p : [cell] The interpolation points.
% 				(M x N)
%
% 			Am, An : [matrix] a matrix defining the linear combination.
%
%			m,n : [int] patch
%
% 			i,j : [int] position of the control point within the patch.
% 				 	 	can be 0, 1, 2 or 3.
%
% 			d : [int] for the linear combination
%
% OUTPUT: 	b : [marrix] the control point.
%
% ------------------------------------------------------------
% Authors:
% 	Pierre-Yves Gousenbourger
% 	Paul Striewski
%
% Versions
% 	04/11/2015: First version.
% ------------------------------------------------------------

function b = bijmn(p, Am, An, i, j, m, n, d)
	
	global_variables;
	
	% Parameters
	[M,N] = size(p);
	op 	 = [0,0,1,1];
	pp 	 = cell(size(p));	% Projected p
	
	% Reference point and projection to its tangent space.
	pref = p{m+op(i+1),n+op(j+1)};
	xbound = max(1,m-d-1):min(M,m+d+2);
	ybound = max(1,n-d-1):min(N,n+d+2);
	
	
	for k = xbound
	for l = ybound
		pp{k,l} = geo_log(pref,p{k,l}); 
	end
	end
	
	% Step 1:
	% -------
	% compute [beta_(m,n-d:n+d)]
	% compute [beta_(m,n-d+1:n+d+1)]
	% compute [beta_(m,n-d+1:n+d+1)]
	% compute [beta_(m,n-d+1:n+d+1)]
	%
	% Step 2:
	% ------- 
	% use this to compute
	% a(m,n), a(m,n+1), a(m+1,n), a(m+1,n+1)
	beta = cell(N,1);
	alpha = cell(2,2);
	for s = 0:1
	for r = 0:1
		bound = max(1,n-d+r-1):min(N,n+d+r+1);
		% Step 1
		for k = bound
			beta{k} = compute_alpha(pp(:,k), Am, m+s, d);
		end
		% Step 2
		alpha{s+1,r+1} = compute_alpha(beta, An, n+r, d);
	end
	end
	
	bij = 	(((3-i)/3)*((3-j)/3))	.*alpha{1,1} ...
		+ 	(((3-i)/3)*(j/3))		.*alpha{1,2} ...
		+ 	((i/3)*((3-j)/3))		.*alpha{2,1} ...
		+ 	((i/3)*(j/3))			.*alpha{2,2} ;
		
	b 	= geo_exp(pref, bij);
end
