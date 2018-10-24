% control_points_simple_generation_2d(pb): 
% 		 Computes the control points for a Bezier surface on M,
% 		 with the property that the curves are C^2 on R^n
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	pb : [struct] The structure of the problem containing:
% 					- pb.interp
% 					- pb.m
% 					- pb.n
%
% OUTPUT: 	pb : [struct] The structure is updated with pb.control
% ------------------------------------------------------------
% Authors:
% 	Pierre-Yves Gousenbourger
% 	Paul Striewski
%
% Versions
% 	04/11/2015: First version.
% ------------------------------------------------------------

function pb = control_points_simple_generation_2d(pb)

	global_variables;
	
	% Data
	p = pb.interp;
	M = pb.n;
	N = pb.m;
	
	assert(M > 2, 'There must be at least 3 interpolation points in the x direction.');
	assert(N > 2, 'There must be at least 3 interpolation points in the y direction.');

	% Parameters
	d  = max(M,N); % the number of elements around the diag of Am, An.
	
	Am = inv((1/6)*(diag(4*ones(1,M-2),0) + diag(ones(1,M-3),1) + diag(ones(1,M-3),-1)));
	An = inv((1/6)*(diag(4*ones(1,N-2),0) + diag(ones(1,N-3),1) + diag(ones(1,N-3),-1)));

	% Preparation of structures
	b 	 	= cell(3*M - 2, 3*N - 2);
	
	
	
	% ==== Computation in 2d ===========================================
	% Interior points
	for m = 1:M-1
	for n = 1:N-1
		for i = 0:3
		for j = 0:3
			b{3*(m-1) + i + 1, 3*(n-1) + j + 1} = bijmn(p, Am, An, i, j, m, n, d);
		end
		end
	end
	end
	% ==== End computation in 2d =======================================
	
	
	%Interpolation points
	b(1:3:end,1:3:end) = p;
	
	% Control points at the interfaces
	% x direction
	for m = 1:M-2
		b(3*m+1, 2:3:end) = geo_map(b(3*m, 2:3:end), b(3*m + 2, 2:3:end), 0.5);
		b(3*m+1, 3:3:end) = geo_map(b(3*m, 3:3:end), b(3*m + 2, 3:3:end), 0.5);
	end
		
	% y direction
	for n = 1:N-2
		b(2:3:end,3*n+1) = geo_map(b(2:3:end,3*n), b(2:3:end,3*n+2), 0.5);
		b(3:3:end,3*n+1) = geo_map(b(3:3:end,3*n), b(3:3:end,3*n+2), 0.5);
	end	
		
	pb.control 	= b;
	pb.d 	 	= 3;

end
