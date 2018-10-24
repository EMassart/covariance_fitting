%% CURVE-RECONSTRUCTION-DOUBLE-BEZIER
%       Reconstrucs the Bezier curve based on the interpolation points and 
%       the control points given in the structure pb.
%
%  INPUT: * pb : the structure of the interpolation problem containing
% 		  * method : the method used (1 = HV, 2 = VH)
% 
%  OUTPUT: pb : the structure updated with the curve pb.curve stored in a
%  tensor.

function pb = curve_reconstruction_double_bezier_c1(pb,method)
	global_variables;
	
    % data
    b      	= pb.control;
    t       = pb.t;
    d       = pb.d;
    n 		= pb.n;
    m 		= pb.m;
    
    % The curve is a different object following of the manifold.
    vt 		= linspace(0, 1, t);
    curve 	= zeros(t*(n-1),t*(m-1),pb.dim1,pb.dim2); 	% Curve is a tensor
														% n x m (rectangle parametrization)
														% dim1xdim2 (points)
	
	% Reconstruction function
	switch method
		case 1
			% Horizontal-vertical
			make_curve = @tensorization_surface_hv_c1;
			methodused = 'type2-tensorization-hv-c1';
		case 2
			% Vertical-horizontal
			make_curve = @tensorization_surface_vh_c1;
			methodused = 'type2-tensorization-vh-c1';
		otherwise
			error('wrong method for tensorization');
	end
	
	% Waitbar just because it is fun ^^.
		h = waitbar(0,'Initialization of the De Casteljau algoritm...');
		total_steps = t*t*(n-1)*(m-1);
		perc = 0;
		
	% For each patch (n ~ x, m ~ y)
	for k = 1:n-1
	for l = 1:m-1
		% forward
		dx = (k-1)*t;
		dy = (l-1)*t;
		
		% -- Control points reduction
		% --------------------------------
		
		% The control points are actually inner control points and outer control points
		% ind_x and ind_y give the range of points needed to be taken in b.
		% x
		if k == 1
			ind_x = [(k-1)*d + 1: k*d , k*d + 2]; 
		elseif k == n-1
			ind_x = [(k-1)*d , (k-1)*d + 2: k*d+1]; 
		else
			ind_x = [(k-1)*d , (k-1)*d + 2: k*d , k*d + 2]; 
		end
		% y
		if l == 1
			ind_y = [(l-1)*d + 1: l*d , l*d + 2];
		elseif l == m-1
			ind_y = [(l-1)*d , (l-1)*d + 2: l*d + 1];
		else 
			ind_y = [(l-1)*d , (l-1)*d + 2: l*d , l*d + 2];
		end
		
		% If only one patch
		if m == 2
			ind_y = 1:d+1;
		end
		if n == 2;
			ind_x = 1:d+1;
		end
		
		
		% reduction
		bb = b( ind_x ,  ind_y );
		
		% --- Reconstruction by overall mean
		% ---------------------------------
		
		for i=1:t
		for j=1:t
			curve(dx+i,dy+j,:,:) = make_curve(bb,d,vt(i),vt(j),k,l,n-1,m-1);

			% Waitbar
			perc = perc + 1/total_steps;
			waitbar(perc,h,sprintf(...
				['Processing De Casteljau (patch [',num2str(k),',',num2str(l),'] of [',num2str(n-1),',',num2str(m-1),']). Complete: %.f %%.'],...
				perc*100));
		end
		end
	end
	end
    
    close(h);

    %% Store
    pb.curve = curve;
    pb.method = methodused;
end
