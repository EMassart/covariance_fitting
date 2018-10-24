function f = cost_function_bezier(pb, t1, t2)
% This function reconstructs the value of the Bezier surface for the couple
% of parameters (t_1, t_2), and compares it with the matrix pb.C

	global_variables

    % data
    b      	= pb.control;
    d       = pb.d;
    n 		= pb.n;
    m 		= pb.m;
    
    k = floor(t1)+1;
    l = floor(t2)+1;
    
    if t1 == n-1
        k = t1;
    end
    
    if t2 == m-1
        l = t2;
    end
		
    % -- Control points reduction
    % --------------------------------
    
    % The control points are actually inner control points and outer control points
    % ind_x and ind_y give the range of points needed to be taken in b.
    % x
    if k == 1
        ind_x = [((k-1)*d + 1):(k*d) , (k*d + 2)];
    elseif k == n-1
        ind_x = [((k-1)*d), ((k-1)*d + 2):(k*d+1)];
    else
        ind_x = [((k-1)*d), ((k-1)*d + 2):(k*d) , (k*d + 2)];
    end
    % y
    if l == 1
        ind_y = [((l-1)*d + 1): (l*d), (l*d + 2)];
    elseif l == m-1
        ind_y = [((l-1)*d), ((l-1)*d + 2):(l*d + 1)];
    else
        ind_y = [((l-1)*d), ((l-1)*d + 2):(l*d), (l*d + 2)];
    end
		
    % If only one patch
    if m == 2
        ind_y = 1:d+1;
    end
    if n == 2;
        ind_x = 1:d+1;
    end
		
    bb = b( ind_x ,  ind_y );
		
    % --- Reconstruction by overall mean
    % ---------------------------------

    curve = tensorization_surface_hv_c1(bb,d,t1-k+1,t2-l+1,k,l,n-1,m-1);
 

% return the value of f

    f = norm(curve*curve' - pb.C,'fro')^2;

end