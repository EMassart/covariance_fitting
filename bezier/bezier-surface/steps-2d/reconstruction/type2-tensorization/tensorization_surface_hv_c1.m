%TENSORIZATION_SURFACE Compute a Bezier surface using the tensorisation
%method: Horizontal Bezier curves are drawn through the rows of the control
%point grid. These curves are evaluated at a time t1, as a result, a list
%of points is obtained. A vertical Bezier curve is drawn through these
%points and evaluated in t2. This gives the point (t1, t2) on the surface.

% Input: b       : the control points (cell d+1 x d+1).
%        t1,t2   : the evaluated time on the square [0,1]x[0,1].
%        d       : the degree of the curve.
%
% Output: the Bezier curve, in a matrix [dim1 x dim2]

function y = tensorization_surface_hv_c1( b,d,t1,t2,n,m,nMax,mMax)

horizontal_curve = cell(d+1,1);

for i=1:d+1
    current_points = b(:,i);
    
    % First and last point starting iteration
	if n > 1
		first = de_casteljau_1d([current_points(1) current_points(2)],1/2);
		current_points{1} = first;
	end
	if n < nMax	
		last  = de_casteljau_1d([current_points(end-1) current_points(end)],1/2);
    	current_points{end} = last;
    end
    
    % Generate the fictionnal control points
    horizontal_curve{i} = de_casteljau_1d(current_points, t1);
end


% First and last point starting iteration
if m > 1
	first = de_casteljau_1d([horizontal_curve(1) horizontal_curve(2)],1/2);
	horizontal_curve{1} = first;
end
if m < mMax	
	last  = de_casteljau_1d([horizontal_curve(end-1) horizontal_curve(end)],1/2);
	horizontal_curve{end} = last;
end
y = de_casteljau_1d(horizontal_curve, t2);

end

