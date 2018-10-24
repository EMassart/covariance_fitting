% DE_CASTELJAU
%       Computes the B�zier function based on the B�zier points p,
%       evaluated at times t in [0,1]. The number of B�zier points (m)
%       determines the order of the B�zier function.
%       The B�zier function is computed on a given manifold.
% 
% Input: manifold: the manifold on which the points are defined;
%        p       : the B�zier points (cell {n x 1} in a dimension
%                   [dim1,dim2]).
%        t       : the time at which the function is computed (double)
%
% Output: the B�zier curve, in a tensor [dim1 x dim 2 x m]

function y = de_casteljau_1d(p,t,discr,fig,perc,pieces)
    global_variables;
    
    if nargin < 2; error('Not enough arguments');
		elseif nargin < 3; discr = 20; fig = -1; perc = -1; pieces = -1; 
		elseif nargin < 4; fig = -1; perc = -1; pieces = -1;
    end
    assert(t <= 1 && t >= 0 , 't must be between [0,1]');
       
    np           	= length(p);
    [dim1,dim2] 	= size(p{1});
    y           	= zeros(dim1,dim2);
    
	pp = p;
	for j = 1:np-1       % For each step in the de casteljau process
		pp = geo_map(pp(1:end-1),pp(2:end),t,discr);
	end
	assert(length(pp) == 1,'The De Casteljau algorithm cannot run to the end');
	y = pp{1};
end
