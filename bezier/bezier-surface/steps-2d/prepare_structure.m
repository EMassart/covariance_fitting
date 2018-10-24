%% PREPARE_STRUCTURE -
%       Prepares the structure PB of an interpolation problem based on an
%       input DATA and a manifold MANIFOLD.
%
%  INPUT: * data : the points of interpolation;
%         * manifold: the manifold on which interpolation is done.
%         * sampling: the number of samples in the curve;
% 		  * sample: the type of data points
% 		  * nint: the number of intermediary points within the method of averaging
% 		  * object: SO3 object type
%
%         Data points must be given in a {n x 1} cell. Each element of the
%         cell is:
%           euclidean   |  vector of dimension dim
%           spherical   |  normalized vector of dimension dim
%           SO(dim)     |  orthogonal matrices of dimension dim x dim
%           shapes      |  matrices of dimension [dim x m], where dim is
%                       |   the dimension of the shape and m the
%                       |   discretization.
% 
%  OUTPUT: pb : the structure of the problem. It is composed of:
%           pb.interp         |  interpolation points in cells
%           pb.manifold       |  the manifold
%           pb.n              |  the number of interpolation points
%           pb.dim            |  the dimension of the manifold
%           pb.control        |  the control points of the path
%           pb.velocity       |  the norm of the velocity of the curve
%           pb.velocities     |  the velocities direction on internal data
%           pb.velocitytype   |  the type of velocity used (orthog or
%                             |   parallel)
%           pb.acceleration   |  the norm of the acceleration of the curve
%           pb.t              |  the discretization of each bezier segment
%                             |   (default: 50)
%           pb.disp           |  boolean for displaying of not certain
%                             |   messages
% 			pb.method         |  the method

function pb = prepare_structure(data, manifold, sampling, sample, nint, object)
	if nargin < 6
		method = 1;
		if nargin < 5
			nint = 0;
			if nargin < 4
				sample = '';
			end
		end
	end
    
    
    %dimension
    [dim1,dim2] = size(data{1,1});
     
    % Structure
    pb.interp       = data;				% Interpolation points
    pb.type 		= sample;
    pb.manifold     = manifold;			% Manifold
    pb.n 			= size(data,1);		% Number of interp (t1)
    pb.m 			= size(data,2);		% Number of interp (t2)
    pb.dim1         = dim1;				% Dimension of the space
    pb.dim2 		= dim2;
    pb.d 			= [];				% Degree of the curve
    pb.t            = sampling;
    if exist('nint','var')
        pb.nint 		= nint;
    end
    if exist('object','var')
        pb.objectSO3	= object;
    end
end
