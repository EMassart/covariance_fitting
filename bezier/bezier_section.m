function gamma = bezier_section(data, sampling, options)
% Run the B�zier algorithm to fit a suface to the data points. We are here
% using Bezier interpolation in a section of the quotient manifold C = YY'.
% 
% Input: - data : cell array containing the Y factors of the covariance 
%                 matrices to be used as data.
%        - sampling : sampling at which the surface will be evaluated
%
% Three possibilities are then considered for the choice of the section:
%   - data : the section is based at one of the points 
%   - arithm : the section is based at the arithmetic mean of the points.
%   Since the arithmetic mean of a set of PSD matrices of rank r has
%   usually a rank larger than r, we truncate the rank afterwards.
%   - inductive : the section is based at the inductive mean of the points.
% 
% Author: E. Massart
% Last modification: October 24, 2018


% definition of the parameters
[m,n] = size(data);
manifold = 'euclidean';         % since in the section we are working as in an euclidean space
global_variables
geo_functions

% project on the section
if strcmp(options.foot, 'data')
    Y_anchor_section = data{1,1};
    for i = 1:m
        for j = 1:n
            if (i~=1 || j~=1)
                Q = orth_pol(Y_anchor_section'*data{i,j});
                data{i,j} = data{i,j}*Q';
            end
        end
    end
elseif strcmp(options.foot, 'arithm')
    Y_anchor_section = m_arithm(data);
    for i = 1:m
        for j = 1:n
            Q = orth_pol(Y_anchor_section'*data{i,j});
            data{i,j} = data{i,j}*Q';
        end
    end
elseif strcmp(options.foot, 'inductive')
    Y_anchor_section = m_ind(data);
    for i = 1:m
        for j = 1:n
            Q = orth_pol(Y_anchor_section'*data{i,j});
            data{i,j} = data{i,j}*Q';
        end
    end
else
    fprintf('Warning: no such choice of the section defined \n');
end


problem  = prepare_structure(data,manifold,sampling);

% compute control points
problem = control_points_simple_generation_2d(problem);

% reconstruct the surface
problem = curve_reconstruction_double_bezier_c1(problem,1); % horizontal-vertical

% recover the curve
gamma_mat = problem.curve;

% remove uninformative rows and columns
s = size(gamma_mat);

indx_to_remove_x = sampling+1 : sampling : s(1);
indx_to_remove_y = sampling+1 : sampling : s(2);

for i = 1:length(indx_to_remove_x)
    gamma_mat(indx_to_remove_x(length(indx_to_remove_x) - i+1),:,:,:) = [];
end
for j = 1:length(indx_to_remove_y)
    gamma_mat(:,indx_to_remove_y(length(indx_to_remove_y) - j+1),:,:) = [];
end

[m_out, n_out, ~, ~] = size(gamma_mat);
gamma = cell(m_out,n_out);
s = size(data{1,1});
for i = 1:m_out
    for j = 1:n_out
        gamma{i,j} = reshape(gamma_mat(i,j,:,:),s(1),s(2));
    end
end


end


