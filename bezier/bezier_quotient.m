function gamma = bezier_quotient(data, sampling)
% Run the Bézier algorithm to fit a suface to the data points. We are here
% using the quotient geometry C = YY', where the total space is endowed with
% the Euclidean metric.
% 
% Input: - data : cell array containing the Y factors of the covariance 
%                 matrices to be used as data.
%        - sampling : sampling at which the surface will be evaluated
% 
% Author: E. Massart
% Last modification: October 24, 2018

manifold = 'psd_quotient';
global_variables
geo_functions
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


