%% This script is the main script to reproduce the results regarding optimization over a Bezier surface 

% (This Bezier surface was built based on a section of the quotient manifold 
% (in contrast with the approach using directly the exponential and log maps for
% the quotient geometry).

% Three possibilities are considered for the choice of the section:
%   - data : the section is based at one of the points 
%   - arithm : the section is based at the arithmetic mean of the points.
%   Since the arithmetic mean of a set of PSD matrices of rank r has
%   usually a rank larger than r, we truncate the rank afterwards.
%   - inductive : the section is based at the inductive mean of the points.

% This code actually does a bit more than running SD to optimize over the
% Bezier surface. A sanity check is also provided, to check the consistency
% of the results: we compute the Bezier surface on a grid and compute, for
% each test point, the minimum of the distance over the grid. 
% The parts of the code that are associated with this sanity check are
% between dashed comment lines, please comment them if you don't want to
% use the sanity check. 
% The sanity check can roughly double the computation time, passing from 5
% hours on a Windows 7 platform with 8 cores at 3.60 GHz and 16 GB ram to 8
% hours.

% Author : E. Massart
% Last modification: October 24, 2018

clear all; close all; clc;
set_path;

% Load the training data
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data


%-----------------------------------------------------------------------------------------------------------------------------------
% % Compute the Bezier surface on these data (not required to solve the
% % optimization problem)
% sampling = 11;
% options.foot = 'data';
% gamma_bezier_S_first = bezier_section(train_data, sampling, options);
% options.foot = 'arithm';
% gamma_bezier_S_arithm = bezier_section(train_data, sampling, options);
% options.foot = 'inductive';
% gamma_bezier_S_ind = bezier_section(train_data, sampling, options);
% fprintf('Piecewise geodesic surfaces constructed :-) \n');
%-----------------------------------------------------------------------------------------------------------------------------------


% Load the test data
[test_data, params] = data_points_all_but_external(r);      % training data
[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Data have non appropriate size :-/');
end
l_test = s1*s2;


%-----------------------------------------------------------------------------------------------------------------------------------
% [m_tot, n_tot] = size(gamma_bezier_S_first);
% error_grid = zeros(3, m_tot, n_tot, l_test);
% error_grid_rec = zeros(3, l_test);
%-----------------------------------------------------------------------------------------------------------------------------------


error_SD_rec = zeros(3, l_test);
sol_SD_rec = zeros(2, 3, l_test);

for i_test = 1:l_test
    fprintf('------------------------------------------------------- TEST NUMBER %d \n', i_test);
    
    C_ref = test_data{i_test}*test_data{i_test}';
    
    
    %-----------------------------------------------------------------------------------------------------------------------------------
%     for i_eval = 1:m_tot
%         fprintf('i_eval = %d \n', i_eval);
%         for j_eval = 1:n_tot
%             C = gamma_bezier_S_first{i_eval,j_eval}*gamma_bezier_S_first{i_eval,j_eval}';
%             error_grid(1, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;
%             C = gamma_bezier_S_arithm{i_eval,j_eval}*gamma_bezier_S_arithm{i_eval,j_eval}';
%             error_grid(2, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;
%             C = gamma_bezier_S_ind{i_eval,j_eval}*gamma_bezier_S_ind{i_eval,j_eval}';
%             error_grid(3, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;                    
%         end
%     end
%     error_grid_rec(1,i_test) = min(min(error_grid(1,:,:,i_test)));
%     error_grid_rec(2,i_test) = min(min(error_grid(2,:,:,i_test)));
%     error_grid_rec(3,i_test) = min(min(error_grid(3,:,:,i_test)));
%     fprintf('---------------------------------- Errors grid : %4.2e, %4.2e, %4.2e \n', error_grid_rec(1,i_test), error_grid_rec(2,i_test), error_grid_rec(3,i_test));
    %-----------------------------------------------------------------------------------------------------------------------------------

    
    fprintf('---------------------------------- Run SD data \n');
    options.foot = 'data';
    [x_end, info] = steepest_descent_for_bezier_section(train_data, C_ref, options);
    error_SD_rec(1,i_test) = info.f(info.k);
    sol_SD_rec(:,1,i_test) = x_end;
    
    fprintf('---------------------------------- Run SD arithm \n');
    options.foot = 'arithm';
    [x_end, info] = steepest_descent_for_bezier_section(train_data, C_ref, options);
    error_SD_rec(2,i_test) = info.f(info.k);
    sol_SD_rec(:,2,i_test) = x_end;
    
    fprintf('---------------------------------- Run SD inductive \n');
    options.foot = 'inductive';
    [x_end, info] = steepest_descent_for_bezier_section(train_data, C_ref, options);
    error_SD_rec(3,i_test) = info.f(info.k);
    sol_SD_rec(:,3,i_test) = x_end;    
    
end

save('error_bezier_section_SD.mat', 'error_SD_rec', 'sol_SD_rec', 'r');


%-----------------------------------------------------------------------------------------------------------------------------------
% save('error_bezier_section_Grid.mat', 'error_grid', 'error_grid_rec', 'sampling');
%-----------------------------------------------------------------------------------------------------------------------------------






