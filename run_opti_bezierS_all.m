
% Compares the cost function when fitting the surface with Bézier, with the
% cost function obtained when fitting the surface by patches


clear all; close all; clc;
set_path;
rng(0);

% load the data
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data
sampling = 11;

options.input = 'first';
gamma_bezier_S_first = bezier_S(train_data, sampling, options);
options.input = 'arithm';
gamma_bezier_S_arithm = bezier_S(train_data, sampling, options);
options.input = 'inductive';
gamma_bezier_S_ind = bezier_S(train_data, sampling, options);

fprintf('Piecewise geodesic surfaces constructed :-) \n');

[test_data, params] = data_points_all_but_external(r);      % training data

[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Oups, data have non appropriate size');
end
l_test = s1*s2;
[m_tot, n_tot] = size(gamma_bezier_S_first);
error_grid = zeros(3, m_tot, n_tot, l_test);
error_grid_rec = zeros(3, l_test);
error_SD_rec = zeros(3, l_test);
sol_SD_rec = zeros(2, 3, l_test);

for i_test = 1:l_test
    fprintf('------------------------------------------------------- TEST NUMBER %d \n', i_test);
    
    C_ref = test_data{i_test}*test_data{i_test}';
    for i_eval = 1:m_tot
        fprintf('i_eval = %d \n', i_eval);
        for j_eval = 1:n_tot
            C = gamma_bezier_S_first{i_eval,j_eval}*gamma_bezier_S_first{i_eval,j_eval}';
            error_grid(1, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;
            C = gamma_bezier_S_arithm{i_eval,j_eval}*gamma_bezier_S_arithm{i_eval,j_eval}';
            error_grid(2, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;
            C = gamma_bezier_S_ind{i_eval,j_eval}*gamma_bezier_S_ind{i_eval,j_eval}';
            error_grid(3, i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;                    
        end
    end
    
    error_grid_rec(1,i_test) = min(min(error_grid(1,:,:,i_test)));
    error_grid_rec(2,i_test) = min(min(error_grid(2,:,:,i_test)));
    error_grid_rec(3,i_test) = min(min(error_grid(3,:,:,i_test)));
    fprintf('---------------------------------- Errors grid : %4.2e, %4.2e, %4.2e \n', error_grid_rec(1,i_test), error_grid_rec(2,i_test), error_grid_rec(3,i_test));

    fprintf('---------------------------------- Run SD first \n');
    options.input = 'first';
    [x_end, info] = steepest_descent_for_bezierS(train_data, C_ref, options);
    error_SD_rec(1,i_test) = info.f(info.k);
    sol_SD_rec(:,1,i_test) = x_end;
    
    fprintf('---------------------------------- Run SD arithm \n');
    options.input = 'arithm';
    [x_end, info] = steepest_descent_for_bezierS(train_data, C_ref, options);
    error_SD_rec(2,i_test) = info.f(info.k);
    sol_SD_rec(:,2,i_test) = x_end;
    
    fprintf('---------------------------------- Run SD inductive \n');
    options.input = 'inductive';
    [x_end, info] = steepest_descent_for_bezierS(train_data, C_ref, options);
    error_SD_rec(3,i_test) = info.f(info.k);
    sol_SD_rec(:,3,i_test) = x_end;    
    
end

save('error_bezier_S_other_points.mat', 'error_grid', 'error_grid_rec', 'error_SD_rec', 'sol_SD_rec');