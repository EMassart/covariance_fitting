
% Compares the cost function when fitting the surface with Bézier, with the
% cost function obtained when fitting the surface by patches


clear all; close all; clc;
rng(0);

% load the data
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data
sampling = 11;
gamma_bezier_G = bezier_G(train_data, sampling);
fprintf('bezier_G done... \n');

[test_data, params] = data_points_all_but_external(r);      % training data

[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Oups, data have non appropriate size');
end
l_test = s1*s2;
[m_tot, n_tot] = size(gamma_bezier_G);
error_grid = zeros(m_tot, n_tot, l_test);
error_grid_rec = zeros(1, l_test);
error_SD_rec = zeros(1, l_test);
sol_SD_rec = zeros(2, l_test);

for i_test = 1:l_test
    fprintf('------------------------------------------------------- TEST NUMBER %d \n', i_test);
    
    C_ref = test_data{i_test}*test_data{i_test}';
    for i_eval = 1:m_tot
        fprintf('i_eval = %d \n', i_eval);
        for j_eval = 1:n_tot
            C = gamma_bezier_G{i_eval,j_eval}*gamma_bezier_G{i_eval,j_eval}';
            error_grid(i_eval, j_eval, i_test) = norm(C - C_ref,'fro')^2;
        end
    end
    error_grid_rec(i_test) = min(min(error_grid(:,:,i_test)));
    fprintf('Best value obtained from grid = %4.2e \n', error_grid_rec(i_test));
    [x_end, info] = steepest_descent_for_bezier(train_data, C_ref, struct());
    error_SD_rec(i_test) = info.f(info.k);
    sol_SD_rec(:,i_test) = x_end';
end

save('error_bezier_Oct6.mat', 'error_grid', 'error_grid_rec', 'error_SD_rec', 'sol_SD_rec', 'r', 'sampling');
