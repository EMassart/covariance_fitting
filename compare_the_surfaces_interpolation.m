%% Comparison between the Bézier approach and the piecewise geodesic approaches
% for surface fitting to Antoni's data.
% Two metrics are also considered : the quotient metric from the
% factorization C = YY' and the section approach for this quotient.

% Author : E. Massart
% Version : 31/7/18

clear all; close all; clc;

% load the data
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data
[m,n] = size(train_data);
[test_data, params, patch] = data_points_all_but_external(r);      % training data
sampling = 3;

[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Oups, data have non appropriate size');
end
l_test = s1*s2;

error = zeros(l_test,12);

% Bezier methods
options_bezier_section.foot = 'data';
gamma_bezier_S_data = bezier_section(train_data, sampling, options_bezier_section);
fprintf('bezier_S1 done... \n');

options_bezier_section.foot = 'arithm';
gamma_bezier_S_ar = bezier_section(train_data, sampling, options_bezier_section);
fprintf('bezier_S2 done... \n');

options_bezier_section.foot = 'inductive';
gamma_bezier_S_ind = bezier_section(train_data, sampling, options_bezier_section);
fprintf('bezier_S3 done... \n');

gamma_bezier_G = bezier_quotient(train_data, sampling);
fprintf('bezier_G done... \n');

% check the sampling
for i_test = 1:l_test
    p = params{i_test};
    C_test =  test_data{i_test}*test_data{i_test}';
    p_new = [(p(1)-4)/3, p(2)];     % such that all the indices are on the grid: [4:1.5:13]x[0:0.5:4] -> [0:0.5:3]x[0:0.5:4]
    indx = 2*p_new + 1;             % the indices on the grid : [1:1:7]x[1:1:9]
    error(i_test,1) = norm(gamma_bezier_S_data{indx(1),indx(2)}*gamma_bezier_S_data{indx(1),indx(2)}' - C_test,'fro');
    error(i_test,2) = norm(gamma_bezier_S_ar{indx(1),indx(2)}*gamma_bezier_S_ar{indx(1),indx(2)}' - C_test,'fro');
    error(i_test,3) = norm(gamma_bezier_S_ind{indx(1),indx(2)}*gamma_bezier_S_ind{indx(1),indx(2)}' - C_test,'fro');
    error(i_test,4) = norm(gamma_bezier_G{indx(1),indx(2)}*gamma_bezier_G{indx(1),indx(2)}' - C_test,'fro');
end

d_ref = zeros(1,l_test);

% Patchwise methods
for i_test = 1:l_test
    
    C_test =  test_data{i_test}*test_data{i_test}';
    
    % define the corners of the patch
    Y = zeros(3024,r,4);
    Y(:,:,1) = train_data{patch{i_test}(1),patch{i_test}(2)};
    Y(:,:,2) = train_data{patch{i_test}(1)+1,patch{i_test}(2)};
    Y(:,:,3) = train_data{patch{i_test}(1),patch{i_test}(2)+1};
    Y(:,:,4) = train_data{patch{i_test}(1)+1,patch{i_test}(2)+1};
    
    patch_display = patch{i_test}
    p = params{i_test}
    
    p_new = [(p(1)-4)/3, p(2)]     % such that all the indices are on the grid: [4:1.5:13]x[0:0.5:4] -> [0:0.5:3]x[0:0.5:4]
    
    if p_new(1) < m-1
        weights(1) = p_new(1) - floor(p_new(1)); % such that all the indices are on the grid: [4:1.5:13]x[0:0.5:4] -> [0:0.5]x[0:0.5]
    else
        weights(1) = 1;
    end
    
    if p_new(2) < n-1
        weights(2) = p_new(2) - floor(p_new(2)); % such that all the indices are on the grid: [4:1.5:13]x[0:0.5:4] -> [0:0.5]x[0:0.5]
    else
        weights(2) = 1;
    end
    
    weights = weights
    % patchwise sectional method - based at one of the data points
    Y_projected = Y;
    for i = 2:4
        Q = orth_pol(Y(:,:,1)'*Y(:,:,i));
        Y_projected(:,:,i) = Y(:,:,i)*Q';
    end
    YS = (1-weights(1))*(1-weights(2))*Y_projected(:,:,1) + weights(1)*(1-weights(2))*Y_projected(:,:,2) + (1-weights(1))*weights(2)*Y_projected(:,:,3) + weights(1)*weights(2)*Y_projected(:,:,4);
    error(i_test,5) = norm(YS*YS' - C_test,'fro');
    
    % patchwise sectional method - based at the arithmetic mean of the data points
    Y_projected = Y;
    Y_arithm = m_arithm(Y);
    for i = 1:4
        Q = orth_pol(Y_arithm'*Y(:,:,i));
        Y_projected(:,:,i) = Y(:,:,i)*Q';
    end
    YS = (1-weights(1))*(1-weights(2))*Y_projected(:,:,1) + weights(1)*(1-weights(2))*Y_projected(:,:,2) + (1-weights(1))*weights(2)*Y_projected(:,:,3) + weights(1)*weights(2)*Y_projected(:,:,4);
    error(i_test,6) = norm(YS*YS' - C_test,'fro');
    
    % patchwise sectional method - based at the inductive mean of the data points
    Y_projected = Y;
    Y_ind = m_ind(Y);
    for i = 1:4
        Q = orth_pol(Y_ind'*Y(:,:,i));
        Y_projected(:,:,i) = Y(:,:,i)*Q';
    end
    YS = (1-weights(1))*(1-weights(2))*Y_projected(:,:,1) + weights(1)*(1-weights(2))*Y_projected(:,:,2) + (1-weights(1))*weights(2)*Y_projected(:,:,3) + weights(1)*weights(2)*Y_projected(:,:,4);
    error(i_test,7) = norm(YS*YS' - C_test,'fro');
    
    % patchwise geodesic method
    orth12 = orth_pol(Y(:,:,1)'*Y(:,:,2));
    orth34 = orth_pol(Y(:,:,3)'*Y(:,:,4));
    Y12 = (1-weights(1))*Y(:,:,1) + weights(1)*Y(:,:,2)*orth12';
    Y34 = (1-weights(1))*Y(:,:,3) + weights(1)*Y(:,:,4)*orth34';
    orth = orth_pol(Y12'*Y34);
    Ygamma = (1-weights(2))*Y12 + weights(2)*Y34*orth';
    error(i_test,8) = norm(Ygamma*Ygamma' - C_test, 'fro');
    
    d_ref(i_test) = 0;
    for i = 1:4
        d_ref(i_test) = d_ref(i_test) + norm(Y(:,:,i)*Y(:,:,i)' - C_test, 'fro')^2 / 4;
    end
        
end

    

