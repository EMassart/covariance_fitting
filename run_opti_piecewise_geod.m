clear all; close all; clc;
set_path;

% Compute the optimal results for the piecewise bilinear geodesic and
% section approaches.
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data
sampling = 11;
[test_data, params] = data_points_all_but_external(r);      % training data

[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Oups, data have non appropriate size');
end
l_test = s1*s2;
% [m_tot, n_tot] = size(test_data{1,1});

t1 = linspace(0,1,11);
t2 = linspace(0,1,11);
error_grid = zeros(4, length(t1), length(t2), l_test);
error_grid_rec = zeros(4, l_test);
error_SD_rec = zeros(4, l_test);
sol_SD_rec = zeros(2, 4, l_test);


% information about the data, required to deduce the patches
heading_type1 = 0.5:1:3.5;
magn_type1 = 4:1.5:13;

heading_type2 = 0:1:4;
magn_type2 = 5.5:3:11.5;

n1 = length(heading_type1);
m1 = length(magn_type1);

n2 = length(heading_type2);
m2 = length(magn_type2);

patch = cell(1,43);
patch(1:3) = {[1,1]};
patch(4:5) = {[1,2]};
patch(6:7) = {[1,3]};
patch(8:10) = {[2,1]};
patch(11:12) = {[2,2]};
patch(13:14) = {[2,3]};
patch(15:17) = {[3,1]};
patch(18:19) = {[3,2]};
patch(20:21) = {[3,3]};
patch(22:24) = {[4,1]};
patch(25:26) = {[4,2]};
patch(27:28) = {[4,3]};
patch(29) = {[1,1]};
patch(30) = {[1,2]};
patch(31) = {[1,3]};
patch(32) = {[2,1]};
patch(33) = {[2,2]};
patch(34) = {[2,3]};
patch(35) = {[3,1]};
patch(36) = {[3,2]};
patch(37) = {[3,3]};
patch(38) = {[4,1]};
patch(39) = {[4,2]};
patch(40) = {[4,3]};
patch(41) = {[4,1]};
patch(42) = {[4,2]};
patch(43) = {[4,3]};

m = zeros(1,l_test);
for i_test = 1:l_test
    fprintf('------------------------------------------------------- TEST NUMBER %d \n', i_test);
    
    C_ref = test_data{i_test}*test_data{i_test}';
    % select the patch and build the sur
    
    Y = zeros(3024,r,4);
    Y(:,:,1) = train_data{patch{i_test}(2),patch{i_test}(1)};
    Y(:,:,2) = train_data{patch{i_test}(2)+1,patch{i_test}(1)};
    Y(:,:,3) = train_data{patch{i_test}(2),patch{i_test}(1)+1};
    Y(:,:,4) = train_data{patch{i_test}(2)+1,patch{i_test}(1)+1};
    
    % piecewise geodesic interpolation of the corners of the patch
    orth12 = orth_pol(Y(:,:,1)'*Y(:,:,2));
    orth34 = orth_pol(Y(:,:,3)'*Y(:,:,4));
    
    for i = 1:length(t1)
        for j = 1:length(t2)
            Y12 = (1-t1(i))*Y(:,:,1) + t1(i)*Y(:,:,2)*orth12';
            Y34 = (1-t1(i))*Y(:,:,3) + t1(i)*Y(:,:,4)*orth34';
            orth = orth_pol(Y12'*Y34);
            Ygamma = (1-t2(j))*Y12 + t2(j)*Y34*orth';
            error_grid(1,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
        end
    end
    
    % bilinear interpolation in the section
    
    % section based on the lower_left anchor point
    Y_section = zeros(size(Y));
    Y_anchor_section = Y(:,:,1);
    Y_section(:,:,1) = Y(:,:,1); 
    for i_loc = 2:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
    for i = 1:length(t1)
        for j = 1:length(t2)
            Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
            error_grid(2,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
        end
    end
    
    Y_section = zeros(size(Y));
    Y_anchor_section =  m_arithm(Y);
    for i_loc = 1:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
    for i = 1:length(t1)
        for j = 1:length(t2)
            Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
            error_grid(3,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
        end
    end
    
        
    Y_section = zeros(size(Y));
    Y_anchor_section =  m_ind(Y);
    for i_loc = 1:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
    for i = 1:length(t1)
        for j = 1:length(t2)
            Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
            error_grid(4,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
        end
    end
    
            
    for i = 1:4
        error_grid_rec(i,i_test) = min(min(error_grid(i,:,:,i_test)));
    end
        
    fprintf('------------------------------------------ Errors on grid: %4.2e, %4.2e, %4.2e, %4.2e \n',  error_grid_rec(1,i_test),  error_grid_rec(2,i_test),  error_grid_rec(3,i_test),  error_grid_rec(4,i_test));

    %     n = zeros(1,4);
%     n(1) = norm(C_ref - Y(:,:,1)*Y(:,:,1)','fro')^2;
%     n(2) = norm(C_ref - Y(:,:,2)*Y(:,:,2)','fro')^2;
%     n(3) = norm(C_ref - Y(:,:,3)*Y(:,:,3)','fro')^2;
%     n(4) = norm(C_ref - Y(:,:,4)*Y(:,:,4)','fro')^2;
%     n = n
%     m(i_test) = mean(n);
%     
    % implementation of the two optimization algorithms
    [~, info] = variable_projection_final(Y, C_ref);
    sol_SD_rec(:, 1, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(1, i_test) = info.cost(info.nIter);
    fprintf('----------Variable projection finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(1,i_test), sol_SD_rec(1, 1, i_test), sol_SD_rec(2, 1, i_test));
    
    options.input = 'first';
    [~, info] = sectional_approach_final(Y, C_ref, options);
    sol_SD_rec(:, 2, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(2, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach first finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(2,i_test), sol_SD_rec(1, 2, i_test), sol_SD_rec(2, 2, i_test));
    
    options.input = 'arithm';
    [~, info] = sectional_approach_final(Y, C_ref, options);
    sol_SD_rec(:, 3, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(3, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach arithm finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(3,i_test), sol_SD_rec(1, 3, i_test), sol_SD_rec(2, 3, i_test))
    
    options.input = 'inductive';
    [~, info] = sectional_approach_final(Y, C_ref, options);
    sol_SD_rec(:, 4, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(4, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach ind finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(4,i_test), sol_SD_rec(1, 4, i_test), sol_SD_rec(2, 4, i_test));
   
end

save('Results_piecewise_geod.mat', 'error_grid', 'error_grid_rec', 'error_SD_rec', 'sol_SD_rec', 'params', 'patch');
    