%% This script is the main script to reproduce the results regarding optimization over a piecewise bilinear surface

% This code performs optimization over a patchwise bilinear surface. The
% code handle both the generalization of biliner interpolation to the
% manifold case, using successive geodesics in the quotient, and the
% bilinear interpolation in the section.

% Three possibilities are then considered for the choice of the section:
%   - data : the section is based at one of the points 
%   - arithm : the section is based at the arithmetic mean of the points.
%   Since the arithmetic mean of a set of PSD matrices of rank r has
%   usually a rank larger than r, we truncate the rank afterwards.
%   - inductive : the section is based at the inductive mean of the points.

% For each test point, the code returns the closest point in the patch of
% the surface that is associated to the test_point.

% This code actually does a bit more than running optimization algorithms
% to optimize over the surface. A sanity check is also provided, to check 
% the consistency of the results: we compute the Bezier surface on a grid and
% compute, for each test point, the minimum of the distance over the grid. 
% The parts of the code that are associated with this sanity check are
% between dashed comment lines, please comment them if you don't want to
% use the sanity check. 
% The sanity check can roughly double the computation time, passing from 3
% hours on a Windows 7 platform with 8 cores at 3.60 GHz and 16 GB ram to 5
% hours.

% Author : E. Massart
% Last modification: October 24, 2018

clear all; close all; clc;
set_path;

% Compute the optimal results for the piecewise bilinear geodesic and
% section approaches.
r = 20;                                          % rank to which the data are truncated
train_data = data_points_external_frame(r);      % training data
[test_data, params, patch] = data_points_all_but_external(r);      % training data

[s1, s2] = size(test_data);
if s1 ~= 1 && s2 ~= 1
    disp('Oups, data have non appropriate size');
end
l_test = s1*s2;


%START SANITY CHECK PART-----------------------------------------------------------------------------------------------------------------------------------
% t1 = linspace(0,1,11);
% t2 = linspace(0,1,11);
% error_grid = zeros(4, length(t1), length(t2), l_test);
% error_grid_rec = zeros(4, l_test);
%STOP SANITY CHECK PART------------------------------------------------------------------------------------------------------------------------------------

error_SD_rec = zeros(4, l_test);
sol_SD_rec = zeros(2, 4, l_test);

m = zeros(1,l_test);
for i_test = 1:l_test
    fprintf('------------------------------------------------------- TEST NUMBER %d \n', i_test);
    
    C_ref = test_data{i_test}*test_data{i_test}';
    
    % select the corners of the patch
    Y = zeros(3024,r,4);
    Y(:,:,1) = train_data{patch{i_test}(2),patch{i_test}(1)};
    Y(:,:,2) = train_data{patch{i_test}(2)+1,patch{i_test}(1)};
    Y(:,:,3) = train_data{patch{i_test}(2),patch{i_test}(1)+1};
    Y(:,:,4) = train_data{patch{i_test}(2)+1,patch{i_test}(1)+1};
    
    %START SANITY CHECK PART-----------------------------------------------------------------------------------------------------------------------------------

%     % interpolation based on geodesics
%     orth12 = orth_pol(Y(:,:,1)'*Y(:,:,2));
%     orth34 = orth_pol(Y(:,:,3)'*Y(:,:,4));
%     
%     for i = 1:length(t1)
%         for j = 1:length(t2)
%             Y12 = (1-t1(i))*Y(:,:,1) + t1(i)*Y(:,:,2)*orth12';
%             Y34 = (1-t1(i))*Y(:,:,3) + t1(i)*Y(:,:,4)*orth34';
%             orth = orth_pol(Y12'*Y34);
%             Ygamma = (1-t2(j))*Y12 + t2(j)*Y34*orth';
%             error_grid(1,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
%         end
%     end
%     
%     % bilinear interpolation in the section
%         % section based on the lower_left anchor point
%     Y_section = zeros(size(Y));
%     Y_anchor_section = Y(:,:,1);
%     Y_section(:,:,1) = Y(:,:,1); 
%     for i_loc = 2:4
%         Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
%         Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
%     end
%     for i = 1:length(t1)
%         for j = 1:length(t2)
%             Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
%             error_grid(2,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
%         end
%     end
%     
%         % section based on the arithmetic mean
%     Y_section = zeros(size(Y));
%     Y_anchor_section =  m_arithm(Y);
%     for i_loc = 1:4
%         Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
%         Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
%     end
%     for i = 1:length(t1)
%         for j = 1:length(t2)
%             Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
%             error_grid(3,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
%         end
%     end
%     
%         % section based on the inductive mean
%     Y_section = zeros(size(Y));
%     Y_anchor_section =  m_ind(Y);
%     for i_loc = 1:4
%         Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
%         Y_section(:,:,i_loc) = Y(:,:,i_loc)*Q';
%     end
%     for i = 1:length(t1)
%         for j = 1:length(t2)
%             Ygamma = (1-t1(i))*(1-t2(j))*Y_section(:,:,1) + t1(i)*(1-t2(j))*Y_section(:,:,2) + (1-t1(i))*t2(j)*Y_section(:,:,3) + t1(i)*t2(j)*Y_section(:,:,4);
%             error_grid(4,i,j,i_test) = norm(Ygamma*Ygamma' - C_ref, 'fro')^2;
%         end
%     end
%     
%             
%     for i = 1:4
%         error_grid_rec(i,i_test) = min(min(error_grid(i,:,:,i_test)));
%     end
%         
%     fprintf('------------------------------------------ Errors on grid: %4.2e, %4.2e, %4.2e, %4.2e \n',  error_grid_rec(1,i_test),  error_grid_rec(2,i_test),  error_grid_rec(3,i_test),  error_grid_rec(4,i_test));
    %STOP SANITY CHECK PART------------------------------------------------------------------------------------------------------------------------------------
   
    % implementation of the optimization algorithms
    [~, info] = variable_projection(Y, C_ref);
    sol_SD_rec(:, 1, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(1, i_test) = info.cost(info.nIter);
    fprintf('----------Variable projection finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(1,i_test), sol_SD_rec(1, 1, i_test), sol_SD_rec(2, 1, i_test));
    
    options.foot = 'data';
    [~, info] = SD_patch_section(Y, C_ref, options);
    sol_SD_rec(:, 2, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(2, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach data finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(2,i_test), sol_SD_rec(1, 2, i_test), sol_SD_rec(2, 2, i_test));
    
    options.foot = 'arithm';
    [~, info] = SD_patch_section(Y, C_ref, options);
    sol_SD_rec(:, 3, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(3, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach arithm. finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(3,i_test), sol_SD_rec(1, 3, i_test), sol_SD_rec(2, 3, i_test))
    
    options.foot = 'inductive';
    [~, info] = SD_patch_section(Y, C_ref, options);
    sol_SD_rec(:, 4, i_test) = [info.t1(info.nIter), info.t2(info.nIter)] ;
    error_SD_rec(4, i_test) = info.cost(info.nIter);
    fprintf('----------Section approach ind. finished, optimal cost = %4.2e, sol = %4.2e, %4.2e \n', error_SD_rec(4,i_test), sol_SD_rec(1, 4, i_test), sol_SD_rec(2, 4, i_test));
   
end

save('Results_piecewise_geod.mat', 'error_SD_rec', 'sol_SD_rec', 'params', 'patch', 'r');

    %START SANITY CHECK PART-----------------------------------------------------------------------------------------------------------------------------------
%     save('Results_patchwise_grid.mat', 'error_grid', 'error_grid_rec');
    %STOP SANITY CHECK PART------------------------------------------------------------------------------------------------------------------------------------
