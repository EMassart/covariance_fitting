function [x_end, info] = steepest_descent_for_bezierS(train_data, C_ref, options)

% Parameters for Bézier
[n_dim,m_dim] = size(train_data);
manifold = 'euclidean';
global_variables
geo_functions

% projection of the train_data on the section
if strcmp(options.input, 'first')
    Y_anchor_section = train_data{1,1};
    for i = 1:n_dim
        for j = 1:m_dim
            if (i~=1 || j~=1)
                Q = orth_pol(Y_anchor_section'*train_data{i,j});
                train_data{i,j} = train_data{i,j}*Q';
            end
        end
    end
elseif strcmp(options.input, 'arithm')
    Y_anchor_section = m_arithm(train_data);
    for i = 1:n_dim
        for j = 1:m_dim
            Q = orth_pol(Y_anchor_section'*train_data{i,j});
            train_data{i,j} = train_data{i,j}*Q';
        end
    end
elseif strcmp(options.input, 'inductive')
    Y_anchor_section = m_ind(train_data);
    for i = 1:n_dim
        for j = 1:m_dim
            Q = orth_pol(Y_anchor_section'*train_data{i,j});
            train_data{i,j} = train_data{i,j}*Q';
        end
    end
else
    fprintf('Warning: no such choice of the section defined \n');
end

problem  = prepare_structure(train_data,manifold,11);
problem = control_points_simple_generation_2d(problem);
problem.C = C_ref;

% Params finite differences
tol = 1e-6;
% tolFD = 1e-6;
% stepFD = 0.1;
% factorFD = 0.5;
% maxiterFD = 20;

% Params Armijo    
alpha = 1;
beta = 0.8;
sigma = 0.5;
m_max = 75;

% other initializations
stop = 0;
maxiter = 1000;
x_end = zeros(2,maxiter+1);
f = zeros(1, maxiter+1);
grad = zeros(2, maxiter);
k = 1;
x_end(:,k) = [0.5, 0.5]';


while ~stop;
    
    f(k) = cost_function_bezier(problem, x_end(1,k), x_end(2,k));
    grad(:,k) = grad_analytic(problem, x_end(1,k), x_end(2,k));
    
%     % estimation for dfdt1
%     kFD = 1;
%     gapFD = zeros(1, maxiterFD); gapFD(1) = 1;
%     stepFD_cur = stepFD;  % difference forward
%     if x_end(1,k)+factorFD^kFD*stepFD_cur > n_dim-1 
%         stepFD_cur = - stepFD;  % difference backward
%     end
%     dfdx = zeros(1,maxiterFD);  dfdx(1) = (cost_function_bezier(problem, x_end(1,k)+factorFD^kFD*stepFD_cur, x_end(2,k)) - f(k)) /(factorFD^kFD*stepFD_cur);
%     while gapFD(kFD) > tolFD && kFD < maxiterFD
%         kFD = kFD+1;
%         dfdx(kFD) = (cost_function_bezier(problem, x_end(1,k)+factorFD^kFD*stepFD_cur, x_end(2,k)) - f(k)) / (factorFD^kFD*stepFD_cur);
%         gapFD(kFD) = abs(dfdx(kFD) - dfdx(kFD-1));
%     end
%     gapFD = gapFD(2:end);  dfdx = dfdx(1:kFD);
%     figure;
% %     plot(1:kFD, dfdx, '.-b', 1:kFD, grad_ana(1)*ones(1,kFD),'-r');
%     plot(1:kFD, dfdx, '.-b');
% 
%     % estimation for dfdt2
%     kFD = 1;
%     gapFD = zeros(1, maxiterFD); gapFD(1) = 1;
%     stepFD_cur = stepFD;
%     if x_end(2,k)+factorFD^kFD*stepFD_cur > m_dim-1  % difference forward
%         stepFD_cur = - stepFD;
%     end
%     dfdy = zeros(1,maxiterFD);  dfdy(1) = (cost_function_bezier(problem, x_end(1,k), x_end(2,k)+factorFD^kFD*stepFD_cur) - f(k)) / (factorFD^kFD*stepFD_cur);
%     while gapFD(kFD) > tolFD && kFD < maxiterFD
%         kFD = kFD+1;
%         dfdy(kFD) = (cost_function_bezier(problem, x_end(1,k), x_end(2,k)+factorFD^kFD*stepFD_cur) - f(k)) / (factorFD^kFD*stepFD_cur);
%         gapFD(kFD) = abs(dfdy(kFD) - dfdy(kFD-1));
%     end
%     gapFD = gapFD(2:end); dfdy = dfdy(1:kFD);
%     figure;
% %     plot(1:kFD, dfdy, '.-b', 1:kFD, grad_ana(2)*ones(1,kFD),'-r');
%     plot(1:kFD, dfdy, '.-b');

%     grad(:,k) = [dfdx(end); dfdy(end)];

    fprintf('k = %d, f = %4.2e, grad = %4.2e,  %4.2e, norm_grad = %4.2e, x = %4.2e, %4.2e \n', k, f(k), grad(1,k), grad(2,k), norm(grad(:,k)), x_end(1,k), x_end(2,k));

    % perform the update using Armijo backtracking
    gapA = -5;
    m = 0;
    while gapA <= 0 && m < m_max
        x_new = x_end(:,k) - alpha*beta^m*grad(:,k);
        
        % check whether we are still in the domain, project otherwise
        if x_new(1) < 0
            x_new(1)  = 0;
        end
        if x_new(2) < 0
            x_new(2)  = 0;
        end
        if x_new(1) >= n_dim-1
            x_new(1)  = n_dim-1;
        end
        if x_new(2) >= m_dim-1
            x_new(2)  = m_dim-1;
        end
        
        stepA = x_end(:,k) - x_new;
        gapA = f(k) - cost_function_bezier(problem, x_new(1), x_new(2)) - sigma*grad(:,k)'*stepA;
        m = m+1;
    end
    

    if m == m_max disp('Maxiter number of iterations Armjo reached'); end
    
    if k >= 3
        stop = abs(f(k-1) - f(k-2)) <= tol;
        fprintf('Stopping criterion = %4.2e, m = %d \n',  abs(f(k-1) - f(k-2)), m);
    end
    
    if k == maxiter
        stop = 1;
    end
    
    k = k+1;
    x_end(:,k) = x_new;
end

f(k) = cost_function_bezier(problem, x_end(1,k), x_end(2,k));
info.f = f;
info.x_end = x_end;
info.k = k;
info.grad = grad;
x_end = x_end(:,k);

end