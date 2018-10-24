function [x_end, info] = steepest_descent_for_bezier_quotient(train_data, C_ref, options)

% Parameters for Bézier
manifold = 'psd_quotient';
global_variables
geo_functions
problem  = prepare_structure(train_data,manifold,11);
problem = control_points_simple_generation_2d(problem);
problem.C = C_ref;
[n_dim,m_dim] = size(train_data);

% Params finite differences
tol = 1e-6;
tolFD = 1e-6;
stepFD = 0.1;
factorFD = 0.5;
maxiterFD = 20;

% Params Armijo    
alpha = 0.01;
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
    
    % estimation for dfdt1
    kFD = 1;
    gapFD = zeros(1, maxiterFD); gapFD(1) = 1;
    stepFD_cur = stepFD;  % difference forward
    if x_end(1,k)+factorFD^kFD*stepFD_cur > n_dim-1 
        stepFD_cur = - stepFD;  % difference backward
    end
    dfdx_prev = (cost_function_bezier(problem, x_end(1,k)+factorFD^kFD*stepFD_cur, x_end(2,k)) - f(k)) / (factorFD^kFD*stepFD_cur);
    while gapFD(kFD) > tolFD && kFD < maxiterFD
        kFD = kFD+1;
        dfdx = (cost_function_bezier(problem, x_end(1,k)+factorFD^kFD*stepFD_cur, x_end(2,k)) - f(k)) / (factorFD^kFD*stepFD_cur);
        gapFD(kFD) = abs(dfdx - dfdx_prev);
        dfdx_prev = dfdx;
    end
    gapFD = gapFD(2:end);
%     figure;
%     plot(1:length(gapFD), gapFD);
    
    % estimation for dfdt2
    kFD = 1;
    gapFD = zeros(1, maxiterFD); gapFD(1) = 1;
    stepFD_cur = stepFD;
    if x_end(2,k)+factorFD^kFD*stepFD_cur > m_dim-1  % difference forward
        stepFD_cur = - stepFD;
    end
    dfdy_prev = (cost_function_bezier(problem, x_end(1,k), x_end(2,k)+factorFD^kFD*stepFD_cur) - f(k)) / (factorFD^kFD*stepFD_cur);
    while gapFD(kFD) > tolFD && kFD < maxiterFD
        kFD = kFD+1;
        dfdy = (cost_function_bezier(problem, x_end(1,k), x_end(2,k)+factorFD^kFD*stepFD_cur) - f(k)) / (factorFD^kFD*stepFD_cur);
        gapFD(kFD) = abs(dfdy - dfdy_prev);
        dfdy_prev = dfdy;
    end
    gapFD = gapFD(2:end);
%     figure;
%     plot(1:length(gapFD), gapFD);
    
    grad(:,k) = [dfdx; dfdy];

    fprintf('k = %d, f = %4.2e, grad = %4.2e, x = %4.2e, %4.2e \n', k, f(k), norm(grad(:,k)), x_end(1,k), x_end(2,k));

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