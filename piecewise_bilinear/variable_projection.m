function [x_sol, info] = variable_projection(Y, C_hat, options)

% Apply a variable projection method to the computation of the point of a
% surface that is the closest to a given data point C_hat.
%
% Inputs : Y : 3D array, that contains the PSD matrices Y_i, along the
%              third dimension (i.e., Y is of size n x r x 4). Those 
%              matrices will be the anchor points of the surface.
%          C_hat : PSD matrix obtained exprimentally. The goal is to find
%                  the point of the surface that is the closest to C_hat
%
% Output : x_sol : vector x_sol = [t1_opt, t2_opt], such that the point of the surface 
%                  that is the closest from C_hat lies on the surface at
%                  parameters (t1_opt, t2_opt).
%          info : additional information, mainly related to the convergence of the method
% 
% Author : E.Massart
% Last modification: October 24, 2018                

% ------------------------------------------------------------------------------
% ---------------------------------   INITIALIZATION ---------------
% ------------------------------------------------------------------------------


if  ~exist('options','var')
    x0 = 0.5;
else 
    x0 = options.x0;
end

% -------------------------------  Other parameters
maxiter = 1000;                % maximal number of iterations allowed
t1 = zeros(1,maxiter+1);        
t2 = zeros(1,maxiter+1);
grad = zeros(1,maxiter);
cost = zeros(1,maxiter+1);
stop = 0;
tol = 1e-12;                     % tolerance on the stopping criterion                 
k = 1;                          % iteration counter
t1(1) = x0(1);                    % initialization for t1. The variable t2 will be initialized at t2^*(t_1), 
                                % i.e., the optimal value of t_2, for that value of t_1.
% Params Armijo    
alpha = 0.005;
beta = 0.5;
m_max = 50;
factor = 1;

% ------------------------------- Computes the current location on the surface
orth12 = orth_pol(Y(:,:,1)'*Y(:,:,2));
orth34 = orth_pol(Y(:,:,3)'*Y(:,:,4));

Y12 = (1-t1(k))*Y(:,:,1) + t1(k)*Y(:,:,2)*orth12';
Y34 = (1-t1(k))*Y(:,:,3) + t1(k)*Y(:,:,4)*orth34';

[t2(1),Q,H] = t2_opt_Var_Proj(Y12, Y34, C_hat);
Y_gamma = (1-t2(1))*Y12 + t2(1)*Y34*Q';
gamma = Y_gamma*Y_gamma';
cost(1) = norm(C_hat - gamma,'fro')^2 ;

% ---------------------------------------------------------------------------------------------
% ---------------------------------   ITERATIONS OF GRADIENT DESCENT START HERE ---------------
% ---------------------------------------------------------------------------------------------

while ~stop 

    if mod(k,10)==0, fprintf('--------------------------------------------------------------------  Iteration number %d \n',k); end

    % -----------------------------------------------------  Computes the gradient
    dot_Y12 = -Y(:,:,1) + Y(:,:,2)*orth12';
    dot_Y34 = -Y(:,:,3) + Y(:,:,4)*orth34';
    
    % derivative of orthogonal factor Q
    dotM = dot_Y12'*Y34 + Y12'*dot_Y34;
    Omega = sylvester(H,H,dotM*Q' - Q*dotM');
    dot_Q = Omega*Q;
    
    dot_Y_gamma = (1-t2(k))*dot_Y12 + t2(k)*dot_Y34*Q' + t2(k)*Y34*dot_Q';
        
    % final derivatives
    dot_gamma = dot_Y_gamma*Y_gamma' + Y_gamma*dot_Y_gamma';  
    grad(k) = 2*trace(dot_gamma*(gamma - C_hat)');
        
    % ----------------------------------------------------  Computes the new iterate + evaluates the cost function
    t1(k+1) = t1(k) - alpha*factor^k*grad(k);
    
    Y12 = (1-t1(k+1))*Y(:,:,1) + t1(k+1)*Y(:,:,2)*orth12';
    Y34 = (1-t1(k+1))*Y(:,:,3) + t1(k+1)*Y(:,:,4)*orth34';
    [t2(k+1),Q,H] = t2_opt_Var_Proj(Y12, Y34, C_hat);
    Y_gamma = (1-t2(k+1))*Y12 + t2(k+1)*Y34*Q';
    gamma = Y_gamma*Y_gamma';
    cost(k+1) = norm(C_hat - gamma,'fro')^2;      
    fprintf('Cost = %4.2e, t1 = %4.2e, t2 = %4.2e, grad = %4.2e \n', cost(k), t1(k+1), t2(k+1), grad(k));
    
    if (cost(k) - cost(k+1)) < 0, disp('Cost function has increased -> decrease stepsize...'); end
    m = 0;
    while (cost(k) - cost(k+1)) < 0 && m < m_max
        m = m+1;
        t1(k+1) = t1(k) - alpha*factor^k*beta^m*grad(k);
        Y12 = (1-t1(k+1))*Y(:,:,1) + t1(k+1)*Y(:,:,2)*orth12';
        Y34 = (1-t1(k+1))*Y(:,:,3) + t1(k+1)*Y(:,:,4)*orth34';
        [t2(k+1),Q,~] = t2_opt_Var_Proj(Y12, Y34, C_hat);
        Y_gamma = (1-t2(k+1))*Y12 + t2(k+1)*Y34*Q';
        gamma = Y_gamma*Y_gamma';
        cost(k+1) = norm(C_hat - gamma,'fro')^2;
    end
    
    fprintf('Step after decreasing : %4.2e \n', alpha*factor^k*beta^m*grad(k));
    if m == m_max, disp('However, this does not allows to decrese the cost function'); end


    % ----------------------------------------------------  Evaluates the stopping criterion   
    if k == maxiter, stop = 1; disp('STOPPING : Maximum number of iterations reached'); end
    if abs(cost(k+1)-cost(k)) < tol
        stop = 1;
        fprintf('STOPPING : decreasing of the cost function is too small: %4.2e \n', abs(cost(k+1)-cost(k)));
    end
        
    k = k+1;
    
end

info = struct();
info.t1 = t1;
info.t2 = t2;
info.grad = grad;
info.cost = cost;
info.nIter = k;

x_sol = [t1(k), t2(k)];

end

   