function [x_sol, info] = sectional_approach_final(Y, C_hat, options)

% Apply a steepest descent to the computation of the point of a sectional
% surface that is the closest to a given data point C_hat.
%
% Inputs : A : 3D array, that contains the PSD matrices A_i, along the
%              third dimension (i.e., A is of size n x n x 4, where n is the size of the
%              PSD matrices to consider). Those matrices will be the anchor
%              points of the surface.
%          r : Rank of the PSD matrices A_i
%          C_hat : PSD matrix obtained exeprimentally. The goal is to find
%                  the point of the surface that is the closest to C_hat
%
% Output : x_sol : vector x_sol = [t1_opt, t2_opt], such that the point of the surface 
%                  that is the closest from C_hat lies on the surface at
%                  parameters (t1_opt, t2_opt).
%          info : additional information, mainly related to the convergence of the method
% 
% Author : E.Massart
% Version : March 13, 2018
                

% ------------------------------------------------------------------------------
% ---------------------------------   INITIALIZATION ---------------
% ------------------------------------------------------------------------------

if ~isfield(options, 'x0')
    x0 = [0.5, 0.5];
else
    x0 = options.x0;
end
if ~isfield(options, 'input')
    options.input = 'inductive';
end

% -------------------------------  Other parameters
maxiter = 1000;                % maximal number of iterations allowed
t1 = zeros(1,maxiter+1);        
t2 = zeros(1,maxiter+1);
grad = zeros(2,maxiter);
cost = zeros(1,maxiter+1);
stop = 0;
tol = 1e-10;            
k = 1;                          % iteration counter
t1(1) = x0(1);                    % initialization for t1.
t2(1) = x0(2);                    % initialization for t2.
m_max = 75;
beta = 0.8;
sigma = 0.5;
alpha = 1;

% project on the section
Y_gamma = zeros(size(Y));
if strcmp(options.input, 'first')
    Y_anchor_section = Y(:,:,1);
    Y_gamma(:,:,1) =  Y(:,:,1);
    for i_loc = 2:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_gamma(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
elseif strcmp(options.input, 'arithm')
    Y_anchor_section = m_arithm(Y);
    for i_loc = 1:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_gamma(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
elseif strcmp(options.input, 'inductive')
    Y_anchor_section = m_ind(Y);
    for i_loc = 1:4
        Q = orth_pol(Y_anchor_section'*Y(:,:,i_loc));
        Y_gamma(:,:,i_loc) = Y(:,:,i_loc)*Q';
    end
else
    fprintf('Warning: no such choice of the section defined \n');
end

Y = Y_gamma;

% We are able to compute Y_gamma
Y_gamma = (1-t1(1))*(1-t2(1))*Y(:,:,1) + t1(1)*(1-t2(1))*Y(:,:,2) + (1-t1(1))*t2(1)*Y(:,:,3) + t1(1)*t2(1)*Y(:,:,4);
gamma = Y_gamma*Y_gamma';

% evaluate cost function for those values of t1 and many values of t2
cost(1) = norm(gamma - C_hat,'fro')^2;


% ---------------------------------------------------------------------------------------------
% ---------------------------------   ITERATIONS OF GRADIENT DESCENT START HERE ---------------
% ---------------------------------------------------------------------------------------------

while ~stop 

    if mod(k,10)==0, fprintf('--------------------------------------------------------------------  Iteration number %d \n',k); end

    % ---------------------------------  Computes the gradient
    deltaYGamma_t1 = t2(k)*(Y(:,:,1) - Y(:,:,2) - Y(:,:,3) + Y(:,:,4)) + (Y(:,:,2) - Y(:,:,1));
    deltaYGamma_t2 = t1(k)*(Y(:,:,1) - Y(:,:,2) - Y(:,:,3) + Y(:,:,4)) + (Y(:,:,3) - Y(:,:,1));
    
    % final derivatives
    
    delta_gamma_t1 = deltaYGamma_t1 * Y_gamma' + Y_gamma * deltaYGamma_t1';
    delta_gamma_t2 = deltaYGamma_t2 * Y_gamma' + Y_gamma * deltaYGamma_t2';

    grad(:,k) = [2*trace(delta_gamma_t1*(gamma - C_hat)'); 2*trace(delta_gamma_t2*(gamma - C_hat)')];
 
    
    %start Armijo here
    gapA = -5;
    m = 0;
    while gapA <= 0 && m < m_max
        t1(k+1) = t1(k) - alpha*beta^m*grad(1,k);
        t2(k+1) = t2(k) - alpha*beta^m*grad(2,k);
        
        % check whether we are still in the domain, project otherwise
        if t1(k+1) < 0
            t1(k+1)  = 0;
        end
        if t2(k+1) < 0
            t2(k+1)  = 0;
        end
        if t1(k+1) >= 1
            t1(k+1)  = 1;
        end
        if t2(k+1) >= 1
            t2(k+1)  = 1;
        end
        
        Y_gamma = (1-t1(k+1))*(1-t2(k+1))*Y(:,:,1) + t1(k+1)*(1-t2(k+1))*Y(:,:,2) + (1-t1(k+1))*t2(k+1)*Y(:,:,3) + t1(k+1)*t2(k+1)*Y(:,:,4);
        gamma = Y_gamma*Y_gamma';
        cost(k+1) = norm(C_hat - gamma,'fro')^2;
    
        stepA = [t1(k) - t1(k+1), t2(k) - t2(k+1)]';
        gapA = cost(k) - cost(k+1) - sigma*grad(:,k)'*stepA;
        m = m+1;
    end
    
    fprintf('Cost = %4.2e, t1 = %4.2e, t2 = %4.2e, norm grad = %4.2e \n', cost(k+1), t1(k+1), t2(k+1), norm(grad(:,k)));
    fprintf('Stepsize after decreasing :  %4.2e \n', alpha*beta^m);
    if m == m_max, disp('However, this does not allows to decrese the cost function'); end
    
    % ----------------------------------------------------  Evaluates the stopping criterion   
    if k == maxiter, stop = 1; disp('STOPPING : Maximum number of iterations reached'); end
    fprintf('Stopping criterion: %4.2e \n', abs(cost(k+1)-cost(k)));
    if abs(cost(k+1)-cost(k))< tol
        stop = 1;
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


    % -------------------------------  Implement Armijo conditions for stepsize
%     m = 0;
%     m_max = 150;
%     sigma = 0.0001;
%     beta = 0.5;
%     alpha = 0.0001;
%     t1(k+1) = t1(k) - alpha*beta^m*grad(1,k);
%     t2(k+1) = t2(k) - alpha*beta^m*grad(2,k);
%     Y_gamma = (1-t1(k+1))*(1-t2(k+1))*Y(:,:,1) + t1(k+1)*(1-t2(k+1))*Y(:,:,2) + (1-t1(k+1))*t2(k+1)*Y(:,:,3) + t1(k+1)*t2(k+1)*Y(:,:,4);
%     gamma = Y_gamma*Y_gamma';
%     cost(k+1) = norm(C_hat - gamma,'fro')^2;
%     while (cost(k) - cost(k+1)) < sigma*alpha*beta^m*norm(grad(:,k))^2 && m < m_max
%         m = m+1;
%         t1(k+1) = t1(k) - alpha*beta^m*grad(1,k);
%         t2(k+1) = t2(k) - alpha*beta^m*grad(2,k);
%         Y_gamma = (1-t1(k+1))*(1-t2(k+1))*Y(:,:,1) + t1(k+1)*(1-t2(k+1))*Y(:,:,2) + (1-t1(k+1))*t2(k+1)*Y(:,:,3) + t1(k+1)*t2(k+1)*Y(:,:,4);
%         gamma = Y_gamma*Y_gamma';
%         cost(k+1) = norm(C_hat - gamma,'fro')^2;
%     end
%     
%     if m == m_max,
%         fprintf('Maximum m allowed in armijo backtracking \n');
%         m = 0;
%         m_max = 150;
%         beta = 0.5;
%         alpha = 1;
%         t1(k+1) = t1(k) - alpha*beta^m*grad(1,k);
%         t2(k+1) = t2(k) - alpha*beta^m*grad(2,k);
%         Y_gamma = (1-t1(k+1))*(1-t2(k+1))*Y(:,:,1) + t1(k+1)*(1-t2(k+1))*Y(:,:,2) + (1-t1(k+1))*t2(k+1)*Y(:,:,3) + t1(k+1)*t2(k+1)*Y(:,:,4);
%         gamma = Y_gamma*Y_gamma';
%         cost(k+1) = norm(C_hat - gamma,'fro')^2;
%         while (cost(k) - cost(k+1)) < 0 && m < m_max
%             m = m+1;
%             t1(k+1) = t1(k) - alpha*beta^m*grad(1,k);
%             t2(k+1) = t2(k) - alpha*beta^m*grad(2,k);
%             Y_gamma = (1-t1(k+1))*(1-t2(k+1))*Y(:,:,1) + t1(k+1)*(1-t2(k+1))*Y(:,:,2) + (1-t1(k+1))*t2(k+1)*Y(:,:,3) + t1(k+1)*t2(k+1)*Y(:,:,4);
%             gamma = Y_gamma*Y_gamma';
%             cost(k+1) = norm(C_hat - gamma,'fro')^2;
%         end 
%     end

