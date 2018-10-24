function [t2,Q,H] = t2_opt_Var_Proj(Y12, Y34, C_hat)
% Computes t_2^*(t_1) for the variable projection method used to compute the closest point from a patch to a point C_hat.
% The patch is built based on geodesics from the quotient geometry.
% The computation of t_2^*(t_1) requires to solve a cubic equation S1 t_2^3 +
% S2 t_2^2 + S3 t_2 + S4, where the coefficients S1, S2, S3, S4 depend on
% t_1.

% Author : E. Massart
% Last modified on October 24, 2018

[Q, H] = orth_pol(Y12'*Y34);
A = Y12*Y12';
B = Y34*Y34';
C = Y12*Q*Y34' + Y34*Q'*Y12';      
tilde_A = A+B-C;
tilde_B = C-2*A;
tilde_C = A - C_hat;

S1 = 2*sum(sum(tilde_A.^2));
S2 = 3*sum(sum(tilde_A.*tilde_B));
S3 = 2*sum(sum(tilde_A.*tilde_C)) + sum(sum(tilde_B.^2));
S4 = sum(sum(tilde_B.*tilde_C));
sol = roots([S1, S2, S3, S4]);

% The output of roots will be of length three, since it also contains the
% complex roots. We remove the complex roots from the solution, and choose,
% in the case for which there are several real roots, the one leading to
% the smallest value of the cost function f(t_1,t_2) = || \gamma(t_1,t_2) -
% C_hat||_F^2.

to_take = zeros(1,length(sol));
for i = 1:length(sol)
    to_take(i) = (abs(imag(sol(i))) <= 1e-6);
end

s = sum(to_take);

if s < 1
    fprintf('Less than one real solution to the cubic equation \n sol = %14.8e, %14.8e, %14.8e \n', sol(1), sol(2), sol(3));
    fprintf('Code is stopping here..... \n');
elseif s > 1
    fprintf('More than one real solution to the cubic equation \n sol = %14.8e, %14.8e, %14.8e \n', sol(1), sol(2), sol(3));
    fprintf('Choose the one leading to smallest cost function... \n');
    indx = (to_take==1);
    t2_cur = zeros(1,s);
    cost_cur = zeros(1,s);
    for i_cur = 1:s
        t2_cur(i_cur) = real(sol(indx(i_cur)));
        Y_gamma = (1-t2_cur(i_cur))*Y12 + t2_cur(i_cur)*Y34*Q';
        gamma = Y_gamma*Y_gamma';
        cost_cur(i_cur) = norm(C_hat - gamma,'fro');
    end
    [~,indx_opt] = min(cost_cur);
    t2 = t2_cur(indx_opt);
else
    t2 = real(sol(to_take==1));
end




end