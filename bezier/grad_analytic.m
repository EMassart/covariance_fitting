function grad = grad_analytic(pb, t1, t2)

% Input : -  control : cell of 4 x 4 control points

% Compute the associated control points
b      	= pb.control;
d       = pb.d;
n 		= pb.n;
m 		= pb.m;

k = floor(t1)+1;
l = floor(t2)+1;

if t1 == n-1,  k = t1; end
if t2 == m-1,  l = t2; end

ind_x = ((k-1)*d + 1): (k*d +1);
ind_y = ((l-1)*d + 1): (l*d +1);

bb = b( ind_x ,  ind_y );

t1 = t1-k+1;
t2 = t2-l+1;

%define the Bernstein polynomial and its derivative
bernstein_poly = @(j, t) (factorial(d)/(factorial(j)*factorial(d-j)))*t^j*(1-t)^(d-j);
bernstein_deriv = @b_deriv;

function dB = b_deriv(j,t)

    coeff = factorial(d)/(factorial(j)*factorial(d-j));
    if j == 0
        dB = -coeff*d*(1-t)^(d-1);
    elseif j == d
        dB = coeff*d*t^(d-1);
    else
        dB = coeff*(j*t^(j-1)*(1-t)^(d-j) - t^j*(d-j)*(1-t)^(d-j-1));
    end

end

s = size(b{1,1});
dbdt1 = zeros(s);
dbdt2 = zeros(s);
b = zeros(s);

% computation of the gradient
for i = 1:d+1
    for j = 1:d+1
        dbdt1 = dbdt1 + bb{i,j}*bernstein_deriv(i-1,t1)*bernstein_poly(j-1,t2);
        dbdt2 = dbdt2 + bb{i,j}*bernstein_poly(i-1,t1)*bernstein_deriv(j-1,t2);
        b = b + bb{i,j}*bernstein_poly(i-1,t1)*bernstein_poly(j-1,t2);
    end
end

dfdt1 = -2*trace((pb.C - b*b')'*(dbdt1*b' + b*dbdt1'));
dfdt2 = -2*trace((pb.C - b*b')'*(dbdt2*b' + b*dbdt2'));
grad = [dfdt1 dfdt2]; 

end