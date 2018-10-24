function Ym = m_ind(Y_data)
% Returns the inductive mean of the set of PSD matrices contained in C
% Author: E. Massart
% Version: October 4, 2018

if iscell(Y_data)
    data = cat(3,Y_data{:});
else
    data = Y_data;
end
s = size(data);

log_riem = @(X,Y) Y*orth_pol(X'*Y)'-X;
exp_riem = @(X,eta) X + eta;

Ym = exp_riem(data(:,:,1),0.5*log_riem(data(:,:,1),data(:,:,2)));
for i = 3:s(3)
    Ym = exp_riem(Ym, log_riem(Ym,data(:,:,i))./i);
end

end