function Ym = m_arithm(Y_data)
% Returns the arithmetic mean of the set of PSD matrices contained in C
% Author: E. Massart
% Version: October 4, 2018
if iscell(Y_data)
    data = cat(3,Y_data{:});
else
    data = Y_data;
end

s = size(data);
M = zeros(s(1), s(1));
for i = 1:s(3)
    M = M + data(:,:,i)*data(:,:,i)';
end
M = M./s(3);

[U, D] = eig(M);
[d,indx] = sort(diag(D),'descend');
Ym = U(:,indx)*diag(sqrt(d));
Ym = Ym(:,1:s(2));

end