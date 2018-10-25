function [data,params] = data_points_external_frame(r)
% This code loads the data corresponding to the external frame (the
% training set that we used).
% The input r is the estimation of the rank of the covariance matrices (we
% are going to work on the manifold of PSD matrices of rank r)
% Author : E. Massart
% Last modification: October 24, 2018

% Values of the parameters that will correspond to the external frame
heading = 0:1:4;
magn = 4:3:13;

% Load the data, and truncate the Y matrices to r columns. 
n = length(heading);
m = length(magn);

data = cell(m,n);
l_tot = m*n;
params = cell(1,l_tot);

count = 0;
for i = 1:m
    for j = 1:n
        count = count +1;
        f = strcat('../data_surface_extracted/Y',sprintf('%2.1f',heading(j)),'_',sprintf('%04.1f',magn(i)));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i,j} = Y(:,1:r);        
        params{count} = [heading(j), magn(i)];
    end
end

end