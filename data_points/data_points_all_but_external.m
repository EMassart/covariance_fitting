function [data, params,patch] = data_points_all_but_external(r)
% This code loads the data, complementary to the code data_points_external_frame.m
% (the test set that we used).
% The input r is the estimation of the rank of the covariance matrices (we
% are going to work on the manifold of PSD matrices of rank r)
% Author : E. Massart
% Last modification: October 24, 2018

% Values of the parameters that we will consider in our test set
magn_type1 = 4:1.5:13;
heading_type1 = 0.5:1:3.5;

magn_type2 = 5.5:3:11.5;
heading_type2 = 0:1:4;


% Load the data, and truncate the Y matrices to r columns
m1 = length(magn_type1);
n1 = length(heading_type1);

m2 = length(magn_type2);
n2 = length(heading_type2);

l_tot = m1*n1 + m2*n2;

data = cell(1,l_tot);
params = cell(1,l_tot);

for i = 1:l_tot
    if i <= m1*n1
        i_heading = floor(i/m1)+1;     % indx of the heading
        if (i-floor(i/m1)*m1 < 1)
           i_heading = floor(i/m1);
        end
        i_magn = i - (i_heading-1)*m1;
        
        f = strcat('../data_surface_extracted/Y',num2str(heading_type1(i_heading),'%2.1f'),'_',num2str(magn_type1(i_magn),'%04.1f'));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i} = Y(:,1:r); 
        params{i} = [magn_type1(i_magn),heading_type1(i_heading)];
    else
        i = i-m1*n1;
        i_heading = floor(i/m2)+1;     % indx of the heading
        if (i-floor(i/m2)*m2 < 1)
            i_heading = floor(i/m2);
        end
        i_magn = i - (i_heading-1)*m2;
        
        f = strcat('../data_surface_extracted/Y',num2str(heading_type2(i_heading),'%2.1f'),'_',num2str(magn_type2(i_magn),'%04.1f'));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i+m1*n1} = Y(:,1:r); 
        params{i+m1*n1} = [magn_type2(i_magn),heading_type2(i_heading)]; 
    end  
end


% If required, the patch to which each test point belongs 
if nargout == 3
    patch = cell(1,43);
%     patch(1:3) = {[1,1]};
%     patch(4:5) = {[1,2]};
%     patch(6:7) = {[1,3]};
%     patch(8:10) = {[2,1]};
%     patch(11:12) = {[2,2]};
%     patch(13:14) = {[2,3]};
%     patch(15:17) = {[3,1]};
%     patch(18:19) = {[3,2]};
%     patch(20:21) = {[3,3]};
%     patch(22:24) = {[4,1]};
%     patch(25:26) = {[4,2]};
%     patch(27:28) = {[4,3]};
%     patch(29) = {[1,1]};
%     patch(30) = {[1,2]};
%     patch(31) = {[1,3]};
%     patch(32) = {[1,1]};
%     patch(33) = {[1,2]};
%     patch(34) = {[1,3]};
%     patch(35) = {[2,1]};
%     patch(36) = {[2,2]};
%     patch(37) = {[2,3]};
%     patch(38) = {[3,1]};
%     patch(39) = {[3,2]};
%     patch(40) = {[3,3]};
%     patch(41) = {[4,1]};
%     patch(42) = {[4,2]};
%     patch(43) = {[4,3]};


%     patch(1:2) = {[1,1]};
%     patch(3:4) = {[1,2]};
%     patch(5:7) = {[1,3]};
%     patch(8:9) = {[2,1]};
%     patch(10:11) = {[2,2]};
%     patch(12:14) = {[2,3]};
%     patch(15:16) = {[3,1]};
%     patch(17:18) = {[3,2]};
%     patch(19:21) = {[3,3]};
%     patch(22:23) = {[4,1]};
%     patch(24:25) = {[4,2]};
%     patch(26:28) = {[4,3]};
%     patch(29) = {[1,1]};
%     patch(30) = {[1,2]};
%     patch(31) = {[1,3]};
%     patch(32) = {[2,1]};
%     patch(33) = {[2,2]};
%     patch(34) = {[2,3]};
%     patch(35) = {[3,1]};
%     patch(36) = {[3,2]};
%     patch(37) = {[3,3]};
%     patch(38) = {[4,1]};
%     patch(39) = {[4,2]};
%     patch(40) = {[4,3]};
%     patch(41) = {[4,1]};
%     patch(42) = {[4,2]};
%     patch(43) = {[4,3]};

    patch(1:2) = {[1,1]};
    patch(3:4) = {[2,1]};
    patch(5:7) = {[3,1]};
    patch(8:9) = {[1,2]};
    patch(10:11) = {[2,2]};
    patch(12:14) = {[3,2]};
    patch(15:16) = {[1,3]};
    patch(17:18) = {[2,3]};
    patch(19:21) = {[3,3]};
    patch(22:23) = {[1,4]};
    patch(24:25) = {[2,4]};
    patch(26:28) = {[3,4]};
    patch(29) = {[1,1]};
    patch(30) = {[2,1]};
    patch(31) = {[3,1]};
    patch(32) = {[1,2]};
    patch(33) = {[2,2]};
    patch(34) = {[3,2]};
    patch(35) = {[1,3]};
    patch(36) = {[2,3]};
    patch(37) = {[3,3]};
    patch(38) = {[1,4]};
    patch(39) = {[2,4]};
    patch(40) = {[3,4]};
    patch(41) = {[1,4]};
    patch(42) = {[2,4]};
    patch(43) = {[3,4]};
end
    
    
end