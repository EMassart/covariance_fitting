function Y_gamma = find_piecewise_geodesic_G(data, sampling)
% Builds a surface based on the data points data, and returns its value at
% the intermediary points.
%
% Author : E.Massart
% Version : July 31, 2018

[m,n] = size(data);
t1_out = linspace(0,1,sampling);
t2_out = linspace(0,1,sampling);

m_out = (length(t1_out)-1)*(m-1)+1;
n_out = (length(t2_out)-1)*(n-1)+1;
Y_gamma = cell(m_out,n_out);

for i = 1:m-1
    for j = 1:n-1
        x_anchors = [i, i+1, i, i+1];
        y_anchors = [j, j, j+1, j+1];
        
        Y_patch = cell(1,4);
        for i_loc = 1:4        
            Y_patch{i_loc} = data{x_anchors(i_loc),y_anchors(i_loc)};
        end

        % fill in gamma               
        for i_loc = 1:length(t1_out)
            for j_loc = 1:length(t2_out)
                x_out = (i-1)*(length(t1_out)-1)+ i_loc;
                y_out = (j-1)*(length(t2_out)-1)+ j_loc;
                Q12 = orth_pol(Y_patch{1}'*Y_patch{2});
                Q34 = orth_pol(Y_patch{3}'*Y_patch{4});
                Y12 = (1-t1_out(i_loc))*Y_patch{1} + t1_out(i_loc)*Y_patch{2}*Q12';
                Y34 = (1-t1_out(i_loc))*Y_patch{3} + t1_out(i_loc)*Y_patch{4}*Q34';
                Q = orth_pol(Y12'*Y34);
                Y_gamma{x_out,y_out} = (1-t2_out(j_loc))*Y12 + t2_out(j_loc)*Y34*Q';
            end
        end
        
    end
    
end

end


