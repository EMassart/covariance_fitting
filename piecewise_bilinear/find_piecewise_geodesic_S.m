function Y_gamma = find_piecewise_geodesic_S(data, sampling, options)
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
        
        % project on the section
        if strcmp(options.input, 'first')
            Y_anchor_section = Y_patch{1};
            for i_loc = 2:4
                Q = orth_pol(Y_anchor_section'*Y_patch{i_loc});
                Y_patch{i_loc} = Y_patch{i_loc}*Q';
            end
        elseif strcmp(options.input, 'arithm')
            Y_anchor_section = m_arithm(Y_patch);
            for i_loc = 1:4
                Q = orth_pol(Y_anchor_section'*Y_patch{i_loc});
                Y_patch{i_loc} = Y_patch{i_loc}*Q';
            end
        elseif strcmp(options.input, 'inductive')
            Y_anchor_section = m_ind(Y_patch);
            for i_loc = 1:4
                Q = orth_pol(Y_anchor_section'*Y_patch{i_loc});
                Y_patch{i_loc} = Y_patch{i_loc}*Q';
            end
        else
            fprintf('Warning: no such choice of the section defined \n');
        end

        % fill in gamma               
        for i_loc = 1:length(t1_out)
            for j_loc = 1:length(t2_out)
                x_out = (i-1)*(length(t1_out)-1)+ i_loc;
                y_out = (j-1)*(length(t2_out)-1)+ j_loc;
                Y_gamma{x_out,y_out} = (1-t1_out(i_loc))*(1-t2_out(j_loc))*Y_patch{1} + t1_out(i_loc)*(1-t2_out(j_loc))*Y_patch{2} + (1-t1_out(i_loc))*t2_out(j_loc)*Y_patch{3} + t1_out(i_loc)*t2_out(j_loc)*Y_patch{4};
            end
        end
        
    end
    
end

end


