function [data, params] = data_points_all_but_external(r)

heading_type1 = 0.5:1:3.5;
magn_type1 = 4:1.5:13;

heading_type2 = 0:1:4;
magn_type2 = 5.5:3:11.5;

n1 = length(heading_type1);
m1 = length(magn_type1);

n2 = length(heading_type2);
m2 = length(magn_type2);

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
        
        f = strcat('../data_surface_extracted/Y',num2str(heading_type1(i_heading),'%2.1f'),'_',num2str(magn_type1(i_magn),'%2.1f'));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i} = Y(:,1:r); 
        params{i} = [heading_type1(i_heading), magn_type1(i_magn)];
    else
        i = i-m1*n1;
        i_heading = floor(i/m2)+1;     % indx of the heading
        if (i-floor(i/m2)*m2 < 1)
            i_heading = floor(i/m2);
        end
        i_magn = i - (i_heading-1)*m2;
        
        f = strcat('../data_surface_extracted/Y',num2str(heading_type2(i_heading),'%2.1f'),'_',num2str(magn_type2(i_magn),'%2.1f'));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i+m1*n1} = Y(:,1:r); 
        params{i+m1*n1} = [heading_type2(i_heading), magn_type2(i_magn)]; 
    end  
end

end