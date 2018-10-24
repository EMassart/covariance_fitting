function data = data_points_external_frame(r)

heading = 0:1:4;
magn = 4:3:13;

m = length(magn);
n = length(heading);

data = cell(m,n);
for i = 1:m
    for j = 1:n
        f = strcat('../data_surface_extracted/Y',num2str(heading(j),'%2.1f'),'_',num2str(magn(i),'%2.1f'));
        indx = strfind(f,'.');
        n_loc = length(indx);
        for i_loc = 1:n_loc
            f(indx(n_loc - i_loc + 1)) = '';
        end
        load(f);
        data{i,j} = Y(:,1:r);        
    end
end

end