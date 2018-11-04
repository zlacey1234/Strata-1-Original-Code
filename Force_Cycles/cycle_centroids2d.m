function out = cycle_centroids2d(c_list,coords)

num_cycle = length(c_list(:,1));
c_length = length(c_list(1,:));
xlength = range(coords(:,1));
ylength = range(coords(:,2));

centroid_list = zeros(num_cycle, 2);

for i = 1:num_cycle
    r_pos = 0;
    r_list = zeros(c_length,length(coords(1,:)));
    for j = 1:c_length
        r_pos = r_pos + 1;
        r_list(r_pos,:) = coords(c_list(i,j),1:2);
    end
    % Adjustment for periodic boundary conditions; maybe replace with
    % actual dimensions of the box
    if range(r_list(:,1)) > 0.75*xlength
        r_list(r_list(:,1)<0,1) = r_list(r_list(:,1)<0,1) + range(r_list(:,1));
    end
    if range(r_list(:,2)) > 0.75*ylength
        r_list(r_list(:,2)<0,2) = r_list(r_list(:,2)<0,2) + range(r_list(:,2));
    end
    centroid_list(i,1:2) = mean(r_list,1);
   
end

out = centroid_list;