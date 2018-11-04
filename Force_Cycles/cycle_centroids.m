function out = cycle_centroids(c_list,shearzone,t)

num_cycle = length(c_list(:,1));
c_length = length(c_list(1,:));
r_t = shearzone(shearzone(:,9) == t,:);

centroid_list = zeros(num_cycle, 3);

for i = 1:num_cycle
    r_list = [];
    for j = 1:c_length
        r_list = [r_list; r_t(r_t(:,10) == c_list(i,j),1:3)];
    end
    centroid_list(i,1:3) = mean(r_list,1);
end

out = centroid_list;