function out = link_centroids(c_list)

num_contacts = length(c_list(:,1));
centroid_list = zeros(num_contacts, 3);

for i = 1:num_contacts
    centroid_list(i,1:3) = (c_list(i,3:5)+c_list(i,6:8))/2;
end

out = unique(centroid_list,'rows');