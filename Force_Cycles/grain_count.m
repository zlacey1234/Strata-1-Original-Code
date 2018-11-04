function out = grain_count(c3_list,c4_list,c5_list,c6_list)

min_id = min([min(min(c3_list)) min(min(c4_list)) min(min(c5_list)) min(min(c6_list))]);
max_id = max([max(max(c3_list)) max(max(c4_list)) max(max(c5_list)) max(max(c6_list))]);

count = zeros(max_id,1);

for i = min_id:max_id
    if i < max(max(c3_list))
        count(i) = count(i) + length(find(c3_list == i));
    end
    if i < max(max(c4_list))
        count(i) = count(i) + length(find(c4_list == i));
    end
    if i < max(max(c5_list))
        count(i) = count(i) + length(find(c5_list == i));
    end
    if i < max(max(c6_list))
        count(i) = count(i) + length(find(c6_list == i));
    end
end

out = count;