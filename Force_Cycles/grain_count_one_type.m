function out = grain_count_one_type(c_list)

min_id = min(min(min(c_list)));
max_id = max(max(max(c_list)));

count = zeros(max_id,1);

for i = min_id:max_id
    count(i) = count(i) + length(find(c_list == i));
end

out = count;