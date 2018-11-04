function out = loop_neighbors(c_list,A)

num_cycles = length(c_list(:,1));
cyc_length = length(c_list(1,:));
ext_links = zeros(num_cycles,100);

for i = 1:num_cycles
    list = zeros(1,100);
    count = 0;
    for j = 1:cyc_length
        part = c_list(i,j);
        links = find(A(part,:)~=0);
        list(count+1:count+length(links)) = links;
        count = count + length(links);
    end
    list = sort(unique(list(list>0)));
    num_links = length(list);
    for k = 1:num_links
        if ~isempty(intersect(list(k),c_list(i,:)))
            list(k) = 0;
        end
    end
    list = list(list > 0);
    num_links = length(list);
    ext_links(i,1:num_links) = list;
end

stop = 0;
n = 0;
while ~stop
    n = n + 1;
    if isempty(find(
end


out = ext_links;