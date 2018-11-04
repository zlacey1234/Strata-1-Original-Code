function nu = fund_cycles(A)

n = length(A(:,1));
T = []; % running list of all nodes considered
G = zeros(size(A)); % running list of connected subgraphs
num_subgraph = 0;

for i = 1:n
    if sum(T == i) > 0
        continue;
    end
    num_subgraph = num_subgraph + 1;
    C = i;
    T = sort(unique([T i]));
    ni = find(A(i,:)~=0);
    C = sort([C ni]);
    T = sort(unique([T ni]));
    while ~isempty(ni)
        C_old = C;
        for j = ni
            nj = find(A(j,:)~=0);
            C = sort(unique([C nj]));
            T = sort(unique([T nj]));
        end
        ni = setdiff(C,C_old);
    end
    G(num_subgraph,1:length(C)) = C;
end
G(num_subgraph+1:end,:) = [];
G = sparse(G);

nu = 0;
for c = 1:num_subgraph
    nc = find(G(c,:));
    A_part = A(nc,nc);
    E = 1/2*nnz(A_part);
    N = length(nc);
    nu = nu + E - N + 1;
end