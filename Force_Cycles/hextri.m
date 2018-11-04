function [A,coords] = hextri(xlength,ylength)

[A,coords] = triangular(xlength,ylength);

for x = (-xlength/2):2:(xlength/2)
    parts = transpose(find(coords(:,1)==x));
    if ~isempty(parts)
        for p = parts
            A(p,:) = zeros(1,length(A));
            A(:,p) = zeros(length(A),1);
        end
    end
end