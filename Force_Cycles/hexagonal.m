function [A,coords] = hexagonal(xlength,ylength)

[A,coords] = triangular(xlength,ylength);

for x = (-xlength/2):3:(xlength/2)
    parts = transpose([find(coords(:,1)==x);find(coords(:,1)==x+1.5)]);
    if ~isempty(parts)
        for p = parts
            A(p,:) = zeros(1,length(A));
            A(:,p) = zeros(length(A),1);
        end
    end
end