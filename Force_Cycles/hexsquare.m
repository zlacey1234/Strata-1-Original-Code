function [A,coords] = hexsquare(xlength,ylength)

[A,coords] = triangular(xlength,ylength);

for x = (-xlength/2):2:(xlength/2)
    parts = transpose(find(coords(:,1)==x));
    if ~isempty(parts)
        for p = parts
            A(p,:) = zeros(1,length(A));
            A(:,p) = zeros(length(A),1);
        end
    end
    parts2 = transpose(find(coords(:,1)==x+1/2));
    if ~isempty(parts)
        for q = parts2
            links = find(A(q,:)~=0);
            for l = links
                if coords(l,2) == coords(q,2) && coords(l,1) > coords(q,1)
                    A(l,q) = 0;
                    A(q,l) = 0;
                end
            end
        end
    end
end