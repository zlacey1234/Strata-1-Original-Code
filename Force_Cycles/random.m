function [A,coords] = random(xlength,ylength)

[A,coords] = triangular(xlength,ylength);

for i = 1:length(coords(:,1))
    links = find(A(i,:));
    for l = links
        if rand > 0.6
            A(i,l) = 0;
            A(l,i) = 0;
        end
    end
end