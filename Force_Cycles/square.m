function [A,coords] = square(xlength,ylength)

num_parts = (xlength+1)*(ylength+1);
coords = zeros(num_parts,2);
part = 0;

for x = (-xlength/2):(xlength/2)
    for y = (-ylength/2):(ylength/2)
        part = part + 1;
        coords(part,:) = [x y];
    end
end

A = zeros(length(coords(:,1)));

for i = 1:length(coords(:,1))-1
    for j = i+1:length(coords(:,1))
        dr = coords(i,:) - coords(j,:);
        if norm(dr) < 1.01 || (i <= (ylength+1) && j == (xlength)*(ylength+1)+i) || (rem(i-1,xlength+1)==0 && j == i+ylength)
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end
