function [A,coords] = triangular(xlength,ylength)

coords = [];

for x = (-xlength/2):(xlength/2)
    coords = [coords; x 0];
    y = 0;
    while y < ylength/2
        y = y + sqrt(3)/2;
        coords = [coords; x+1/2 -y; x+1/2 y];
        y = y + sqrt(3)/2;
        coords = [coords; x -y; x y];
    end
end

A = zeros(length(coords(:,1)));

for i = 1:length(coords(:,1))-1
    for j = i+1:length(coords(:,1))
        dr = coords(i,:) - coords(j,:);
        if norm(dr) < 1.01
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end