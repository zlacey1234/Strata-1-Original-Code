function [A,coords] = octagon(xlength,ylength)

coords = [];
x = -xlength/2;
y = -ylength/2;

while y < ylength/2
    while x < xlength/2
        coords = [coords; x y; x y+1; x+1/sqrt(2) y+1+1/sqrt(2); ...
            x+1+1/sqrt(2) y+1+1/sqrt(2); x+1+sqrt(2) y+1; x+1+sqrt(2) y; ...
            x+1+1/sqrt(2) y-1/sqrt(2); x+1/sqrt(2) y-1/sqrt(2)];
        x = x + 2 + sqrt(2);
    end
    x = -xlength/2;
    y = y + 2 + sqrt(2);
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