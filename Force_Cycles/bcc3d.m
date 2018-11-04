function [A,coords] = bcc3d(xlength,ylength,zlength)

coords = [];

lboundx = -(xlength-1)/2;
uboundx = (xlength-1)/2;
lboundy = -(ylength-1)/2;
uboundy = (ylength-1)/2;
lboundz = -(zlength-1)/2;
uboundz = (zlength-1)/2;

z = lboundz;
unit = 1;

while z <= uboundz
    if unit
        for x = lboundx:1:uboundx
            for y = lboundy:1:uboundy
                coords = [coords; x y z];
            end
        end
    else
        for x = (lboundx+0.5):1:(uboundx-0.5)
            for y = (lboundy+0.5):1:(uboundy-0.5)
                coords = [coords; x y z];
            end
        end
    end
    unit = ~unit;
    z = z + 1/sqrt(2);
end

A = zeros(length(coords(:,1)));

for i = 1:length(coords(:,1))-1
    for j = i+1:length(coords(:,1))
        disp(['i = ' num2str(i) ', j = ' num2str(j)]);
        dr = coords(i,:) - coords(j,:);
        if norm(dr) < 1.01 
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end