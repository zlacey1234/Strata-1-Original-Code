function [A,coords] = hcp3d(xlength,ylength,zlength)

coords = [];
lboundx = -(xlength-1)/2;
uboundx = (xlength-1)/2;
lboundy = -(ylength-1)/2;
uboundy = (ylength-1)/2;
lboundz = -(zlength-1)/2;
uboundz = (zlength-1)/2;
z = lboundz;
unit = 0;

while z <= uboundz;
    if unit
        for x = lboundx:1:uboundx
            y = 0;
            coords = [coords; x y z];
            y = sqrt(3)/2;
            while y < uboundy
                coords = [coords; x+1/2 -y z; x+1/2 y z];
                y = y + sqrt(3)/2;
                if y < uboundy
                    coords = [coords; x -y z; x y z];
                    y = y + sqrt(3)/2;
                end
            end
        end
    else
        for x = lboundx:1:uboundx
            y = sqrt(3)/4;
            coords = [coords; x y z];
            y1 = y - sqrt(3)/2;
            y2 = y + sqrt(3)/2;
            while y1 > lboundy
                coords = [coords; x+1/2 y1 z];
                y1 = y1 - sqrt(3)/2;
                if y1 > lboundy
                    coords = [coords; x y1 z];
                    y1 = y1 - sqrt(3)/2;
                end
            end
            while y2 < uboundy
                coords = [coords; x+1/2 y2 z];
                y2 = y2 + sqrt(3)/2;
                if y2 < uboundy
                    coords = [coords; x y2 z];
                    y2 = y2 + sqrt(3)/2;
                end
            end
        end
    end
    z = z + 3/4;
    unit = ~unit;
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