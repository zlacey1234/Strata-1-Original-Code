function [A,coords] = cubic3d(xlength,ylength,zlength)

lboundx = -(xlength-1)/2;
uboundx = (xlength-1)/2;
lboundy = -(ylength-1)/2;
uboundy = (ylength-1)/2;
lboundz = -(zlength-1)/2;
uboundz = (zlength-1)/2;
for x = lboundx:1:uboundx
    for y = lboundy:1:uboundy
        for z = lboundz:1:uboundz
            disp(['x = ' num2str(x) ', y = ' num2str(y) ', z = ' num2str(z)]);
            coords = [coords; x y z];
        end
    end
end

% [coords1,coords2,coords3] = meshgrid(lboundx:uboundx,lboundy:uboundy,lboundz:uboundz);
% coords = [coords1;coords2;coords3];

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