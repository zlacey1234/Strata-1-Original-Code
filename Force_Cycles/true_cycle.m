function out = true_cycle(loop,A,coords,xlength,ylength)

% Input:
% loop = candidate loop, in order that loop is traversed
% A = full adjacency matrix
% coords = coordinates of all particles centered at (0,0) (only for 2D)
% xlength = horizontal length of containing box
% ylength = vertical length of containing box

out = 1; % Assume the candidate loop is a fundamental loop, otherwise proven otherwise

cyc_length = length(loop);

% Basic checks previously done in cycle_finder.m (backtracking and diagonal
% contacts)

if length(unique(loop)) ~= cyc_length % Backtracking
    out = 0; % If 'out' is ever assigned a zero, stop the function immediately
    return;
else
    for k = 1:cyc_length % Diagonal contacts
        part = loop(k);
        contacts = find(A(part,:)~=0);
        if length(intersect(contacts,loop)) ~= 2
            out = 0;
            return;
        end
    end
end

% Final check for interior contacts is only for planar cycles
if ~isempty(coords)
    coords_loop = zeros(cyc_length,length(coords(1,:)));
    for v = 1:cyc_length
        part_v = loop(v);
        coords_loop(v,:) = coords(part_v,:);
    end

    if (length(unique(coords_loop(:,1))) == 1 || length(unique(coords_loop(:,2))) == 1 || sum(size(coords_loop(1,:))==[1 2])==2 || length(unique(coords_loop(:,3))) == 1)
        % Specify vertices of the loop
        vx = zeros(cyc_length,1);
        vy = zeros(cyc_length,1);
        for v = 1:cyc_length
            part_v = loop(v);
            if sum(size(coords_loop(1,:))==[1 2])==2 || length(unique(coords_loop(:,3))) == 1
                vx(v) = coords(part_v,1);
                vy(v) = coords(part_v,2);
            elseif length(unique(coords_loop(:,1))) == 1
                vx(v) = coords(part_v,2);
                vy(v) = coords(part_v,3);
            else
                vx(v) = coords(part_v,1);
                vy(v) = coords(part_v,3);
            end
        end

        % Adjust for periodic boundary conditions
        if max(vx) - min(vx) > 0.75*xlength
            vx(vx < 0) = vx(vx < 0) + xlength;
            coords(coords(:,1)<0,1) = coords(coords(:,1)<0,1) + (xlength+.01);
        end
        if max(vy) - min(vy) > 0.75*ylength
            vy(vy < 0) = vy(vy < 0) + ylength;
            coords(coords(:,2)<0,2) = coords(coords(:,2)<0,2) + (ylength+.01);
        end

        % Loop over all unique pairs of particles and find out if there are any
        % common links inside the loop
        for i = 1:(cyc_length-1)
            part1 = loop(i);
            n1 = find(A(part1,:)~=0);
            for j = (i+1):cyc_length
                part2 = loop(j);
                n2 = find(A(part2,:)~=0);
                commons = intersect(n1,n2); % List of common links
                if ~isempty(commons)
                    for part3 = commons
                        if sum(part3 == loop) > 0
                            continue; % Skip common links that are part of the loop
                        else
                            coords3 = coords(part3,:);
                            out = ~inpolygon(coords3(1),coords3(2),vx,vy);
                        end
                    end
                    if out == 0 % If 'out' is ever assigned a zero, stop the function immediately
                        return;
                    end
                end
            end
        end 
    end
end