function [cycles_new,cycles_dot] = planar_cyc(cycles,persist)

% cycles = list of 4- or 6-cycles
% persist = persistant network with coordinates

cyc_length = length(cycles(1,:));
if cyc_length ~= 4 && cyc_length ~= 6
    error('Improper cycle length')
end
count = 0;
cycles_new = [];
cycles_dot = zeros(length(cycles(:,1)),cyc_length+1);
for i = 1:length(cycles(:,1))
    parts = cycles(i,:);
    r1 = unique(persist(persist(:,1)==parts(1),3:5),'rows');
    r2 = unique(persist(persist(:,1)==parts(2),3:5),'rows');
    r3 = unique(persist(persist(:,1)==parts(3),3:5),'rows');
    plane_vec1 = r2 - r1;
    plane_vec2 = r3 - r1;
    normal = cross(plane_vec1,plane_vec2);
    unit_normal = normal / norm(normal);
    if cyc_length == 4
        r4 = unique(persist(persist(:,1)==parts(4),3:5),'rows');
        test_vec1 = (r4 - r1) / norm(r4 - r1);
        test_vec2 = (r4 - r2) / norm(r4 - r2);
        test_vec3 = (r4 - r3) / norm(r4 - r3);
        if dot(unit_normal,test_vec1) == 0 && dot(unit_normal,test_vec2) == 0 && dot(unit_normal,test_vec3) == 0
            count = count + 1;
            cycles_new(count,:) = cycles(i,:);
        end
        cycles_dot(i,:) = [cycles(i,:) sqrt(1-min([dot(unit_normal,test_vec1) dot(unit_normal,test_vec2) dot(unit_normal,test_vec3)])^2)];
    else
        error('Still working on 6-cycle code')
    end
end