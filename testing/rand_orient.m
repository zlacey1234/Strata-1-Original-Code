function part_list = rand_orient(n,t_min,t_max)
% INPUTS:
% n: number of particles
% t: min and max frame number

% OUTPUTS:
% part_list: list of randomly oriented points on the unit hemisphere 

part_list = zeros(n*(t_max-t_min+1),4);
for t = t_min:t_max
    i = t-t_min+1;
    theta = 2*pi*rand(n,1);
    phi = asin(rand(n,1));
    [part_list(n*(i-1)+1:n*i,1),part_list(n*(i-1)+1:n*i,2),part_list(n*(i-1)+1:n*i,3)]=sph2cart(theta,phi,ones(n,1));
    part_list(n*(i-1)+1:n*i,4)=t*ones(n,1);
end