function [s,T,v] = align_param(varargin)
% INPUTS:
% part_list: particle list including list of unit norm orientations
% t: frame number at which alignment order parameter is computed
% orient_col: column in which x-component of orientaion is listed;
% orientation information is thus stored in columns orient_col,
% orient_col+1, and orient_col+2
% time_col: column in which frame number is stored
% plot: optional; shows visual representation of orientations on unit 
% sphere (0, no; 1, yes)

% OUTPUTS:
% s: alignment order parameter at time t
% T: the associated traceless order tensor
% v: associated eigenvector

% USAGE:
% [s,T,v] = align_param(part_list,t,orient_col,time_col,plot)

% A consistently low value of s~0 in the shear zone would indicate there is
% no tendency for the particles to align as a result of their drilled holes 
% or slight asphericity

part_list = varargin{1};
t = varargin{2};
orient_col = varargin{3};
time_col = varargin{4};
if nargin < 5
    p = 0;
else
    p = varargin{5};
end

part_t = part_list(part_list(:,time_col)==t,:);
part_t = part_t(part_t(:,orient_col)~=0 & part_t(:,orient_col+1)~=0 & part_t(:,orient_col+2)~=0,:);
N = length(part_t(:,1)); % number of particles
T = -1/2*eye(3); % initialize the T tensor, including the Kronecker delta term
for i = 1:3
    for j = i:3
        T(i,j) = T(i,j) + 3/(2*N)*sum(part_t(:,orient_col+i-1).*part_t(:,orient_col+j-1));
        T(j,i) = T(i,j);
    end
end

[VT,ST] = eig(T);
ST = diag(ST);
s = max(ST);
ind = (ST == s);
v = VT(:,ind);

% visual representation of particle orientations on a unit sphere

if p
    figure(5)
    plot3(part_t(:,orient_col),part_t(:,orient_col+1),part_t(:,orient_col+2),'.');
    axis equal
end