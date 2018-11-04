function [sigma,R] = Rotation3D(sigma0,axis,ang_disp)

% sigma0 = column vector with original orientation
% axis = rotation axis (unit vector)
% ang_disp = angular displacement about rotation axis


a = cos(ang_disp/2);
b = axis(1)*sin(ang_disp/2);
c = axis(2)*sin(ang_disp/2);
d = axis(3)*sin(ang_disp/2);

R = 2*[a^2+b^2-1/2 (b*c+a*d) (b*d-a*c);...
    (b*c-a*d) a^2+c^2-1/2 (c*d+a*b);...
    (b*d+a*c) (c*d-a*b) a^2+d^2-1/2];

sigma = R*sigma0;