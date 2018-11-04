function [sigma,R] = Rotation3D(sigma0,axis,ang_disp)

% sigma0 = column vector with original orientation
% axis = rotation axis (unit vector)
% ang_disp = angular displacement about rotation axis


a = cos(ang_disp/2);
b = axis(1)*sin(ang_disp/2);
c = axis(2)*sin(ang_disp/2);
d = axis(3)*sin(ang_disp/2);

R = [a^2+b^2-c^2-d^2 2*(b*c+a*d) 2*(b*d-a*c);...
    2*(b*c-a*d) a^2+c^2-b^2-d^2 2*(c*d+a*b);...
    2*(b*d+a*c) 2*(c*d-a*b) a^2+d^2-b^2-c^2];

sigma = R*sigma0;