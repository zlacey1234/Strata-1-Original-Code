function [vr,vs] = rollingSlidingDecomp(r1,r2,v1,v2,omega1,omega2,R1,R2)

% Takes in info of two grains in contact and computes their rolling and
% sliding components, based on Wang, et al Particuology 2015

% INPUTS:
% r1,r2: The positions of grains 1 and 2, respectively, in units of pixels
% v1,v2: The translational velocity of grains 1 and 2, respectively, in
% units of pixels / frame, and in the form of a 3D vector
% omega1,omega2: The angular velocity of grains 1 and 2, respectively, in
% units of radians / frame, and in the form of a 3D vector that points
% along rotation axis, with norm equal to the angular velocity
% R1,R2: The radii of grains 1 and 2, respectively, in units of pixels; for
% monodisperse systems, these should be set to be the same

% OUTPUTS:
% vr: the rolling velocity
% vs: the sliding velocity
% both in pixels / frame

% Generate unit vector from grain 1 to 2
dr = r2 - r1;
n = dr/norm(dr);

% Calculate components of v1 and v2 that are tangential to point of contact
v1t = v1 - dot(v1,n)*n;
v2t = v2 - dot(v2,n)*n;

% Objective velocities of the grains
s1 = R1*cross(omega1,n) - R1*(v2t - v1t)/(R1 + R2);
s2 = -R2*cross(omega2,n) + R2*(v2t - v1t)/(R1 + R2);

% Compute the sliding and rolling components of these
s1r = (R2*s1 + R1*s2) / (R1 + R2);
% s2r = s1r;
s1s = -R1/(R1 + R2)*(s2 - s1);
s2s =  R2/(R1 + R2)*(s2 - s1);

vr = s1r;
vs = s2s - s1s;