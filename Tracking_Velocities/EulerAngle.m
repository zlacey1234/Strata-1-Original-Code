% Given two frames of orientation vectors sigma1 and sigma2, determines the
% complete Euler rotation of the grain from the first frame to the second.
% If only information about sigma1 is given, 2 of the 3 angles are
% computed.

function [phi,theta,psi] = EulerAngle(varargin)

% Inputs: sigma: orientation row vector (x y z) of grain; sigma1 is primary
% hole, sigma2 is secondary hole; initial and final states of sigma are
% passed, respectively (see usage for precise definition); over-rotation
% should be already taken care of in tracking algorithm
%
% Outputs: 
%
% Usage: EulerAngle(sigma1_init,sigma1_final,sigma2_init,sigma2_final)
%        EulerAngle(sigma1_init,sigma1_final)

% !! This needs to be edited to match notation for Euler angles on
% MathWorld !!

if nargin < 2
    error('Initial and final sigma1 orientations required');
end
sigma1_init = varargin{1};
sigma1_final = varargin{2};
% Rotate reference frame so that sigma1_init is in the +x-direction; then,
% theta and phi and simply given by sigma1_final coords in the new frame
% In this function, theta is azimuthal rotation, and phi is polar rotation
[theta1_init,phi1_init,~] = cart2sph(sigma1_init(1),sigma1_init(2),sigma1_init(3)); %[azimuth, elevation, ~]
% Rotate sigma1_init and sigma1_final by -theta1_init and -phi1_init
R_theta = [cos(theta1_init) sin(theta1_init) 0; -sin(theta1_init) cos(theta1_init) 0; 0 0 1];
R_phi = [cos(phi1_init) 0 sin(phi1_init); 0 1 0; -sin(phi1_init) 0 cos(phi1_init)];
sigma1_i_rotate = R_phi*R_theta*transpose(sigma1_init); % This should be [1;0;0], up to machine precision; Euler angles are given in x-notation
sigma1_f_rotate = R_phi*R_theta*transpose(sigma1_final);
[phi,psi,~] = cart2sph(sigma1_f_rotate(1),sigma1_f_rotate(2),sigma1_f_rotate(3));
% theta and phi range: [-pi, pi]

if nargin > 2
    if nargin == 3
        error('Initial and final sigma2 orientations required');
    elseif nargin > 4
        error('Too many inputs');
    end
    % If sigma2 is passed, put it into the same reference frame
    sigma2_init = varargin{3};
    sigma2_final = varargin{4};
    sigma2_i_rotate = R_phi*R_theta*sigma2_init';
    sigma2_f_rotate = R_phi*R_theta*sigma2_final';
    % Project sigma2_rotate by rotating it by -phi and -psi (angular
    % displacements of sigma1)
    R_azi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
    R_pol = [cos(psi) 0 sin(psi); 0 1 0; -sin(psi) 0 cos(psi)];
    sigma2_i_project = R_pol*R_azi*sigma2_i_rotate; % This should be ~[0;0;1] or ~[0;1;0]
    sigma2_f_project = R_pol*R_azi*sigma2_f_rotate;
    % Find azimuthal angle between projection and rotated initial sigma2
    % The sign is relative to direction of sigma1_i_rotate
    % Azimuthal angle is found to account for the fact that the notch
    % might not be perfectly perpendicular to long axis in general
    [theta,~,~] = EulerAngle(sigma2_i_project',sigma2_f_project');
else
    theta = 0; % rotational information about main hole is not accessible
end