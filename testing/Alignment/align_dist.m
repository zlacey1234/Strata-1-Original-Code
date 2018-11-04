function [dist_x,dist_y,dist_z,dist_th] = align_dist(varargin)
% INPUTS:
% part_list: particle list including list of unit norm orientations
% t: frame number(s) at which alignment order parameter is computed; if t
% is a vector, cumulative distributions are computed
% orient_col: column in which x-component of orientaion is listed;
% orientation information is thus stored in columns orient_col,
% orient_col+1, and orient_col+2
% time_col: column in which frame number is stored
% plot: optional; generate plots of distributions (0, no; 1, yes)

% OUTPUTS:
% dist: distribution of x, y, z, and theta components of orientation
% vectors; note: since |R^2 d\theta d(cos(\phi))| is the area element of a
% spherical surface, distribution of cos(phi) is proper to report, which is
% the same as the distribution of z

% USAGE:
% [dist_x,dist_y,dist_z,dist_th,dist_ph] = align_dist(part_list,t,orient_col,time_col,p)

% A consistently low value of s~0 in the shear zone would indicate there is
% no tendency for the particles to align as a result of their drilled holes 
% or slight asphericity

part_list = varargin{1};
t = varargin{2};
t = unique(t);
orient_col = varargin{3};
time_col = varargin{4};
if nargin < 5
    p = 0;
else
    p = varargin{5};
end

fstart = 6;

part_t = [];
for i = 1:length(t)
    part_t = [part_t; part_list(part_list(:,time_col)==t(i),:)];
end

bin_width = 0.01;
bin_w_th = pi/100;
bins_xy = -1:bin_width:1;
bins_z = 0:(bin_width/2):1;
bins_th = -pi:bin_w_th:pi;

hist_x = hist(part_t(:,orient_col),bins_xy);
hist_y = hist(part_t(:,orient_col+1),bins_xy);
hist_z = hist(part_t(:,orient_col+2),bins_z);
dist_x = hist_x/sum(hist_x)/bin_width;
dist_y = hist_y/sum(hist_y)/bin_width;
dist_z = hist_z/sum(hist_z)/(bin_width/2);


[th_list,~,~] = cart2sph(part_t(:,orient_col),part_t(:,orient_col+1),part_t(:,orient_col+2));
hist_th = hist(th_list,bins_th);
dist_th = hist_th/sum(hist_th)/bin_w_th;

if p
    figure(fstart);
    plot(bins_xy,dist_x,'r');
    hold on
    plot(bins_xy,mean(dist_x(2:end-1))*ones(size(bins_xy)),'k-','LineWidth',2);
    hold off
    axis square
    figure(fstart+1);
    plot(bins_xy,dist_y,'g');
    hold on
    plot(bins_xy,mean(dist_y(2:end-1))*ones(size(bins_xy)),'k-','LineWidth',2);
    hold off
    axis square
    figure(fstart+2);
    plot(bins_z,dist_z,'b');
    hold on
    plot(bins_z,mean(dist_z(2:end-1))*ones(size(bins_z)),'k-','LineWidth',2);
    hold off
    axis square
    figure(fstart+3);
    plot(bins_th,dist_th,'r');
    hold on
    plot(bins_th,mean(dist_th(2:end-1))*ones(size(bins_th)),'k-','LineWidth',2);
    hold off
    axis square
end