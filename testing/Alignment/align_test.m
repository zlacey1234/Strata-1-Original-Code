function [s,dist_s,z,dist_z,th,dist_th] = align_test(varargin)
% INPUTS:
% part_list: particle list including list of unit norm orientations
% orient_col: column in which x-component of orientaion is listed;
% orientation information is thus stored in columns orient_col,
% orient_col+1, and orient_col+2
% t: min and max time over which to calculate alignment and determine
% distributions of mean z and mean theta
% time_col: column in which frame number is stored
% plot: optional; generate plots of distributions (0, no; 1, yes)

% OUTPUTS:
% 
% dist: distribution of S over time, as well as distributions of the mean
% of z and theta, as measured in each individual frame

% USAGE:
% [s,dist_s,z,dist_z,th,dist_th] = align_test(part_list,orient_col,t_min,t_max,time_col,p)

% A consistently low value of S~0 in the shear zone would indicate there is
% no tendency for the particles to align as a result of their drilled holes 
% or slight asphericity

part_list = varargin{1};
orient_col = varargin{2};
t_min = varargin{3};
t_max = varargin{4};
time_col = varargin{5};
if nargin == 6
    p = varargin{6};
else
    p = 0;
end

s = zeros(t_max-t_min+1,1);
z = zeros(t_max-t_min+1,1);
th = zeros(t_max-t_min+1,1);

for t = t_min:t_max
    disp(['time = ' num2str(t)]);
    [s(t-t_min+1),~,~] = align_param(part_list,t,orient_col,time_col);
    %x(t) = mean(orient_list(:,1));
    %y(t) = mean(orient_list(:,2));
    z(t-t_min+1) = mean(part_list(part_list(:,time_col)==t,orient_col+2));
    th(t-t_min+1) = mean(cart2pol(part_list(part_list(:,time_col)==t,orient_col),part_list(part_list(:,time_col)==t,orient_col+1),part_list(part_list(:,time_col)==t,orient_col+2)));
end

bins_s = 0:0.0006:0.03;
hist_s = hist(s,bins_s);
dist_s = hist_s/sum(hist_s)/0.0006;
mu_s = mean(s);
std_s = std(s);
bins_z = 0.45:0.001:0.55;
hist_z = hist(z,bins_z);
dist_z = hist_z/sum(hist_z)/0.001;
bins_th = -0.2:0.008:0.2;
hist_th = hist(th,bins_th);
dist_th = hist_th/sum(hist_th)/0.008;

if p
    figure(6)
    plot(t_min:t_max,s,'b-');
    hold on
    plot(t_min:t_max,mean(s)*ones(t_max-t_min+1,1),'k-','LineWidth',2);
    hold off
    xlabel('Trial');
    ylabel('S');
    axis square

    figure(7)
    plot(bins_s,dist_s,'r-');
    hold on
    plot(mu_s*ones(1,2),[0 max(dist_s)+20],'k-','LineWidth',2);
    plot((mu_s+std_s)*ones(size(0:10:max(dist_s)+20)),0:10:max(dist_s)+20,'k.','LineWidth',2);
    plot((mu_s-std_s)*ones(size(0:10:max(dist_s)+20)),0:10:max(dist_s)+20,'k.','LineWidth',2);
%     mode_s = bins_s(dist_s == max(dist_s));
%     ind_s = find(dist_s == max(dist_s));
%     std_s1 = std(s(s <= mode_s));
%     std_s2 = std(s(s >= mode_s));
%     plot(bins_s(1:ind_s),(sqrt(2/pi)/(std_s1+std_s2))*exp(-(bins_s(1:ind_s)-mode_s).^2./(2*std_s1.^2)),'k-');
%     plot(bins_s(ind_s:end),(sqrt(2/pi)/(std_s1+std_s2))*exp(-(bins_s(ind_s:end)-mode_s).^2./(2*std_s2.^2)),'k-');
    hold off
    xlabel('S');
    ylabel('p(S)');
    axis square

    figure(8)
    plot(bins_z,dist_z,'b-');
    hold on
    plot(mean(z)*ones(1,2),[0 max(dist_z)+20],'k-','LineWidth',2);
    plot((mean(z)+std(z))*ones(size(0:10:max(dist_z))),0:10:max(dist_z),'k.','LineWidth',2);
    plot((mean(z)-std(z))*ones(size(0:10:max(dist_z))),0:10:max(dist_z),'k.','LineWidth',2);
    plot(bins_z,(1/(sqrt(2*pi)*std(z)))*exp(-(bins_z-mean(z)).^2./(2*std(z).^2)),'k-');
    hold off
    xlabel('z_{bar}');
    ylabel('p(z_{bar})');
    axis([0.45 0.55 0 max(dist_z)+20]);
    axis square
    
    figure(9)
    plot(bins_th,dist_th,'r-');
    hold on
    plot(mean(th)*ones(1,2),[0 max(dist_th)+5],'k-','LineWidth',2);
    plot((mean(th)+std(th))*ones(size(0:2:max(dist_th))),0:2:max(dist_th),'k.','LineWidth',2);
    plot((mean(th)-std(th))*ones(size(0:2:max(dist_th))),0:2:max(dist_th),'k.','LineWidth',2);
    plot(bins_th,(1/(sqrt(2*pi)*std(th)))*exp(-(bins_th-mean(th)).^2./(2*std(th).^2)),'k-');
    hold off
    xlabel('\theta_{bar}');
    ylabel('p(\theta_{bar})');
    axis([-0.2 0.2 0 max(dist_th)+5]);
    axis square

% bins_xy = -1:0.001:1;
% hist_x = hist(x,bins_xy);
% dist_x = hist_x/sum(hist_x)/0.001;
% figure(f+2)
% plot(bins_xy,dist_x,'r-');
% hold on
% plot(mean(x)*ones(1,2),[0 max(dist_x)+20],'k-','LineWidth',2);
% plot((mean(x)+std(x))*ones(size(0:10:max(dist_x)+20)),0:10:max(dist_x)+20,'k.','LineWidth',2);
% plot((mean(x)-std(x))*ones(size(0:10:max(dist_x)+20)),0:10:max(dist_x)+20,'k.','LineWidth',2);
% hold off
% xlabel('x_{bar}');
% ylabel('p(x_{bar})');
% axis([-0.05 0.05 0 max(dist_x)+20]);
% 
% hist_y = hist(y,bins_xy);
% dist_y = hist_y/sum(hist_y)/0.001;
% figure(f+3)
% plot(bins_xy,dist_y,'g-');
% hold on
% plot(mean(y)*ones(1,2),[0 max(dist_y)+20],'k-','LineWidth',2);
% plot((mean(y)+std(y))*ones(size(0:10:max(dist_y)+20)),0:10:max(dist_y)+20,'k.','LineWidth',2);
% plot((mean(y)-std(y))*ones(size(0:10:max(dist_y)+20)),0:10:max(dist_y)+20,'k.','LineWidth',2);
% hold off
% xlabel('y_{bar}');
% ylabel('p(y_{bar})');
% axis([-0.05 0.05 0 max(dist_y)+20]);
end