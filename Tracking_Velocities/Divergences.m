function [div,x,y,t]=Divergences(vtracks,framerange,noisy)
% Usage: [div,x,y,t]=Divergences(vtracks,[framerange],[noisy])
% Working from the velocity tracks in "vtracks", Divergences plots the
% divergence of each particle in frames "framerange". Divergence values are 
% returned in the vector "div"; corresponding positions and times are 
% returned in "x", "y", and "t". If noisy==0, no plot is produced. The input 
% must be a struct array of the form produced by PredictiveTracker.m. 
% Specify "framerange" as a two-element vector: [starttime endtime]. See
% also Velocities.m, Vorticities.m. This file can be downloaded from 
% http://leviathan.eng.yale.edu/software.

% Written 17 May 2011 by Douglas H. Kelley.
% Fixed frame range in title when using default 5 August 2011. 

framerangedefault = [-inf inf]; % all frames
noisydefault=1;
satgoal=1; % colormap limits chosen for 1% saturation
minparticles = 10; % do not attempt calcs w/o at least this many particles

if nargin<1
    error(['Usage: [div,x,y,t] = ' mfilename ...
        '(vtracks,[framerange],[noisy]'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerangedefault;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end
fn=fieldnames(vtracks);
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'U')) || ~any(strcmp(fn,'V')) || isempty(vtracks)
    error('Sorry, the input does not appear to contain tracks.')
end
u=[]; % initialize arrays
v=[];
x=[];
y=[];
t=[];
for ii=1:numel(vtracks) % assemble arrays from vtracks structure
    ind = (vtracks(ii).T>=min(framerange)) & (vtracks(ii).T<=max(framerange));
    u=[u vtracks(ii).U(ind)];
    v=[v vtracks(ii).V(ind)];
    x=[x vtracks(ii).X(ind)];
    y=[y vtracks(ii).Y(ind)];
    t=[t vtracks(ii).T(ind)];
end
[t,ind]=sort(t); % sort by time
u=u(ind);
v=v(ind);
x=x(ind);
y=y(ind);

[tlist,ends,junk]=unique(t);
begins=circshift(ends,1)+1;
begins(1)=1;
div=NaN(size(x));
for ii=1:numel(tlist)
    ind=begins(ii):ends(ii);
    if numel(ind)<minparticles
        t(ind)=NaN;
        continue
    end
    tri=delaunayn([x(ind)' y(ind)']);
    [ux,junk]=pdegrad([x(ind) ; y(ind)],tri',u(ind)');
    [junk,vy]=pdegrad([x(ind) ; y(ind)],tri',v(ind)');
    div(ind)=pdeprtni([x(ind) ; y(ind)],tri',ux+vy);
end
ind=~isnan(t); % remove frames with too few particles
if ~sum(ind)
    error(['Sorry, no frames between ' num2str(framerange(1)) ' and ' ...
        num2str(framerange(2)) ' contain enough particles.'])
end
div=div(ind);
x=x(ind);
y=y(ind);
t=t(ind);
div=div(:); % output column vectors
x=x(:);
y=y(:);
t=t(:);

% If requested, plot the particles
if noisy
    figure;
    scatter(x,y,[],div,'o','filled','displayname','divergence');
    if which('redblue')
        set(gcf,'colormap',redblue);
    end
    colorbar
    set(gca,'dataaspectratio',[1 1 1], ...
        'clim',prctile(abs(div),100-satgoal)*[-1 1]);
    axis tight
    if framerange(1)==framerange(2)
        title(['divergence at ' num2str(numel(x)) ' particles in frame ' ...
            num2str(framerange(1))]);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title(['divergence at ' num2str(numel(x)) ' particles in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange)) ]);
    end
end

