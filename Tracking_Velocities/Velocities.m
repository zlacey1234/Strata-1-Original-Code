function [u,v,x,y,t]=Velocities(vtracks,framerange,noisy)
% Usage: [u,v,x,y,t]=Velocities(vtracks,[framerange],[noisy])
% Working from the velocity tracks in "vtracks", Velocities plots the
% velocity of each particle in frames "framerange". Velocity values are 
% returned in vectors "u" and "v"; corresponding positions and times are 
% returned in "x", "y", and "t". If noisy==0, no plot is produced. The input 
% must be a structure of the form produced by PredictiveTracker.m. Specify 
% "framerange" as a two-element vector: [starttime endtime]. This file can be 
% downloaded from http://leviathan.eng.yale.edu/software.

% Written 13 April 2011 by Douglas H. Kelley.
% Fixed frame range in title when using default 1 August 2011. 

framerangedefault = [-inf inf]; % all frames
noisydefault=1;

if nargin<1
    error(['Usage: [u,v,x,y,t] = ' mfilename '(vtracks,[framerange],[noisy]'])
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
[t,ind]=sort(t(:)); % sort by time
u=u(ind);
v=v(ind);
x=x(ind);
y=y(ind);
u=u(:); % output column vectors
v=v(:);
x=x(:);
y=y(:);

% If requested, plot the particles
if noisy
    figure;
    quiver(x,y,u,v);
    set(gca,'dataaspectratio',[1 1 1]);
    axis tight
    if framerange(1)==framerange(2)
        title([num2str(numel(x)) ' particles in frame ' ...
            num2str(framerange(1))]);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title([num2str(numel(x)) ' particles in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange)) ]);
    end
end
