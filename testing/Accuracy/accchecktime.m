%accchecktime
clc
close all

largetracks = result;
%frames=[4 max(largetracks(:,end-1))]; % time(s) you want to check 
z_list=80; % height(s) you want to check
%time_col = 7; % column that lists the time
%tol_s = 6.5; % tolerance in height--how far particle centers can be from the given height
tol_l = 9;
scale = 0.203;
r = 2.5/scale;
%amp = 40;
testscan = 1;

close all
for z = z_list
    %for time=frames
        %cycle = floor(time/amp)+1;
        %frame = time - (cycle-1)*amp;
        %if frame == 0
%             cycle = cycle - 1;
%             frame = amp;
%         end
        % load image
        if testscan
            im=im2double(imread(['L:\testscan_28SEP2015\test' num2str(z,'%05.0f') '.tif']));
        else
            %im=im2double(imread(['N:\segregation_steadyshear_2deg_09JAN12\frame' num2str(time,'%04.0f') '\shear' num2str(time,'%04.0f') num2str(z,'%05.0f') '.tif']));
        end
        figure;
        imagesc(im); % show original grayscale image
        axis equal;
        axis off;
        colormap gray;
        hold on
        set(gca,'position',[0 0 1 1],'units','normalized')
        %set(1,'Position',[80,80,900,1200]);
        % Small particles
        % find particles at given time within a certain height range
%         idx=find( (smalltracks(:,3) >= z-tol_s) & (smalltracks(:,3) <= z+tol_s) & (smalltracks(:,time_col)==time) );
%    
%         for i=1:length(idx)
%                 % draw a circle around the particle
%                 diffz=2*abs(z-smalltracks(idx(i),3));
%                 angle=sin(diffz/7);
%                 [cx,cy]=jcircle(7*cos(angle));
%                 plot(smalltracks(idx(i),1)+cx-1,smalltracks(idx(i),2)+cy-1,'b');
%         end    

        %Large particles
        % repeat for a second species, if necessary
        %idx=find( (largetracks(:,3) >= z-tol_l) & (largetracks(:,3) <= z+tol_l) & (largetracks(:,time_col)==time) );
        idx=find( (largetracks(:,3) >= z-tol_l) & (largetracks(:,3) <= z+tol_l));
        for i=1:length(idx)
            diffz=2*abs(1*(z-largetracks(idx(i),3)));
            angle=sin(diffz/11);
            [cx,cy]=jcircle(r*cos(angle));
            plot(largetracks(idx(i),1)+cx-1,largetracks(idx(i),2)+cy-1,'r');
        end
        
        imw=getframe;
        % write to file or make a movie
        %filename = ['C:\students\matt\WolfgangsMovie_4040\check' num2str(time,'%04.0f') 'z' num2str(z,'%03.0f') '.tif'];
        %print ('-dtiff', filename);
        %hold off
        %close all
   % end
end