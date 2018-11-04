clc
close all

cyclestart = 0;
cycleend = 0;
framestart=56;
frameend=56;
z_list=190;
f = 0;
radius = 15;

for z = z_list
    f = f + 1;
    for cyc = cyclestart:cycleend
        for fr=framestart:1:frameend
            im=im2double(imread(['H:\RotationsExperiments\steady2deg_rotations09JAN13\frame' num2str(fr,'%04.0f') '\shear' num2str(fr,'%04.0f') num2str(z,'%05.0f') '.tif']));
            figure(f);
            imagesc(im);
            axis equal;
            colormap('gray');
            hold on
            set(f,'Position',[80,80,1200,900]);
            idx=find((large_t(:,3)+11 >= z-11) & (large_t(:,3)+11 <= z+11) & (large_t(:,14) == fr));
            % Plot circles showing particle centers
            for i=1:length(idx)
                diffz=2*abs(z-large_t(idx(i),3));
                angle=sin(diffz/11);
                [cx,cy]=jcircle(11*cos(angle));
                plot(large_t(idx(i),1)+cx-1,large_t(idx(i),2)+cy-1,'r');
            end
%           Plot arrows showing orientation         
            dir = radius*large_t(:,9:11);
            for j = 1:length(idx)
                endpts = [large_t(idx(j),1)-dir(idx(j),1) large_t(idx(j),2)-dir(idx(j),2) large_t(idx(j),3)-dir(idx(j),3); large_t(idx(j),1)+dir(idx(j),1) large_t(idx(j),2)+dir(idx(j),2) large_t(idx(j),3)+dir(idx(j),3)];
                xpts = linspace(endpts(1,1),endpts(2,1),50);
                ypts = linspace(endpts(1,2),endpts(2,2),50);
                zpts = linspace(endpts(1,3),endpts(2,3),50);
                index = find((zpts < z+5) & (zpts > z-5));
                plot(xpts(index),ypts(index),'g-');
            end
        end
    end
    hold off
end