function out_as=orientation_finder_mjh(IMSCr,local_ind,local_sph_IND,xbar,ybar,zbar,radius,hole)
% MJH's method for locating orientations using the two-axis beads, building
% off of what ML started

% My idea is to replace the method of using three small regions, which
% don't put any constrains on the shape of the cavity

% After ML's method finds the two points that should correspond to the
% cylinder, I add in a test that is performed on the rods data, where a
% candidate cylinder is drawn with prescribed length (distance between the 
% two points) and hole radius (input: hole).  Then, once this refined
% measurement for the long axis orientation is found, these pixels are
% blacked out so that another round of pca can be performed on the most
% prominent bright region left (the notch).

%INPUT:

%OUTPUT: principle components organized into a matrix
    %% find cube around particle in Imscr

%bead_box_shift = [ybar-round(ybar) xbar-round(xbar) zbar-round(zbar)];
bead_box=IMSCr(round(ybar)-round(radius):round(ybar)+round(radius),round(xbar)-round(radius):round(xbar)+round(radius),round(zbar)-round(radius):round(zbar)+round(radius));%this will produce box (2*r+1)^3
bead_box_new = bead_box;
bead_box_new(local_sph_IND) = 0;
bead_box_reg = bead_box_new;
bead_box_reg(local_ind) = 0;
    %now bead box is zeros outside the shell
tagged_box=bwlabeln(bead_box_reg);%output is an array with size of pksbw where all touching pixels in the 3d array have the same id number, an integer)
%% erase small regions from tagged box
area_strc=regionprops(tagged_box,'Area');%[NOTE L is array of TAGGED regions]; creates structure Resultunf with one 1x1 matricies(in a col) that are the areas of the tagged regions (sequentially by tag #) 
if isempty(area_strc)
    out_as = [0 0 0 0 0 0];
    return;
end
area_list_unsorted=[area_strc.Area]';
area_list=sortrows([area_strc.Area]',-1);%list of areas of regions decending order
if numel(area_list)>=4 
    axis_pts_idx=[find(area_list_unsorted==area_list(1)); find(area_list_unsorted==area_list(2)); find(area_list_unsorted==area_list(3)); find(area_list_unsorted==area_list(4))];%index of all regions with highest 3 area
    axis_pts_idx = unique(axis_pts_idx,'stable');
    axis_pts_idx=axis_pts_idx(1:4);
    bead_box_reg=ismember(tagged_box,axis_pts_idx);% it converts the three largest regions to all 1's and other regions to zeros
    tagged_box=bwlabeln(bead_box_reg);%now retags regions with region IDs ie 1, 2, 3, ...
else
    out_as=[0 0 0 0 0 0];%disp('program terminated: found less than three light regions initially');
    return
end
%% get region properties
%disp('a_s: Determining weighted centroid locations and orientations');
sor=regionprops(tagged_box,'PixelIdxList', 'PixelList');%sor is a struct that holds structs for each tagged 
%region. the 2nd level struct holds two matricies: pxlidlist is the linear indicies
%of the nonzero pxls in that region. pxllist is the
%coordinates of each pxl in that region. NOTE:these indicies apply to the
%bandpassed image
sleng = numel(sor);
if sleng~=4
    %disp('program terminated: found not three light regions (post removal of small regions)');
    out_as=[0 0 0 0 0 0];
    return
end
%% find region centroids
axis_pts=zeros(4,3);
for n = 1:4;%#elements in s (#regions or particles) SHOULD BE 3
    idxor = sor(n).PixelIdxList;%lin index of all points in region n
    pixel_values = double(bead_box(idxor));%list of values of the pixels in convol which has size of idx
    sum_pixel_values = sum(pixel_values);   
    x = sor(n).PixelList(:, 1);%the list of x-coords of all points in the region n WITH RESPECT TO the beadbox
    y = sor(n).PixelList(:, 2);
    z = sor(n).PixelList(:, 3);
    Axbar = sum(x .* pixel_values)/sum_pixel_values;%
    Aybar = sum(y .* pixel_values)/sum_pixel_values;%
    Azbar = sum(z .* pixel_values)/sum_pixel_values;%
    axis_pts(n,:)=[Axbar-radius-1, Aybar-radius-1, Azbar-radius-1];%coordinates relative to bead box center. beadbox has center at (radius+1)x3
end
%plot3([ones(3,1)*0 axis_pts(:,1)],[ones(3,1)*0 axis_pts(:,2)],[ones(3,1)*0 axis_pts(:,3)],'ro')
%%

%% find axes
%axis_pts=[axis_pts(1,:)/norm(axis_pts(1,:)); axis_pts(2,:)/norm(axis_pts(2,:)); axis_pts(3,:)/norm(axis_pts(3,:));zeros(1,3)];%normalizes vectors
% Take the cross product of all pairs of patch centroids, pick the two that
% are most co-linear with the center, but also check that they don't share
% any common points; this provides a first guess for the cavity
% orientations (to be refined)
xes = [norm(cross(axis_pts(1,:),axis_pts(2,:))); norm(cross(axis_pts(1,:),axis_pts(3,:))); norm(cross(axis_pts(1,:),axis_pts(4,:))); norm(cross(axis_pts(2,:),axis_pts(3,:))); norm(cross(axis_pts(2,:),axis_pts(4,:))); norm(cross(axis_pts(3,:),axis_pts(4,:)));];
key=find(xes == min(xes),1,'first');%id for the combination of two vects that form the long (bead) axis
switch key
    case 1
        %[~, rod,~]=fit_3D_data(axis_pts([1,2,4],1),axis_pts([1,2,4],2),axis_pts([1,2,4],3),'line','off','off');
        rod_guess_1 = axis_pts(1,:)-axis_pts(2,:);
        if xes(6) ~= min(xes(2:6))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(3,:)-axis_pts(4,:);
    case 2
        %[~, rod,~]=fit_3D_data(axis_pts([1,3,4],1),axis_pts([1,3,4],2),axis_pts([1,3,4],3),'line','off','off');
        rod_guess_1 = axis_pts(1,:)-axis_pts(3,:);
        if xes(5) ~= min(xes([1 3:6]))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(2,:)-axis_pts(4,:);
    case 3
        %[~, rod,~]=fit_3D_data(axis_pts([2,3,4],1),axis_pts([2,3,4],2),axis_pts([2,3,4],3),'line','off','off');
        rod_guess_1 = axis_pts(1,:)-axis_pts(4,:);
        if xes(4) ~= min(xes([1:2 4:6]))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(2,:)-axis_pts(3,:);
    case 4
        rod_guess_1 = axis_pts(2,:)-axis_pts(3,:);
        if xes(3) ~= min(xes([1:3 5:6]))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(1,:)-axis_pts(4,:);
    case 5
        rod_guess_1 = axis_pts(2,:)-axis_pts(4,:);
        if xes(2) ~= min(xes([1:4 6]))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(1,:)-axis_pts(3,:);
    otherwise
        rod_guess_1 = axis_pts(3,:)-axis_pts(4,:);
        if xes(1) ~= min(xes(1:5))
            out_as = [0 0 0 0 0 0];
            return;
        end
        Nperp_cand = axis_pts(1,:)-axis_pts(4,:);
end

[theta,phi,~] = cart2sph(rod_guess_1(2),rod_guess_1(1),rod_guess_1(3));
phi = pi/2 - phi;
rod_test = newCylinder(3/2*hole,2*radius,theta,phi,1,0,size(bead_box,1));
correlate = bead_box.*rod_test;
long_axis_IND = find(correlate);
if length(long_axis_IND) < 200 % condition for truly finding cylinder could be adjusted
    out_as = [0 0 0 0 0 0];
    return;
end
[ylist,xlist,zlist] = ind2sub(size(correlate),long_axis_IND);
%[~,Nlong,~] = fit_3D_data(xlist,ylist,zlist,'line','off','off');
cent_guess = [mean(ylist) mean(xlist) mean(zlist)];
N = length(xlist);
covar = zeros(3);
covar(1,1) = sum((xlist-cent_guess(2)).^2)/N;
covar(1,2) = sum((xlist-cent_guess(2)).*(ylist-cent_guess(1)))/N;
covar(1,3) = sum((xlist-cent_guess(2)).*(zlist-cent_guess(3)))/N;
covar(2,1) = covar(1,2);
covar(2,2) = sum((ylist-cent_guess(1)).^2)/N;
covar(2,3) = sum((ylist-cent_guess(1)).*(zlist-cent_guess(3)))/N;
covar(3,1) = covar(1,3);
covar(3,2) = covar(2,3);
covar(3,3) = sum((zlist-cent_guess(3)).^2)/N;
COEFF = pcacov(covar);
Nlong = COEFF(:,1);
Nlong = (Nlong / norm(Nlong))';

[theta_perp,phi_perp,~] = cart2sph(Nperp_cand(2),Nperp_cand(1),Nperp_cand(3));
phi_perp = pi/2 - phi_perp;
rod_test_2 = newCylinder(7/4*hole,2*radius,theta_perp,phi_perp,1,0,size(bead_box,1));
correlate_2 = bead_box.*rod_test_2;
long_axis2_IND = find(correlate_2);
if length(long_axis2_IND) < 200 % condition for truly finding cylinder could be adjusted
    out_as = [0 0 0 0 0 0];
    return;
end
[ylist,xlist,zlist] = ind2sub(size(correlate_2),long_axis2_IND);
%[~,Nlong,~] = fit_3D_data(xlist,ylist,zlist,'line','off','off');
cent_guess = [mean(ylist) mean(xlist) mean(zlist)];
N = length(xlist);
covar = zeros(3);
covar(1,1) = sum((xlist-cent_guess(2)).^2)/N;
covar(1,2) = sum((xlist-cent_guess(2)).*(ylist-cent_guess(1)))/N;
covar(1,3) = sum((xlist-cent_guess(2)).*(zlist-cent_guess(3)))/N;
covar(2,1) = covar(1,2);
covar(2,2) = sum((ylist-cent_guess(1)).^2)/N;
covar(2,3) = sum((ylist-cent_guess(1)).*(zlist-cent_guess(3)))/N;
covar(3,1) = covar(1,3);
covar(3,2) = covar(2,3);
covar(3,3) = sum((zlist-cent_guess(3)).^2)/N;
COEFF = pcacov(covar);
Nperp = COEFF(:,1);
Nperp = (Nperp / norm(Nperp))';


%figure;
%IND = find(bead_box ==1);
%[X,Y,Z] = ind2sub(size(bead_box),IND);
%plot3(X,Y,Z,'k.')
%axis equal
%hold on
%IND = find(rod_test ==1);
%[X,Y,Z] = ind2sub(size(bead_box),IND);
%plot3(X,Y,Z,'r.')
%IND = find(rod_test_2 ==1);
%[X,Y,Z] = ind2sub(size(bead_box),IND);
%plot3(X,Y,Z,'g.')
%plot3([radius-radius*Nlong(2) radius+radius*Nlong(2)],[radius-radius*Nlong(1) radius+radius*Nlong(1)],[radius-radius*Nlong(3) radius+radius*Nlong(3)],'r-','LineWidth',8)
%plot3([radius-radius*Nperp(2) radius+radius*Nperp(2)],[radius-radius*Nperp(1) radius+radius*Nperp(1)],[radius-radius*Nperp(3) radius+radius*Nperp(3)],'g-','LineWidth',8)


% Erase the long axis, then just do pca on the leftover notch!!
%bead_box_new = bead_box;
%[theta_f,phi_f,~] = cart2sph(Nlong(2),Nlong(1),Nlong(3));
%phi_f = pi/2-phi_f;
%rod_f = newCylinder(3*hole,2*radius,theta_f,phi_f,1,0,size(bead_box,1));
%bead_box_new(rod_f==1) = 0;
% [X,Y,Z] = meshgrid(1:size(bead_box_new,1),1:size(bead_box_new,2),1:size(bead_box_new,3));
% X = X(1,:,1);
% Y = transpose(Y(:,1,1));
% Z = Z(1,1,:);
% Z = transpose(Z(:));
% local_sph = transpose(combvec(Y,X,Z));
% [~,~,R]=cart2sph((local_sph(:,2)-radius-1),(local_sph(:,1)-radius-1),(local_sph(:,3)-radius-1));
% local_sph=local_sph((R(:,1)>4/5*radius),:);
% clear R;
% local_sph_IND=sub2ind(size(bead_box_new),local_sph(:,1),local_sph(:,2),local_sph(:,3));
% bead_box_new(local_sph_IND) = 0;
% tagged_box=bwlabeln(bead_box_new);
% % 7/29/14: could add minimum threshold for area of notch
% area_strc=regionprops(tagged_box,'Area');%[NOTE L is array of TAGGED regions]; creates structure Resultunf with one 1x1 matricies(in a col) that are the areas of the tagged regions (sequentially by tag #) 
% if isempty(area_strc)
%     out_as = [Nlong' 0 0 0];
%     return;
% end
% area_list_unsorted=[area_strc.Area]';
% area_list=sortrows([area_strc.Area]',-1); 
% axis_pts_idx=find(area_list_unsorted==area_list(1));
% tagged_box=ismember(tagged_box,axis_pts_idx);
% sor=regionprops(tagged_box,'PixelIdxList', 'PixelList');
% idxor = sor.PixelIdxList;
% pixel_values = double(bead_box_new(idxor));
% sum_pixel_values = sum(pixel_values);
% new_axis_pts = zeros(numel(sor),3);
% cand_best = zeros(numel(sor),1);
% for m = 1:numel(sor)
%     x = sor(m).PixelList(:, 1);%the list of x-coords of all points in the region n WITH RESPECT TO the beadbox
%     y = sor(m).PixelList(:, 2);
%     z = sor(m).PixelList(:, 3);
%     Axbar = sum(x .* pixel_values)/sum_pixel_values;%
%     Aybar = sum(y .* pixel_values)/sum_pixel_values;%
%     Azbar = sum(z .* pixel_values)/sum_pixel_values;%
%     %new_axis_pts(m,:) = [Axbar-radius-1 Aybar-radius-1 Azbar-radius-1]-bead_box_shift;
%     new_axis_pts(m,:) = [Axbar-radius-1 Aybar-radius-1 Azbar-radius-1];
%     cand_best(m) = abs(dot(new_axis_pts(m,:),Nperp_cand))/norm(new_axis_pts(m,:))/norm(Nperp_cand);
% end
% ind_best = find(cand_best == max(cand_best));
% if length(ind_best) > 1
%     out_as = [Nlong' 0 0 0];
%     return;
% end
% %     covar = zeros(3);
% %     cent_notch = new_axis_pts(ind_best,:) + (radius+1)*ones(1,3);
% %     xlist = sor(ind_best).PixelList(:,1);
% %     ylist = sor(ind_best).PixelList(:,2);
% %     zlist = sor(ind_best).PixelList(:,3);
% %     N = length(sor(ind_best).PixelList(:,1));
% %     covar(1,1) = sum((xlist-cent_notch(1)).^2)/N;
% %     covar(1,2) = sum((xlist-cent_guess(1)).*(ylist-cent_guess(2)))/N;
% %     covar(1,3) = sum((xlist-cent_guess(1)).*(zlist-cent_guess(3)))/N;
% %     covar(2,1) = covar(1,2);
% %     covar(2,2) = sum((ylist-cent_guess(2)).^2)/N;
% %     covar(2,3) = sum((ylist-cent_guess(2)).*(zlist-cent_guess(3)))/N;
% %     covar(3,1) = covar(1,3);
% %     covar(3,2) = covar(2,3);
% %     covar(3,3) = sum((zlist-cent_guess(3)).^2)/N;
% %     COEFF = pcacov(covar);
% %     Nperp_cand = (Nperp_cand / norm(Nperp_cand))';
% %     Nperp_cand_mat = repmat(Nperp_cand,1,3);
% %     Nperp_test = abs(dot(COEFF,Nperp_cand_mat,1));
% %     Nperp = COEFF(:,Nperp_test==max(Nperp_test));

% Nperp = new_axis_pts(ind_best,:) / norm(new_axis_pts(ind_best,:)); % We want rotations about the sphere CENTER, so use this definition for notch direction
% temp = Nperp(2);
% Nperp(2) = Nperp(1);
% Nperp(1) = temp;
% if dot(Nperp,Nperp_cand) < 0
%     Nperp = -Nperp;
% end

%out_as=[Nlong' Nperp];
out_as = [Nlong Nperp];