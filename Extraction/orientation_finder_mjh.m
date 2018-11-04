function out_as=orientation_finder_mjh(IMSCr,local_ind,local_sph_IND,xbar,ybar,zbar,radius,hole_r)
% MJH's edits to ML's method for finding two hole orientations

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
bead_box=IMSCr(round(ybar)-radius:round(ybar)+radius,round(xbar)-radius:round(xbar)+radius,round(zbar)-radius:round(zbar)+radius);%this will produce box 31x31x31
bead_box_reg = bead_box;
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
if numel(area_list)>=3
    axis_pts_idx=[find(area_list_unsorted==area_list(1)); find(area_list_unsorted==area_list(2)); find(area_list_unsorted==area_list(3))];%index of all regions with highest 3 area
    axis_pts_idx=axis_pts_idx(1:3);
    bead_box_reg=ismember(tagged_box,axis_pts_idx);% it converts the three largest regions to all 1's and other regions to zeros
    tagged_box=bwlabeln(bead_box_reg);%now retaggs regions with region IDs ie 1, 2, 3
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
if sleng~=3
    %disp('program terminated: found not three light regions (post removal of small regions)');
    out_as=[0 0 0 0 0 0];
    return
end
%% find region centroids
axis_pts=zeros(3);
for n = 1:3;%#elements in s (#regions or particles) SHOULD BE 3
    idxor = sor(n).PixelIdxList;%lin index of all points in region n
    pixel_values = double(bead_box(idxor));%list of values of the pixels in convol which has size of idx
    sum_pixel_values = sum(pixel_values);   
    x = sor(n).PixelList(:, 1);%the list of x-coords of all points in the region n WITH RESPECT TO the beadbox
    y = sor(n).PixelList(:, 2);
    z = sor(n).PixelList(:, 3);
    Axbar = sum(x .* pixel_values)/sum_pixel_values;%
    Aybar = sum(y .* pixel_values)/sum_pixel_values;%
    Azbar = sum(z .* pixel_values)/sum_pixel_values;%
    axis_pts(n,:)=[Axbar-radius-1, Aybar-radius-1, Azbar-radius-1];%coordinates relative to bead box center. beadbox has center at 16,16,16 = (radius+1)x3
end
%plot3([ones(3,1)*0 axis_pts(:,1)],[ones(3,1)*0 axis_pts(:,2)],[ones(3,1)*0 axis_pts(:,3)],'ro')
%%

%% find axes
%axis_pts=[axis_pts(1,:)/norm(axis_pts(1,:)); axis_pts(2,:)/norm(axis_pts(2,:)); axis_pts(3,:)/norm(axis_pts(3,:));zeros(1,3)];%normalizes vectors
key=find([norm(cross(axis_pts(2,:),axis_pts(3,:))); norm(cross(axis_pts(1,:),axis_pts(3,:))); norm(cross(axis_pts(1,:),axis_pts(2,:)))]... 
    ==min([norm(cross(axis_pts(2,:),axis_pts(3,:))); norm(cross(axis_pts(1,:),axis_pts(3,:))); norm(cross(axis_pts(1,:),axis_pts(2,:)))]));%id for the combination of two vects that form the long (bead) axis
switch key
    case 3
        %[~, rod,~]=fit_3D_data(axis_pts([1,2,4],1),axis_pts([1,2,4],2),axis_pts([1,2,4],3),'line','off','off');
        rod_guess = axis_pts(1,:)-axis_pts(2,:);
        Nperp_cand = axis_pts(3,:);
    case 2
        %[~, rod,~]=fit_3D_data(axis_pts([1,3,4],1),axis_pts([1,3,4],2),axis_pts([1,3,4],3),'line','off','off');
        rod_guess = axis_pts(3,:)-axis_pts(1,:);
        Nperp_cand = axis_pts(2,:);
    otherwise
        %[~, rod,~]=fit_3D_data(axis_pts([2,3,4],1),axis_pts([2,3,4],2),axis_pts([2,3,4],3),'line','off','off');
        rod_guess = axis_pts(2,:)-axis_pts(3,:);
        Nperp_cand = axis_pts(1,:);
end

[theta,phi,~] = cart2sph(rod_guess(2),rod_guess(1),rod_guess(3));
phi = pi/2 - phi;
rod_test = newCylinder(2*hole_r,8/5*radius,theta,phi,1,0,size(bead_box,1));
correlate = bead_box.*rod_test;
long_axis_IND = find(correlate);
if length(long_axis_IND) < 900 % condition for truly finding cylinder could be adjusted
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

% cent_guess = round(cent_guess);
% shift = (radius+1)*ones(1,3)-cent_guess;
% bead_box_new = ones(size(bead_box));
% yrange = 1:2*radius+1;
% yrange2 = yrange;
% xrange = 1:2*radius+1;
% xrange2 = xrange;
% zrange = 1:2*radius+1;
% zrange2 = zrange;
% if shift(1) > 0
%     yrange(1:shift(1)) = [];
%     yrange2(2*radius+2-shift(1):end) = [];
% else
%     yrange(2*radius+2+shift:end) = [];
%     yrange2(1:-shift(1)) = [];
% end
% if shift(2) > 0
%     xrange(1:shift(2)) = [];
%     xrange2(2*radius+2-shift(2):end) = [];
% else
%     xrange(2*radius+2+shift(2):end) = [];
%     xrange2(1:-shift(2)) = [];
% end
% if shift(3) > 0
%     zrange(1:shift(3)) = [];
%     zrange2(2*radius+2-shift(3):end) = [];
% else
%     zrange(2*radius+2+shift(3):end) = [];
%     zrange2(1:-shift(3)) = [];
% end
% bead_box_new(yrange,xrange,zrange) = bead_box(yrange2,xrange2,zrange2);
% Erase the long axis, then just do pca on the leftover notch!!
bead_box_new = bead_box;
[theta_f,phi_f,~] = cart2sph(Nlong(2),Nlong(1),Nlong(3));
phi_f = pi/2-phi_f;
rod_f = newCylinder(3*hole_r,2*radius,theta_f,phi_f,1,0,size(bead_box,1));
bead_box_new(rod_f==1) = 0;
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
bead_box_new(local_sph_IND) = 0;
tagged_box=bwlabeln(bead_box_new);
% 7/29/14: could add minimum threshold for area of notch
area_strc=regionprops(tagged_box,'Area');%[NOTE L is array of TAGGED regions]; creates structure Resultunf with one 1x1 matricies(in a col) that are the areas of the tagged regions (sequentially by tag #) 
if isempty(area_strc)
    out_as = [Nlong' 0 0 0];
    return;
end
area_list_unsorted=[area_strc.Area]';
area_list=sortrows([area_strc.Area]',-1); 
axis_pts_idx=find(area_list_unsorted==area_list(1));
tagged_box=ismember(tagged_box,axis_pts_idx);
sor=regionprops(tagged_box,'PixelIdxList', 'PixelList');
idxor = sor.PixelIdxList;
pixel_values = double(bead_box_new(idxor));
sum_pixel_values = sum(pixel_values);
new_axis_pts = zeros(numel(sor),3);
cand_best = zeros(numel(sor),1);
for m = 1:numel(sor)
    x = sor(m).PixelList(:, 1);%the list of x-coords of all points in the region n WITH RESPECT TO the beadbox
    y = sor(m).PixelList(:, 2);
    z = sor(m).PixelList(:, 3);
    Axbar = sum(x .* pixel_values)/sum_pixel_values;%
    Aybar = sum(y .* pixel_values)/sum_pixel_values;%
    Azbar = sum(z .* pixel_values)/sum_pixel_values;%
    %new_axis_pts(m,:) = [Axbar-radius-1 Aybar-radius-1 Azbar-radius-1]-bead_box_shift;
    new_axis_pts(m,:) = [Axbar-radius-1 Aybar-radius-1 Azbar-radius-1];
    cand_best(m) = abs(dot(new_axis_pts(m,:),Nperp_cand))/norm(new_axis_pts(m,:))/norm(Nperp_cand);
end
ind_best = find(cand_best == max(cand_best));
if length(ind_best) > 1
    out_as = [Nlong' 0 0 0];
    return;
end
%     covar = zeros(3);
%     cent_notch = new_axis_pts(ind_best,:) + (radius+1)*ones(1,3);
%     xlist = sor(ind_best).PixelList(:,1);
%     ylist = sor(ind_best).PixelList(:,2);
%     zlist = sor(ind_best).PixelList(:,3);
%     N = length(sor(ind_best).PixelList(:,1));
%     covar(1,1) = sum((xlist-cent_notch(1)).^2)/N;
%     covar(1,2) = sum((xlist-cent_guess(1)).*(ylist-cent_guess(2)))/N;
%     covar(1,3) = sum((xlist-cent_guess(1)).*(zlist-cent_guess(3)))/N;
%     covar(2,1) = covar(1,2);
%     covar(2,2) = sum((ylist-cent_guess(2)).^2)/N;
%     covar(2,3) = sum((ylist-cent_guess(2)).*(zlist-cent_guess(3)))/N;
%     covar(3,1) = covar(1,3);
%     covar(3,2) = covar(2,3);
%     covar(3,3) = sum((zlist-cent_guess(3)).^2)/N;
%     COEFF = pcacov(covar);
%     Nperp_cand = (Nperp_cand / norm(Nperp_cand))';
%     Nperp_cand_mat = repmat(Nperp_cand,1,3);
%     Nperp_test = abs(dot(COEFF,Nperp_cand_mat,1));
%     Nperp = COEFF(:,Nperp_test==max(Nperp_test));

Nperp = new_axis_pts(ind_best,:) / norm(new_axis_pts(ind_best,:)); % We want rotations about the sphere CENTER, so use this definition for notch direction
temp = Nperp(2);
Nperp(2) = Nperp(1);
Nperp(1) = temp;
% if dot(Nperp,Nperp_cand) < 0
%     Nperp = -Nperp;
% end

out_as=[Nlong' Nperp];