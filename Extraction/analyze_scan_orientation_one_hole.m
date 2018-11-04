function out_as=analyze_scan_orientation_one_hole(imagefolder,imageprefix,start_image,end_image,x1,x2,y1,y2,radius,sigma0,AR_z,local_sph_IND,h)
tic
splits = 3;
%KERNEL Radius is larger than physical radius
kx=round(2.5*radius);%these are the matrix dimesions of the box that contains the kernel
ky=round(2.5*radius);
kz=round(2.5*radius/AR_z);%AR_z=1 currently
no_images=end_image-start_image+1;
if (mod(x2-x1,2) == 0 && mod(y2-y1,2) == 0) && (mod(kx,2)==0)
    kx = kx+1;
    ky = ky+1;
    kz = kz+1;
end

disp('***********************************************');
disp('This is the analyze_scan function');
%%
%First create the Gauss_sphere
disp('a_s: Creating (non z-deformed(!)) Gaussian sphere');
Gauss_sph=single(Gauss_sphere(radius,sigma0,kx,ky,kz,AR_z));
%%
%Load all images%creates 3d image array and makes a cropped version to
%line up with final moments
disp('a_s: Loading and pre-processing images');
[IMS,bit]=load_images_simple(start_image,end_image,x1,x2,y1,y2,imagefolder,imageprefix);
%%
%Threshold and invert
disp('a_s: Thresholding and inverting images');
% By threshold, suppress top 5% of bright pixels to the intensity value
% corresponding to the 95th percentile
IMS = thresh_invert(IMS,bit);

%%
% Band-filter all images
disp('a_s: Band-filtering all images');
%bpass_jhw cuts off radius*2 from EACH side in x
%so the x and y dimensions are reduced by !!4*radius!!
Cr=2*radius;%the amount that the original image is cropped by in x and y on either side by the bandpass
IMSbp=single(zeros(size(IMS,1)-2*Cr,size(IMS,2)-2*Cr,no_images));
for b=1:no_images
    IMSbp(:,:,b)=single(bpass_jhw(IMS(:,:,b),1,Cr));
end
IMSCr=max(max(max(IMSbp)))-IMSbp;
%% create IMSCr (thresholded bandpassed image)
%IMSCr is an inverted copy of the bandpassed img since IMSbp gets convolved
%IMSCr has a two-peaked disribution; one that corresponds to the bright
%background and one that corresponds to dark solid grains; the threshold is
%determined by locating the local minima between those two peaks
[hst,bins] = hist(IMSCr(:),100);
dffs = diff(hst);
dffs2 = diff(dffs);
bins = bins(2:end-1);
thres_val = round(bins(find(dffs2 < 0, 1,'last')));
%thres_val = 800; % !! data needs to be clearly distinct between liquid and solid; can be improved w/ index matching, exposure time
IMSCr = IMSCr > thres_val;

%%
%Convolve
disp('a_s: Convolving... This may take a while');
Convol=single(jcorr3d(IMSbp,Gauss_sph,splits));
sizekernel=size(Gauss_sph);
sC=size(Convol);%convolve does change size by adding a Kernel radius on either side of all dimmensions
Convol=Convol( round(sizekernel(1)/2) : round(sC(1) - sizekernel(1)/2) ...
            ,  round(sizekernel(2)/2) : round(sC(2) - sizekernel(2)/2) ...
            ,  round(sizekernel(3)/2) : round(sC(3) - sizekernel(3)/2)  );%crops (radius of the kernel*2) in order to bring Covol back to the same dimensions as the bandpassed image

sIMSCr=size(IMSCr); 
sC=size(Convol);
if sIMSCr~=sC
    disp('error: size of bandpassed image ~= to convolved+cropped image')
    return
end
%% create pkswb (thresholded convolution image)
disp('a_s: Thresholding');
pksbw = Convol > 0;
%%
disp('a_s: Tagging regions');
L=bwlabeln(pksbw);%output is an array with size of pksbw where all touching pixels in the 3d array have the same id number, an integer)
disp('a_s: Imposing Volume minimum of 4');
Resultunf=regionprops(L,'Area');%[NOTE L is array of TAGGED regions]; creates structure Resultunf with one 1x1 matricies(in a col) that are the areas of the tagged regions (sequentially by tag #) 
idx=find([Resultunf.Area]>=4);%index of all regions with nonzero area
L2=ismember(L,idx);%output is array with size L of 1's where elements of L are in the set idx~which is just 1:number of regions. Therefore it converts all tagged regions to all 1's
L3=bwlabeln(L2);% L3 now retaggs (L3=old L2)
%%
disp('a_s: Determining weighted centroid locations and orientations');
s=regionprops(L3,'PixelIdxList', 'PixelList');%s is a struct that holds structs for each tagged 
%region. the 2nd level struct holds two matricies: pxlidlist is the linear indicies
%of the nonzero pxls in that region. pxllist is the
%coordinates of each pxl in that region. NOTE:these indicies apply to the
%bandpassed image
%% 
Result=zeros(numel(s),11);
for k = 1:numel(s)%#elements in s (#regions or particles)
    idx = s(k).PixelIdxList;%lin index of all points in region k
    pixel_values = double(Convol(idx)+.0001);%list of values of the pixels in convol which has size of idx
    sum_pixel_values = sum(pixel_values);   
    x = s(k).PixelList(:, 1);%the list of x-coords of all points in the region k WITH RESPECT TO the bandpassed image
    y = s(k).PixelList(:, 2);
    z = s(k).PixelList(:, 3);
    xbar = sum(x .* pixel_values)/sum_pixel_values + Cr;%PLUS Cr BECAUSE
    ybar = sum(y .* pixel_values)/sum_pixel_values + Cr;%I CUT OFF Cr OF THE IMAGE DURING BANDPASS!(in x and y only) AND cropped kernelradius/2 off each side (in x,y,z)(but it was put back)                                                         %cropped radius/2 of each side
    zbar = sum(z .* pixel_values)/sum_pixel_values;
    x2moment = sum((x - xbar + Cr).^2 .* pixel_values) / sum_pixel_values;%+2*radius is added to translate the x coord(ie xbar has already been translated)
    y2moment = sum((y - ybar + Cr).^2 .* pixel_values) / sum_pixel_values;%these are with respto the translated image(ie the original IMS)
    z2moment = sum((z - zbar).^2 .* pixel_values) / sum_pixel_values;%the pixelvalues and sum of pixvalues are taken from the corresponding points in the bandpassed image. only the location has been translated
    x3moment = sum((x - xbar + Cr).^3 .* pixel_values) / sum_pixel_values;
    y3moment = sum((y - ybar + Cr).^3 .* pixel_values) / sum_pixel_values;
    z3moment = sum((z - zbar).^3 .* pixel_values) / sum_pixel_values;
    xskew = x3moment/(x2moment)^(1.5);
    yskew = y3moment/(y2moment)^(1.5);
    zskew = z3moment/(z2moment)^(1.5);
    Result(k,1:3) = [xbar+x1-1 ybar+y1-1 zbar+start_image-1]; % Set x,y,z coordinates relative to image coordinate system
    Result(k,4)   = max(pixel_values);
    Result(k,5)   = sum_pixel_values;
    Result(k,6:8) = [xskew yskew zskew];
    Result(k,9:11)=zeros(1,3);

    xbar = xbar - Cr; % Translate back to cropped image (now realigned with IMSCr, pksbw, convol,and IMSbp)                     
    ybar = ybar - Cr;
    if round(ybar)-radius<1||round(ybar)+radius>size(IMSCr,1)||round(xbar)-radius<1||round(xbar)+radius>size(IMSCr,2)||round(zbar)-radius<1||round(zbar)+radius>size(IMSCr,3)||(Result(k,5)<80)%if this true, orientation_finder wont be able to find a bead's pixel collection that fits within dimensions of IMSCr
        continue;
    end
    bead_box=IMSCr(round(ybar)-round(radius):round(ybar)+round(radius),round(xbar)-round(radius):round(xbar)+round(radius),round(zbar)-round(radius):round(zbar)+round(radius));%this will produce box (2*r+1)^3, with 1's in voxels that correspond to liquid
    bead_box(local_sph_IND)=0;
    
    
    % Add code so that the connected component through the center is
    % selected
    tagged_box=bwlabeln(bead_box);
    area_strc=regionprops(tagged_box,'Area');
    if isempty(area_strc)
        continue;
    end
    area_list_unsorted=[area_strc.Area]';
    area_list=sortrows([area_strc.Area]',-1);%list of areas of regions decending order
    axis_pts_idx=find(area_list_unsorted==area_list(1));
    bead_box_reg=ismember(tagged_box,axis_pts_idx);
    sor=regionprops(bead_box_reg,'PixelIdxList', 'PixelList');
    long_axis_IND = sor.PixelIdxList;
    
    %long_axis_IND = find(bead_box);
    [ylist,xlist,zlist] = ind2sub(size(bead_box),long_axis_IND);
    cent_guess = [mean(ylist) mean(xlist) mean(zlist)];
    N = length(xlist);
    covar = zeros(3);
    covar(1,1) = sum((xlist - cent_guess(2)).^2)/N; 
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
    Result(k,9:11) = Nlong';
    if Result(k,8) < 0
        Result(k,9:11) = -Result(k,9:11);
    end
    
end  
%%   
% disp('a_s: Removing invalid boundary points and missing orientations');
r= Result(:,1)<x1+radius   | Result(:,1) > x2-radius;
Result(r,:)=[];
r= Result(:,2)<y1+radius   | Result(:,2) > y2-radius;
Result(r,:)=[];
r= Result(:,3)<start_image+radius/AR_z | Result(:,3) > end_image-radius/AR_z;
Result(r,:)=[];
r= Result(:,9)==0 & Result(:,10)==0 & Result(:,11)==0;
Result(r,:)=[];
%%
disp('a_s: DONE, the analyze_scan function has ended.');
disp('***********************************************');
toc
out_as=Result;