function out_as=analyze_scan_TESTorientation(IMS,x1,x2,y1,y2,radius,sigma0,AR_z,local_sph_IND,h)
%tic
splits = 3;
%KERNEL Radius is larger than physical radius
kx=round(2.5*radius);%these are the matrix dimesions of the box that contains the kernel
ky=round(2.5*radius);
kz=round(2.5*radius/AR_z);%AR_z=1 currently
% if (mod(x2-x1,2) == 0 && mod(y2-y1,2) == 0) && (mod(kx,2)==0)
%     kx = kx+1;
%     ky = ky+1;
%     kz = kz+1;
% end
disp('***********************************************');
disp('This is the analyze_scan function');
%%
%First create the Gauss_sphere
disp('a_s: Creating (non z-deformed(!)) Gaussian sphere');
Gauss_sph=single(Gauss_sphere(radius,sigma0,kx,ky,kz,AR_z));

%% 
% Invert the image stack
% Presume the data is 16 bit, as it always is
IMS = 2^16-1-IMS;
%%
% Band-filter all images
Cr=2*radius;%the amount that the original image is cropped by in x and y on either side by the bandpass
IMSbp=single(zeros(size(IMS,1)-2*Cr,size(IMS,2)-2*Cr,size(IMS,3)));
disp('a_s: Band-filtering all images');
for b=1:size(IMS,3)
    IMSbp(:,:,b)=single(bpass_jhw(IMS(:,:,b),0,Cr));
end
IMSCr=max(max(max(IMSbp)))-IMSbp;
%% create IMSCr (thresholded bandpassed image)
%IMSCr is an inverted copy of the bandpassed img since IMSbp gets convolved
%IMSCr has a two-peaked disribution; one that corresponds to the bright
%background and one that corresponds to dark solid grains; the threshold is
%determined by locating the local minima between those two peaks
[hst,bins] = hist(IMSCr(:),100);
dffs = diff(hst);
thres_val = round(bins(find(dffs < 0, 1,'last')+1));
IMSCr = IMSCr > thres_val; %thresholding value may change for different frames

%%
%Convolve
disp('a_s: Convolving... This may take a while');
Convol=single(jcorr3d(IMSbp,Gauss_sph,splits));

sizekernel=size(Gauss_sph);
sC=size(Convol);%convolve does change size by adding a Kernel radius on either side of all dimmensions
% if mod(sizekernel(1),2) == 0
%     endpt = sizekernel(1)/2;
% else
%     endpt = floor(sizekernel(1)/2);
% end
%Convol = Convol(endpt+1:end-endpt,endpt+1:end-endpt,endpt+1:end-endpt);
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
%TWEAK THRESHOLD
pksbw = Convol > 0;
%%
disp('a_s: Tagging regions');
L=bwlabeln(pksbw);%output is an array with size of pksbw where all touching pixels in the 3d array have the same id number, an integer)
disp('a_s: Imposing Volume minimum');
Resultunf=regionprops(L,'Area');%[NOTE L is array of TAGGED regions]; creates structure Resultunf with one 1x1 matricies(in a col) that are the areas of the tagged regions (sequentially by tag #) 
idx=find([Resultunf.Area]>1);%index of all regions with nonzero area
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
Result=zeros(numel(s),8);
for k = 1:numel(s);%#elements in s (#regions or particles)
%for k = 4236
    idx = s(k).PixelIdxList;%lin index of all points in region k
    pixel_values = double(Convol(idx)+.0001);%list of values of the pixels in convol which has size of idx
    sum_pixel_values = sum(pixel_values);   
    x = s(k).PixelList(:, 1);%the list of x-coords of all points in the region k WITH RESPECT TO the bandpassed image
    y = s(k).PixelList(:, 2);
    z = s(k).PixelList(:, 3);
    xbar = sum(x .* pixel_values)/sum_pixel_values + Cr;%PLUS Cr BECAUSE
    ybar = sum(y .* pixel_values)/sum_pixel_values + Cr;%I CUT OFF Cr OF THE IMAGE DURING BANDPASS!(in x and y only) AND cropped kernelradius/2 off each side (in x,y,z)(but it was put back)                                                         %cropped radius/2 of each side
    zbar = sum(z .* pixel_values)/sum_pixel_values     ;
%     x2moment = sum((x - xbar + Cr).^2 .* pixel_values) / sum_pixel_values;%+2*radius is added to translate the x coord(ie xbar has already been translated)
%     y2moment = sum((y - ybar + Cr).^2 .* pixel_values) / sum_pixel_values;%these are with respto the translated image(ie the original IMS)
%     z2moment = sum((z - zbar).^2 .* pixel_values) / sum_pixel_values;%the pixelvalues and sum of pixvalues are taken from the corresponding points in the bandpassed image. only the location has been translated
%     x3moment = sum((x - xbar + Cr).^3 .* pixel_values) / sum_pixel_values;
%     y3moment = sum((y - ybar + Cr).^3 .* pixel_values) / sum_pixel_values;
%     z3moment = sum((z - zbar).^3 .* pixel_values) / sum_pixel_values;
%     xskew = x3moment/(x2moment)^(1.5);
%     yskew = y3moment/(y2moment)^(1.5);
%     zskew = z3moment/(z2moment)^(1.5);
    Result(k,1:3) = [xbar+x1-1 ybar+y1-1 zbar];
    Result(k,4)   = max(pixel_values);
    Result(k,5)   = sum_pixel_values;
    %Result(k,6:8) = zeros(1,3);%[xskew yskew zskew];
    Result(k,6:8)=zeros(1,3);

    xbar = xbar - Cr; % Translate back to cropped image (now realigned with IMSCr, pksbw, convol,and IMSbp)                     
    ybar = ybar - Cr;
    if round(ybar)-radius<1||round(ybar)+radius>size(IMSCr,1)||round(xbar)-radius<1||round(xbar)+radius>size(IMSCr,2)||round(zbar)-radius<1||round(zbar)+radius>size(IMSCr,3)||(Result(k,5)<80)%if this true, orientation_finder wont be able to find a bead's pixel collection that fits within dimensions of IMSCr
        continue;
    end
    bead_box=IMSCr(round(ybar)-round(radius):round(ybar)+round(radius),round(xbar)-round(radius):round(xbar)+round(radius),round(zbar)-round(radius):round(zbar)+round(radius));%this will produce box (2*r+1)^3, with 1's in voxels that correspond to liquid
    bead_box(local_sph_IND)=0;%Both IMSCr and local_list have z indicies which start at 1.
    
    long_axis_IND = find(bead_box);
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
    Result(k,6:8) = Nlong';
    if Result(k,8) < 0
        Result(k,6:8) = -Result(k,6:8);
    end
end  
%%   
disp('a_s: Removing leftover invalid boundary points and missing orientations');%result x y z are IN IMS
% % Manually set these; at what point should there be no grains?
r= Result(:,1)<x1+radius | Result(:,1) > x2-radius;%r has size of Result and is logical
Result(r,:)=[];
r= Result(:,2)<x1+radius | Result(:,2) > y2-radius;
Result(r,:)=[];
r= Result(:,3)< radius | Result(:,3) > (size(IMS,3)-radius/AR_z );
Result(r,:)=[];
r= Result(:,6)==0 & Result(:,7)==0 & Result(:,8)==0;
Result(r,:)=[];

%%
disp('a_s: DONE, the analyze_scan function has ended.');
disp('***********************************************');
%toc
out_as=Result;