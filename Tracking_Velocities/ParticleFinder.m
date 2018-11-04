function [x,y,t,ang] = ParticleFinder(inputnames,threshold,framerange,outputname,bground_name,minarea,invert,noisy)
% Usage: [x,y,t,ang] = ParticleFinder(inputnames,threshold,[framerange],[outputname],[bground_name],[minarea],[invert],[noisy])
% Given a movie of particle motions, ParticleFinder identifies the 
% particles, returning their positions, times, and orientations in x, y, 
% and t, respectively. The movie must be saved as a series of image files, an 
% image stack in .tif or .gif format, or an uncompressed .avi file; specify 
% the movie in "inputnames" (e.g., '0*.png' or 'stack.tif', or 'movie.avi'). 
% To be identified as a particle, a part of the image must have brightness 
% that differs from the background by at least "threshold". If invert==0, 
% ParticleFinder seeks particles brighter than the background; if invert==1, 
% ParticleFinder seeks particles darker than the background; and if 
% invert==-1, ParticleFinder seeks any sort of contrast. The background is 
% read from the file "bground_name"; see BackgroundImage. Frames outside the 
% range specified by the two-element vector "framerange" are ignored. If 
% minarea==1, ParticleFinder seeks single-pixel particles by comparing 
% brightness to adjacent pixels (fast and good for small particles); otherwise ParticleFinder seeks particles 
% having areas larger than "minarea" (in square pixels; this method is 
% better for tracking large particles). If "outputname" is not empty, 
% particle positions are also saved as a binary file of that name. The 
% file begins with a long int giving the frame count, then each frame 
% begins with a long int giving its particle count, and continues with 
% double-precision floats giving the x and y coordinates of each particle. 
% If noisy~=0, the movie is repeated with particle locations overlaid. If 
% noisy>1, each movie frame is also saved to disk as an image. See also 
% BackgroundImage.m and PredictiveTracker.m. Requires 
% read_uncompressed_avi.m for use with .avi movies. 

% Written 20 October 2011 by Doug Kelley, largely based on PredictiveTracker.m.
% Renamed ParticleFinder and incorporated FindParticles function 27 October 2011. 
% Fixed plotting bug in bug plotting 15 November 2011. 
% Added "ang" output (particle orientation) 18 November 2011. 
% Added "framerange" input 28 November 2011. 
% Updated to use weighted centroid 2 December 2011. 
% Included invert==-1 option 22 February 2012.
% Made compatible with tiff & gif stacks and squelched regionprops 
% divide-by-zero warning 7 March 2012. 

% Next: write angle to output file!

% -=- Set defaults -=-----------------------------------------------------
framerange_default = [1 inf]; % by default, all frames
bground_name_default = 'background.tif';
noisy_default=0; % don't plot unless requested
minarea_default=1; 
invert_default=-1; % by default, use absolute contrast
pausetime=1/30; % seconds to pause between frames when plotting
savedirname='particlesmovie';

if nargin<2
    error(['Usage: [x,y,t,ang] = ' mfilename ...
        '(inputnames,threshold,[framerange],[outputname],[bground_name],' ...
        '[minarea],[invert],[noisy])'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerange_default;
end
if ~exist('bground_name','var') || isempty(bground_name)
    bground_name=bground_name_default;
end
if ~exist('minarea','var') || isempty(minarea)
    minarea=minarea_default;
end
if ~exist('invert','var') || isempty(invert)
    invert=invert_default;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisy_default;
end
if ~exist('outputname','var') || isempty(outputname)
    writefile=false;
else
    writefile=true;
end

% -=- Decide whether avi, stack, or images; set up -=---------------------
[filepath,junk,ext]=fileparts(inputnames);
names=dir(inputnames);
if strcmpi(ext,'.avi')
    movtype='avi';
    if isempty(which('read_uncompressed_avi.m')) % check for req'd helper function
        error(['Sorry, reading .avi files requires ' ...
            'read_uncompressed_avi.m.'])
    end
    movinfo=aviinfo(names.name);
    color_depth=movinfo.NumColormapEntries;
    ht=movinfo.Height;
    wd=movinfo.Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) movinfo.NumFrames]);
elseif numel(names)==1 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) % single file, looks like an image stack
    movtype='stack';
    movinfo=imfinfo(fullfile(filepath,names.name));
    color_depth=2^(movinfo(1).BitDepth);
    ht=movinfo(1).Height;
    wd=movinfo(1).Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) movinfo.NumFrames]);
else
    movtype='images';
    movinfo=imfinfo(fullfile(filepath,names(1).name));
    color_depth=2^(movinfo.BitDepth);
    ht=movinfo.Height;
    wd=movinfo.Width;
    for ii=1:numel(names);
        if strcmp(fullfile(filepath,names(ii).name),bground_name)
            names(ii)=[]; % don't try to track the background file
            break
        end
    end
    filenum=NaN(size(names));
    for ii=1:numel(names)
        [junk,myname]=fileparts(names(ii).name);
        filenum(ii)=str2double(myname);
    end
    if any(isnan(filenum)) % at least one filename was not a number, so just make ordinals
        filenum=1:numel(filenum);
    end
    tmin=max([framerange(1) min(filenum)]);
    tmax=min([framerange(2) max(filenum)]);
end
Nf=tmax-tmin+1;

% -=- Pre-compute logarithms for locating particle centers -=-------------
if minarea==1
    logs = 1:color_depth;
    logs = [log(0.0001) log(logs)];
end

% -=- Read background image -=--------------------------------------------
if exist(bground_name,'file')==2
    background = double(imread(bground_name));
else
    warning('MATLAB:ParticleFinder:noBackgroundFile', ...
        ['Cannot find background file ' bground_name ...
        '. Using blank image instead.'])
    if invert==1
        background = color_depth*ones(ht,wd);
    else
        background = zeros(ht,wd);
    end
end % if exist(bground_name,'file')==2
if ndims(background)==3 % check for RGB instead of grayscale
    background=round(mean(background,3));
end

% -=- Set up output file and variables -=---------------------------------
if writefile
    if exist(outputname,'file')
        yn=input(['File ' outputname ' exists. Overwrite (y/n)? '],'s');
        if ~strcmpi(yn(1),'y')
            disp('Not overwritten.')
            writefile=false;
        else
            disp(['Replacing file ' outputname '.'])
            fid=fopen(outputname,'w');
            fwrite(fid,Nf,'int32'); % header is integer frame count
        end
        pause(1)
    else
        fid=fopen(outputname,'w');
        fwrite(fid,Nf,'int32'); % header is integer frame count
    end
end
x=[];
y=[];
t=[];
ang=[];
begins=ones(Nf+1,1);
mark=[];
for tt=tmin:tmax-3
    switch movtype
        case('avi')
            if tt<tmin
                continue
            end
            [im,mark] = read_uncompressed_avi( ...
                fullfile(filepath,names.name),tt,mark);
            im = double(im.cdata);
        case('stack')
            if tt<tmin
                continue
            end
            im = double(imread(fullfile(filepath,names.name),1));
        case('images')
            if filenum(tt)<tmin
                continue
            end
            im = double(imread(fullfile(filepath,names(tt).name)));
    end % switch movtype
    if ndims(im)==3
        im=round(mean(im,3)); % convert to grayscale if necessary
    end
    if invert==1 % dark particles on light background
        im = background - im;
    elseif invert==0 % light particles on dark background
        im = im - background;
    else
        im = abs( im - background ); % seek any contrast, light or dark
    end
    im(im<0) = 0;
    if minarea==1
        pos = FindParticles(im,threshold,logs);
    else
        [pos,ang1] = FindRegions(im,threshold,minarea); % could also keep angle...
    end
    N=size(pos,1);
    x=[x ; pos(:,1)];
    y=[y ; pos(:,2)];
    t=[t ; tt*ones(N,1)];
    if minarea~=1
        ang=[ang ; ang1];
    end
    begins(tt-tmin+2)=begins(tt-tmin+1)+N;
    switch movtype
        case('avi')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' names.name ...
                ' frame ' num2str(tt) ' (' num2str(tt-tmin+1) ' of ' ...
                num2str(Nf) ').'])
        case('stack')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' names.name ...
                ' frame ' num2str(tt) ' (' num2str(tt-tmin+1) ' of ' ...
                num2str(Nf) ').'])
        case('images')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' names(tt).name ...
            ' (' num2str(tt-tmin+1) ' of ' num2str(Nf) ').'])
    end % switch movtype
    if writefile
        fwrite(fid,N,'int32'); % start the frame w/ particle count
        fwrite(fid,pos','float32'); % then write particle locations x1 y1 x2 y2  ...
        % Next: write angle as well!
    end
end % for tt=1:tmax
if writefile
    fclose(fid);
end

% -=- Plot if requested -=-
if noisy
    if isnumeric(noisy) && noisy>1
        disp(['Plotting and saving frames. ' ...
            'Please do not cover the figure window!'])
        if exist(savedirname,'file')~=7
            mkdir(savedirname)
        end
    else
        disp('Plotting...')
    end
    figure;
    axes('nextplot','add','dataaspectratio',[1 1 1],'ydir','reverse', ...
        'xlim',0.5+[0 size(im,2)],'ylim',0.5+[0 size(im,1)]);
    switch movtype
        case('avi')
            xlabel([names.name ', threshold = ' num2str(threshold) ...
                ', ' bground_name ', minarea = ' num2str(minarea), ...
                ', invert = ' num2str(invert)],'interpreter','none');
            [im,mark] = read_uncompressed_avi( ...
                fullfile(filepath,names.name),1);
            im = im.cdata;
        case('stack')
            xlabel([names.name ', threshold = ' num2str(threshold) ...
                ', ' bground_name ', minarea = ' num2str(minarea), ...
                ', invert = ' num2str(invert)],'interpreter','none');
            im = imread(fullfile(filepath,names.name),1);
        case('images')
            xlabel([inputnames ', threshold = ' num2str(threshold) ...
                ', ' bground_name ', minarea = ' num2str(minarea)], ...
                'interpreter','none');
            im = imread(fullfile(filepath,names(1).name));
    end % switch movtype
    hi=imagesc(im);
    colormap(gray); % has no effect if image is color b/c RGB values override
    hl=plot(NaN,NaN,'r.');
    mark=[];
    for ii=1:Nf
        ind=begins(ii):begins(ii+1)-1;
        set(hl,'xdata',x(ind),'ydata',y(ind));
        switch movtype
            case('avi')
                [im,mark] = read_uncompressed_avi( ...
                    fullfile(filepath,names.name),tmin+ii-1,mark);
                set(hi,'cdata',im.cdata);
                set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
                    names.name ' frame ' num2str(tmin+ii-1) ' (' ...
                    num2str(ii) ' of ' num2str(Nf) ')']); 
            case('stack')
                set(hi,'cdata',imread(fullfile(filepath,names.name), ...
                    tmin+ii-1)); 
                set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
                    names.name ' frame ' num2str(tmin+ii-1) ' (' ...
                    num2str(ii) ' of ' num2str(Nf) ')']); 
            case('images')
                set(hi,'cdata',imread(fullfile(filepath, ...
                    names(tmin+ii-1).name))); 
                set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
                    names(tmin+ii-1).name ' (' num2str(ii) ' of ' ...
                    num2str(Nf) ')']); 
        end % switch movtype
        drawnow
        pause(pausetime); 
        if isnumeric(noisy) && noisy>1
            snap=getframe(gca);
            imwrite(snap.cdata, ...
                fullfile(filepath,savedirname,names(ii).name)); 
        end 
    end
end
disp('Done.')

end % function ParticleFinder

% -=- function FindParticles -=-------------------------------------------
function pos = FindParticles(im, threshold, logs)
% Given an image "im", FindParticles finds small particles that are 
% brighter than their four nearest neighbors and also brighter than
% "threshold". Particles are located to sub-pixel accuracy by applying a 
% Gaussian fit in each spatial direction. The input "logs" depends on 
% the color depth and is re-used for speed. Particle locations are 
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second). 

    s = size(im);

    % identify the local maxima that are above threshold  
    maxes = find(im >= threshold & ...
        im > circshift(im,[0 1]) & ...
        im > circshift(im,[0 -1]) & ...
        im > circshift(im,[1 0]) & ...
        im > circshift(im,[-1 0]));

    % now turn these into subscripts
    [x,y] = ind2sub(s, maxes);

    % throw out unreliable maxes in the outer ring
    good = find(x~=1 & y~=1 & x~=s(1) & y~=s(2));
    x = x(good);
    y = y(good);

    % find the horizontal positions

    % look up the logarithms of the relevant image intensities
    z1 = logs(im(sub2ind(s,x-1,y)) + 1)';
    z2 = logs(im(sub2ind(s,x,y)) + 1)';
    z3 = logs(im(sub2ind(s,x+1,y)) + 1)';

    % compute the centers
    xcenters = -0.5 * (z1.*(-2*x-1) + z2.*(4*x) + z3.*(-2*x+1)) ./ ...
        (z1 + z3 - 2*z2);

    % do the same for the vertical position
    z1 = logs(im(sub2ind(s,x,y-1)) + 1)';
    z3 = logs(im(sub2ind(s,x,y+1)) + 1)';
    ycenters = -0.5 * (z1.*(-2*y-1) + z2.*(4*y) + z3.*(-2*y+1)) ./ ...
        (z1 + z3 - 2*z2);

    % make sure we have no bad points
    good = find(isfinite(xcenters) & isfinite(ycenters));

    % fix up the funny coordinate system used by matlab
    pos = [ycenters(good), xcenters(good)];

end % function FindParticles

% -=- function FindRegions -=---------------------------------------------
function [pos,ang] = FindRegions(im,threshold,minarea)
% Given an image "im", FindRegions finds regions that are brighter than
% "thresold" and have area larger than "minarea". Region centroids are 
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second). Region orientations are 
% returned in radians, in the vector "ang".

    s = size(im);
    warnstate=warning('off','MATLAB:divideByZero'); % squelch regionprops divide-by-zero warnings
    props=regionprops(im>threshold,im,'WeightedCentroid','Area', ...
        'Orientation','PixelList');
    pos=reshape([props.WeightedCentroid],2,numel(props))';
    ang=[props.Orientation]';
    if numel(pos)>0
        good = pos(:,1)~=1 & pos(:,2)~=1 & pos(:,1)~=s(2) & ...
            pos(:,2)~=s(1) & [props.Area]'>minarea; % remove regions on edge or too small
        pos=pos(good,:);
        ang=-ang(good)/180*pi; % convert to radians and right-handed coordinates
    end
    warning(warnstate) % put things back the way we found them
end % function FindRegions

