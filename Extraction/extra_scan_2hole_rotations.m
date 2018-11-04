clc
clear
%Do a run
disp('**************************************************');
disp('  Welcome to the particle & orientation extracter');
disp('**************************************************');
disp(' ');
testscan = 1; % set to 1 if you are just doing an initial analysis on a single volume
if ~testscan
    basefolder='Q:\reverse0505_2deg_rotations_15JAN13\'; % folder within which images are saved
    savefolder = 'C:\students\matt\analyses_2holerotation_steady_16JUN14\'; % folder within which temporary results are saved, in addition to full position list
    cycles = 20; % Set to -1 for steady shear
    frames = 40; %
else
    basefolder='O:\testscan_14SEP2015\'; % folder in which testscan is saved
    cycles = 0;
    frames = 1;
end
disp('Gathering information...');
AR_z=1; % aspect ratio in z-direction; experiments should be calibrated so that this is always 1
start_image=0;%these are the images corresponding to each z slice in a single frame (NOTE: z-slices start with index 0)
end_image=269;%last z slice (top)
x1=1; %these are limits of the region of interest
x2=797;
y1=14;
y2=800; 
rlist = 12.5; % Bead radius in pixels; should be set to a whole number or +.5
Nlist = 1e4; % List of how many grains to expect for each size in a single frame
cols = 12;
c1=clock;
%%
for j = rlist %radius -- this loop can be used for bidisperse systems
    N_guess = Nlist(rlist == j);
    if  ~testscan
        row = 1;
        if cycles > 0
            largeparts = zeros(N_guess*cycles*frames,cols+2);%MAY NEED TO CHANGE DIMENSIONS DEPENDING ON HOW LARGE tempresultlarge IS
        else
            largeparts = zeros(N_guess*frames,cols);
        end
    end
    h = ceil(j/4); % Cavity radius in pixels (okay to overestimate)
    % Create shell indicies list aka "local_ind"
    [X,Y,Z] = meshgrid(1:2*j+1,1:2*j+1,1:2*j+1);
    x = X(1,:,1);
    y = transpose(Y(:,1,1));
    z = Z(1,1,:);
    z = transpose(z(:));
    local_list = transpose(combvec(y,x,z));
    [~,~,R]=cart2sph((local_list(:,2)-j-1),(local_list(:,1)-j-1),(local_list(:,3)-j-1));
    local_sph=local_list((R(:,1)>j),:);
    local_list=local_list((R(:,1)>4/5*j)|(R(:,1)<2/3*j),:);%local list is now a (inverse of) shell centered at current centroid
    clear R;
    local_ind=sub2ind(size(X),local_list(:,1),local_list(:,2),local_list(:,3));%indicies of a shell in a (2*r+1)x(2*r+1)x(2*r+1) box
    local_sph_IND=sub2ind(size(X),local_sph(:,1),local_sph(:,2),local_sph(:,3));%indices of a core in a (2*r+1)x(2*r+1)x(2*r+1) box
    clear local_list X Y Z x y z;
    radius=j;
    sigma0=1;
    for k=0:cycles%cycle -- this loop is only used for cyclic experiments
        for i=1:frames %frames -- frame number within a cycle (cyclic shear) or overall frame number (steady shear)
            disp('====================');
            disp('MAIN LOOP');
            disp('Starting extraction for all particles:');
            cml1=clock;
            if ~testscan
                imageprefix=['shear' num2str(i,'%03.0f')];
                if cycles > 0
                    imagefolder=[basefolder 'cycle' num2str(k,'%03.0f') '\frame' num2str(i,'%03.0f') '\'];
                    disp(['Currently j = ' num2str(j) ', k = ' num2str(k) ', i = ' num2str(i) '.']);
                else
                    imagefolder=[basefolder 'frame' num2str(i,'%03.0f') '\'];
                    disp(['Currently j = ' num2str(j) ', i = ' num2str(i) '.']);
                end
            else
                imagefolder = basefolder;
                imageprefix='test';
            end
            result=analyze_scan_orientations(imagefolder,imageprefix,start_image,end_image,x1,x2,y1,y2,radius,sigma0,AR_z,local_ind,local_sph_IND,h);
            % yields position list in current frame
            disp(['>>Found ' num2str(size(result,1)) ' particles.']);
            if ~testscan
                disp('Saving lifeguard txt-files');
                if cycles > 0
                    tempresultlarge=[result k*ones(size(result,1),1) i*ones(size(result,1),1)];
                    dlmwrite([savefolder 'diameter' num2str(fix(2*j),'%03.0f') 'k' num2str(k,'%03.0f') 'i' num2str(i,'%04.0f') '_large.txt'],tempresultlarge,'delimiter','\t');
                else
                    tempresultlarge=[result i*ones(size(result,1),1)];
                    dlmwrite([savefolder 'diameter' num2str(fix(2*j),'%03.0f') 'i' num2str(i,'%04.0f') '_large.txt'],tempresultlarge,'delimiter','\t');
                end
                largeparts(row:row+size(tempresultlarge,1)-1,1:cols+1*(cycles>0)) = tempresultlarge;
                row = row + size(tempresultlarge,1);
                clear result tempresultlarge
            end
            disp('>>Please don`t turn me off!');
            cml2=clock;
            disp(['>>Main loop cycle took ' num2str(etime(cml2,cml1)) ' seconds.']);          
        end %for i
    end %for k
    if ~testscan
        largeparts(row:end,:) = [];
        if cycles > 0
            largeparts(:,end) = frames*largeparts(:,end-2)+largeparts(:,end-1);
        end
        save([savefolder 'results_diameter' num2str(fix(2*j),'%03.0f') '.mat'],'largeparts');
        dlmwrite([savefolder 'results_large_diameter' num2str(fix(2*j),'%03.0f') '.txt'],largeparts,'delimiter','\t');
    end
    c2=clock;
    disp('=========================================');
    disp(['Program has finished in: ' num2str(etime(c2,c1)) ' seconds.']);
end %for j
clear AR_z local_ind local_sph_IND c1 c2 cml1 cml2 i j k sigma0;