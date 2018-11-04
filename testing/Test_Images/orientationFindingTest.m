clear
clc
close all
%reset(RandStream.getDefaultStream,sum(100*clock))
% Parameters of drilled sphere
R = 12; % radius of sphere
x = 0.25; % size ratio of drill hole radius to sphere radius
N = 101; % Size of cubic cell within which we will try to extract sphere
cent = (N+1)/2; % x, y, and z coordinate of centroid
%halfcyl = (H-1)/2;
% Average value and noise of "images"
valin = 1000;
valout = 2000;
noise_in = 0;
noise_out = 0;
n = 1; % Number of spheres to test
% Parameters for analyze_scan
sigma0 = 1;
AR_z = 1;
x1 = 1;
x2 = N;
y1 = 1;
y2 = N;

no_holes = 2;
out_as = zeros(n,5+3*no_holes);
theta = zeros(n,no_holes);
phi = theta;
theta_m = theta;
phi_m = theta;
axes = zeros(n,3,no_holes);
diffs = zeros(no_holes*n,1);
tic
for j = 1:n
    disp(['i = ' num2str(j)]);
    phi(j,1) = pi*rand/2;
    theta(j,1) = 2*pi*rand-pi;
    %phi = [pi/2 0];
    %theta = [pi/2 0];
    
    [axes(j,2,1),axes(j,1,1),axes(j,3,1)] = sph2cart(theta(j,1),pi/2-phi(j,1),1);


     if no_holes == 2
         theta_comp = rand;
         phi_comp = sqrt(1-theta_comp^2);
%     
         axes(j,:,2) = theta_comp*[cos(theta(j,1)) -sin(theta(j,1)) 0] - phi_comp*[cos(phi(j,1))*sin(theta(j,1)) cos(phi(j,1))*cos(theta(j,1)) -sin(phi(j,1))];

         if sum(axes(j,:,1).*axes(j,:,2)) > 1e-15
             pause;
         end
         [theta(j,2),phi(j,2),~] = cart2sph(axes(j,2,2),axes(j,1,2),axes(j,3,2));
         phi(j,2) = pi/2-phi(j,2);
     end
    
    sphtest = newDrilledSphere(R,x*R,theta(j,:),phi(j,:),valin,valout,no_holes);

    
    
    sample = valout*ones(N,N,N);
    sample = sample + noise_out*randn(size(sample));
    sample((N+1)/2-R:(N+1)/2+R,(N+1)/2-R:(N+1)/2+R,(N+1)/2-R:(N+1)/2+R) = sphtest;
    
    
    [X,Y,Z] = meshgrid(1:2*R+1,1:2*R+1,1:2*R+1);
    xcoords = X(1,:,1);
    ycoords = transpose(Y(:,1,1));
    zcoords = Z(1,1,:);
    zcoords = transpose(zcoords(:));
    %sampleIND = sub2ind(size(sphtest),Y,X,Z);
    %intense = sphtest(sampleIND);
    %scatter3(X(:),Y(:),Z(:),1,intense(:));
    %axis equal
    %axis off
    clear sphtest;
    local_list = transpose(combvec(ycoords,xcoords,zcoords));
    [~,~,Rtest]=cart2sph((local_list(:,2)-R-1),(local_list(:,1)-R-1),(local_list(:,3)-R-1));
    h = ceil(x*R);
    if no_holes == 2
        local_sph=local_list((Rtest(:,1)>R),:);
        local_list=local_list((Rtest(:,1)>4/5*R)|(Rtest(:,1)<2/3*R),:);%local list is now a (inverse of) shell centered at current centroid
        clear Rtest;
        local_ind=sub2ind(size(X),local_list(:,1),local_list(:,2),local_list(:,3));
        local_sph_IND=sub2ind(size(X),local_sph(:,1),local_sph(:,2),local_sph(:,3));
        clear local_list X Y Z xcoords ycoords zcoords;
        sph_m = analyze_scan_TESTorientations(sample,x1,x2,y1,y2,R,sigma0,AR_z,local_ind,local_sph_IND,h);
    elseif no_holes == 1
        local_sph = local_list((Rtest(:,1)>4/5*R),:);
        local_sph_IND = sub2ind(size(X),local_sph(:,1),local_sph(:,2),local_sph(:,3)); 
        clear Rtest;
        clear local_list X Y Z xcoords ycoords zcoords;
        sph_m = analyze_scan_TESTorientation(sample,x1,x2,y1,y2,R,sigma0,AR_z,local_sph_IND,h);
    else
        sph_m = analyze_scan_TEST(sample,x1,x2,y1,y2,R,sigma0,AR_z);
    end
    if ~isempty(sph_m)
        out_as(j,:) = sph_m;
        if no_holes > 0
            [theta_m(j,1),phi_m(j,1),~] = cart2sph(out_as(j,7),out_as(j,6),out_as(j,8));
            phi_m(j,1) = pi/2-phi_m(j,1);
            %figure(1)
            %hold on
            %plot3([R-R*out_as(j,6) R+R*out_as(j,6)],[R-R*out_as(j,7)  R+R*out_as(j,7)],[ R-R*out_as(j,8)  R+R*out_as(j,8)],'r-','LineWidth',8);
            if no_holes == 2
                [theta_m(j,2),phi_m(j,2),~] = cart2sph(out_as(j,10),out_as(j,9),out_as(j,11));
                phi_m(j,2) = pi/2-phi_m(j,2);
                %plot3([R-R*out_as(j,9) R+R*out_as(j,9)],[R-R*out_as(j,10) R+R*out_as(j,10)],[R-R*out_as(j,11) R+R*out_as(j,11)],'g-','LineWidth',8);
            end
        end
   
        %theta_m(j) = theta_1(j);
        %phi_m_1(j) = phi_1(j);
        %theta_m_2(j) = theta_2(j);
        %phi_m_2(j) = phi_2(j);
    end
    if no_holes == 2
        diffs_1 = sqrt(sum(([out_as(j,9:11);out_as(j,6:8)]-repmat(axes(j,:,1),2,1)).^2,2));
        [diffs(2*j-1),i_min] = min(diffs_1);
        diffs(2*j) = sqrt(sum((out_as(j,3+3*i_min:5+3*i_min)-axes(j,:,2)).^2,2));
    elseif no_holes == 1
        diffs(j) = sqrt(sum((out_as(j,6:8)-axes).^2,2));
    end
end

%d_r = out_as(:,1:3) - cent*ones(n,3);
%d_th = theta_m - theta;
%d_ph = phi_m - phi;
clear x1 x2 y1 y2 sigma0 AR_z
toc
% figure(1)
% hist(d_r(d_r(:,1)~=0,1));
% figure(2)
% hist(d_r(d_r(:,2)~=0,2));
% figure(3)
% hist(d_r(d_r(:,3)~=0,3));
% figure(4)
% hist(sqrt(sum(d_r(d_r(:,1)~=0,1:3).^2,2)))
% figure(5)
% hist(d_th(d_th~=0))
% figure(6)
% hist(d_ph(d_ph~=0))
