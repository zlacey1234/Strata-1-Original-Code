clear
close all
theta = 2*pi*rand-pi;
phi = pi*rand/2;
window_size = 31;
half = (window_size-1)/2;
cent = [16 16 16];
r = 15;
x = 3.75;
%part_id = n;
%cent = round(out_as(part_id,1:3));
%rod  = out_as(part_id,4:6);
%figure;
%[X,Y,Z] = meshgrid(cent(2)-half:1:cent(2)+half,cent(1)-half:1:cent(1)+half,cent(3)-half:1:cent(3)+half);
[X,Y,Z] = meshgrid(1:31,1:31,1:31);
xlist = squeeze(X(1,:,1));
ylist = squeeze(transpose(Y(:,1,1)));
Ztemp = Z(1,1,:);
zlist = squeeze(transpose(Ztemp(:)));
val_sol = 0;
val_liq = 1;
sphtest = newDrilledSphere(r,x,theta,phi,val_sol,val_liq);
sampleIND = sub2ind(size(sphtest),Y,X,Z);


intense = sphtest(sampleIND);
scatter3(X(:),Y(:),Z(:),1,intense(:));
%axis equal
%axis off
%hold on
%plot3(cent(1)-half*rod(1):rod(1):cent(1)+half*rod(1),cent(2)-half*rod(2):rod(2):cent(2)+half*rod(2),cent(3)-half*rod(3):rod(3):cent(3)+half*rod(3),'g.','MarkerSize',30);
%plot3(cent(1),cent(2),cent(3),'ko','MarkerSize',15,'LineWidth',5);
%view(51,-24);