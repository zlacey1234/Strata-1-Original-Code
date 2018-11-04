function cyltest = newCylinder(d,H,theta,phi,valin,valout,box)
% H > d
% generates a cylinder of arbitrary orientation -- note the convention used
% to remain consistent with image coordinates (row,column,height)

%reset(RandStream.getDefaultStream,sum(100*clock))
cyltest = valout*ones(box,box,box);
cent_xy = (size(cyltest,1)+1)/2;
cent_z = (size(cyltest,3)+1)/2;
[X,Y,Z] = meshgrid(1:size(cyltest,2),1:size(cyltest,1),1:size(cyltest,3));
Xnew = (X-cent_xy).*cos(theta)          -(Y-cent_xy).*sin(theta);
Ynew = (X-cent_xy).*cos(phi).*sin(theta)+(Y-cent_xy).*cos(phi).*cos(theta)-(Z-cent_z).*sin(phi);
Znew = (X-cent_xy).*sin(phi).*sin(theta)+(Y-cent_xy).*sin(phi).*cos(theta)+(Z-cent_z).*cos(phi);           
rho_test = hypot(Xnew,Ynew);
cyltest(rho_test <= d/2 & abs(Znew) <= H/2) = valin;