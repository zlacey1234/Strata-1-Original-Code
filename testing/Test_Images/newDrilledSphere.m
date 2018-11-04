function sphtest = newDrilledSphere(r,x,theta,phi,val_sol,val_liq,no_holes)
% r > x

%reset(RandStream.getDefaultStream,sum(100*clock))
sphtest = val_sol*ones(2*r+1,2*r+1,2*r+1);
cent_xy = (size(sphtest,1)+1)/2;
cent_z = (size(sphtest,3)+1)/2;
[X,Y,Z] = meshgrid(1:size(sphtest,2),1:size(sphtest,1),1:size(sphtest,3));
Xnew = (X-cent_xy).*cos(theta(1))             -(Y-cent_xy).*sin(theta(1));
Ynew = (X-cent_xy).*cos(phi(1)).*sin(theta(1))+(Y-cent_xy).*cos(phi(1)).*cos(theta(1))-(Z-cent_z).*sin(phi(1));
Znew = (X-cent_xy).*sin(phi(1)).*sin(theta(1))+(Y-cent_xy).*sin(phi(1)).*cos(theta(1))+(Z-cent_z).*cos(phi(1));           
rho_test = hypot(Xnew,Ynew);
R_test = hypot(rho_test,Znew);
sphtest(R_test > r) = val_liq;
if no_holes
    sphtest(rho_test <= x & abs(Znew) <= r & R_test <= r) = val_liq;
end

if no_holes == 2
    Xnew2 = (X-cent_xy).*cos(theta(2))             -(Y-cent_xy).*sin(theta(2));
    Ynew2 = (X-cent_xy).*cos(phi(2)).*sin(theta(2))+(Y-cent_xy).*cos(phi(2)).*cos(theta(2))-(Z-cent_z).*sin(phi(2));
    Znew2 = (X-cent_xy).*sin(phi(2)).*sin(theta(2))+(Y-cent_xy).*sin(phi(2)).*cos(theta(2))+(Z-cent_z).*cos(phi(2));           
    rho_test2 = hypot(Xnew2,Ynew2);
    R_test2 = hypot(rho_test2,Znew2);
    sphtest(rho_test2 <= x & abs(Znew) <= r & R_test <= r) = val_liq;
elseif no_holes == 1 || no_holes == 0
    return;
else
    error('Improper number of sphere holes!');
end
