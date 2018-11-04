figure;
[x,y,z] = sphere;
surf(x,y,z,'FaceColor','none');
axis equal
axis off
hold on
grain_id = 301;
orientations = tracks(tracks(:,end)==grain_id,4:9); % Maybe limit this to not the full range of frames
spline_pts = 1; % number of points to interpolate between orientation measurements
ori_plot = zeros(size(orientations,1)*(1+spline_pts)-spline_pts,6);
for n = 1:size(orientations,1)-1
    row = 1+(n-1)*(spline_pts+1);
    ori_plot(row,:) = orientations(n,:);
    axis_1 = -cross(orientations(n,1:3),orientations(n+1,1:3))/norm(cross(orientations(n,1:3),orientations(n+1,1:3)));
    axis_2 = -cross(orientations(n,4:6),orientations(n+1,4:6))/norm(cross(orientations(n,4:6),orientations(n+1,4:6)));
    theta_1 = acos(dot(orientations(n,1:3),orientations(n+1,1:3)))/(spline_pts+1);
    theta_2 = acos(dot(orientations(n,4:6),orientations(n+1,4:6)))/(spline_pts+1);
    for m = 1:spline_pts
        ori_plot(row+m,1:3) = Rotation3D(ori_plot(row+m-1,1:3)',axis_1,theta_1)';
        ori_plot(row+m,4:6) = Rotation3D(ori_plot(row+m-1,4:6)',axis_2,theta_2)';
    end
end
ori_plot(end,:) = orientations(end,:);
leng = size(ori_plot,1);
h1 = cline(ori_plot(:,1),ori_plot(:,2),ori_plot(:,3),leng:-1:1,autumn);
h2 = cline(ori_plot(:,4),ori_plot(:,5),ori_plot(:,6),leng:-1:1,winter);
hold off