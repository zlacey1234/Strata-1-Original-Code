%PARAMETERS -- LOOK THESE UP from data
tfirst=2;
tlast=239;
rmin=0;
rbin=24.5; % Grain diameter
rmax=500;
zmin=40;
zbin=24.5; % Grain diameter
zmax=208;
%cenx and ceny should already be set

% Initiate
max_id=size(tracks_v,2);
zcount=1;
r=rmin;
z=zmax;
data_all_th=[];
data_all_ph=[];
data_out_th=[];
data_out_ph=[];
fig_start = 4;
f1=figure(fig_start); hold on;
f2=figure(fig_start+1); hold on;
f3=figure(fig_start+2); hold on;
f4=figure(fig_start+3); hold on;

while zit<zmax
    r=0;
    dataz_th=[];
    dataz_ph=[];
    while r<rmax
        dth_list=[];
        dph_list=[];
        for t=tfirst:tlast
            disp(['t = ' num2str(t) ', z = ' num2str(z) ', r = ' num2str(r)]);
            tracks_t = reshape(tracks_v(t,:,:),max_id,8);
            dth=transpose(tracks_t(tracks_t(:,2)>=r & tracks_t(:,2)<r+rbin & tracks_t(:,3)<=z & tracks_t(:,3)>z-zbin,7)); % dth is now rotational motion of grain in plane of disk, not to be confused with dth for azimuthal translational motion
            dth=dth(~isnan(dth));
            dth_list=[dth_list dth];
            dph=transpose(tracks_t(tracks_t(:,2)>=r & tracks_t(:,2)<r+rbin & tracks_t(:,3)<=z & tracks_t(:,3)>z-zbin,8)); 
            dph=dph(~isnan(dph));
            dph_list=[dph_list dph];
        end
        %calculate the mean of the list and stdev
        dataz_th=[dataz_th; [r+rbin./2.0 zit+zbin./2.0 nanmean(dth_list) std(dth_list)./sqrt(size(dth_list,2))]];
        dataz_ph=[dataz_ph; [r+rbin./2.0 zit+zbin./2.0 nanmean(dph_list) std(dph_list)./sqrt(size(dph_list,2))]];
        r=r+rbin;
    end
    %calculate the r-graident of the velocity
    dataz_th(:,5)=gradient(dataz_th(:,3),rbin);
    dataz_ph(:,5)=gradient(dataz_large(:,3),rbin);
    data_all_th(zcount,:,:,:,:)=dataz_th;    
    data_out_th=[data_out_th; dataz_th];
    data_all_ph(zcount,:,:,:,:)=dataz_ph;    
    data_out_ph=[data_out_ph; dataz_ph];

    zit=zit+zbin;
    zcount=zcount+1;
    figure(f1);errorbar(dataz_th(:,1),dataz_th(:,3),dataz_th(:,4),'-bo');
    figure(f2);plot(dataz_th(:,1),dataz_th(:,5).*dataz_th(:,1),'-bo');
    figure(f3);errorbar(dataz_ph(:,1),dataz_ph(:,3),dataz_ph(:,4),'-ro');
    figure(f4);plot(dataz_ph(:,1),dataz_ph(:,5).*dataz_ph(:,1),'-ro');
end

rminmesh=min(data_out_th(:,1));
rmaxmesh=max(data_out_th(:,1));
zminmesh=min(data_out_th(:,2));
zmaxmesh=max(data_out_th(:,2));
[mr,mz]=meshgrid(rminmesh:rbin:rmaxmesh,zminmesh:zbin:zmaxmesh);

p_th=mr;
p_ph=mr;
for i=1:length(data_out_th(:,1))
    rval=data_out_th(i,1);
    zval=data_out_th(i,2);
    dthval=data_out_th(i,3);
    dphval=data_out_ph(i,3);
    rr=find(mr(1,:)==rval);
    rz=find(mz(:,1)==zval);
    p_th(rz,rr)=dthval;
    p_ph(rz,rr)=dphval;
end

[FR_th,FZ_th]=gradient(p_th,rbin,zbin);
[FR_ph,FZ_ph]=gradient(p_ph,rbin,zbin);

r_dwdr=mr.*FR_th;
r_dwdz=mr.*FZ_th;
r_dgdr=mr.*FR_ph;
r_dgdz=mr.*FZ_ph;
