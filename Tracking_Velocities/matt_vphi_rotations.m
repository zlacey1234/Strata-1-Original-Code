%mh_vphi
%This program calculates the velocity in the sheer direction (phi) for all
%particles.
%changed to scs_vphi with minor alterations
%changed to matt_vphi with other alterations
%changed to matt_vphi_rotations with inclusion of calculations of
%rotation axis and angular displacement
function tracks=matt_vphi_rotations(trackfile,cenx,ceny,time_col) 
track_raw=trackfile;
max_row=size(track_raw,1);
max_id=track_raw(max_row,time_col+1);
frames = max(track_raw(:,time_col));
tracks=zeros(frames,max_id,10);
for j=1:max_row
   time=track_raw(j,time_col);
   id=track_raw(j,time_col+1);
   [tracks(time,id,1),tracks(time,id,2),tracks(time,id,3)]=cart2pol(track_raw(j,1)-cenx,track_raw(j,2)-ceny,track_raw(j,3));
end
%now, calculate the displacement vectors (translational -- cylindrical;
%rotational -- spherical)
dt = 3;
for time=1:frames-dt
    for id=1:max_id
        disp(['t = ' num2str(time) ', id = ' num2str(id)]);
        %find the t1 and t2 positions of the particle
        th1=tracks(time,id,1);
        th2=tracks(time+dt,id,1);
        r1=tracks(time,id,2);
        r2=tracks(time+dt,id,2);
        z1=tracks(time,id,3);
        z2=tracks(time+dt,id,3);
        %make sure that neither are '0,0,0'
        if (th1==0.0 && r1==0 && z1==0) || (th2==0 && r2==0 && z2==0)
            tracks(time,id,4)=nan;
            tracks(time,id,5)=nan;
            tracks(time,id,6)=nan;
        else
            delth=th2-th1;
            if delth < -3*pi/2
                delth = delth + 2*pi;
            elseif delth > 3*pi/2
                delth = delth - 2*pi;
            end
            delr=r2-r1;
            delz=z2-z1;
            tracks(time,id,4)=delth/dt;
            tracks(time,id,5)=delr/dt;
            tracks(time,id,6)=delz/dt;
        end
        % make sure that these reference the proper columns
        sigma1_init = track_raw(track_raw(:,time_col)==time & track_raw(:,time_col+1)==id,time_col-6:time_col-4);
        sigma1_final = track_raw(track_raw(:,time_col)==time+dt & track_raw(:,time_col+1)==id,time_col-6:time_col-4);
        sigma2_init = track_raw(track_raw(:,time_col)==time & track_raw(:,time_col+1)==id,time_col-3:time_col-1);
        sigma2_final = track_raw(track_raw(:,time_col)==time+dt & track_raw(:,time_col+1)==id,time_col-3:time_col-1);
        if isempty(sigma1_init) || isempty(sigma1_final) || isempty(sigma2_init) || isempty(sigma2_final)
            tracks(time,id,7) = nan;
            tracks(time,id,8) = nan;
            tracks(time,id,9) = nan;
            tracks(time,id,10) = nan;
        else
            [axis,ang_disp,~] = FindRotationAxis(sigma1_init,sigma1_final,sigma2_init,sigma2_final);
            tracks(time,id,7) = axis(1);
            tracks(time,id,8) = axis(2);
            tracks(time,id,9) = axis(3);
            tracks(time,id,10) = ang_disp;
        end
    end
end