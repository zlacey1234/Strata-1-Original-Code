%v_rotate
%This program calculates the particle rotational velocity in the shear direction for all
%particles.

function tracks_v=v_rotate(tracks, cenx, ceny, time_col, orient_col) 

% orient_col: first column that stores orientation information

track_raw=tracks;
max_row=size(track_raw,1);
max_id=max(track_raw(:,time_col+1)); 
max_t=max(track_raw(:,time_col)); 
tracks_v=zeros(max_t,max_id,8); % time , N_part_id , [th r z sigma_x sigma_y sigma_z del_theta del_phi]

for j=1:max_row
   time=track_raw(j,9);
   id=track_raw(j,10);
   [tracks_v(time,id,1),tracks_v(time,id,2),tracks_v(time,id,3)]=cart2pol(track_raw(j,1)-cenx,track_raw(j,2)-ceny,track_raw(j,3));
   tracks_v(time,id,4)=track_raw(j,orient_col);
   tracks_v(time,id,5)=track_raw(j,orient_col+1);
   tracks_v(time,id,6)=track_raw(j,orient_col+2);
end

%now, calculate the rotational displacement vectors, as a first pass we will only use
%two frames difference

for time=2:max_t-1
    for id=1:max_id
        %find the t1 and t2 orientations of the particle
        sigma1_init=tracks_v(time-1,id,4:6);
        sigma1_final=tracks_v(time+1,id,4:6);
        %make sure that neither are '0,0,0'
        if (isnan(sigma1_init(1)) && isnan(sigma1_init(2)) && isnan(sigma1_init(3))) || (isnan(sigma1_final(1)) && isnan(sigma1_final(2)) && isnan(sigma1_final(3)))
            tracks_v(time,id,7)=nan;
            tracks_v(time,id,8)=nan;
        else
            % Note: MATLAB's convention is phi = 0 when vector lies in
            % xy-plane, and positive direction is toward +z-axis; domain:
            % [-pi/2,pi/2]
            [delth,delph,~]=EulerAngle(sigma1_init,sigma1_final);
            tracks_v(time,id,7)=delth./2;
            tracks_v(time,id,8)=delph./2;
        end
    end
end