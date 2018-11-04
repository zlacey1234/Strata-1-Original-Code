id_col = size(tracks,2);
time_col = id_col-1;
folder = 'C:\students\matt\56on tracks\contacts\';
min_time = 56;
max_time = 357;
cosines = zeros(max_time-min_time+1,17000);
c_bin = -0.995:0.01:0.995;
c_dist = c_bin;
bias = zeros(1,max_time-min_time+1);
flat = bias;
t_ind = 0;
for ref_time = min_time:max_time;
    count = 0;
    ref_time
    t_ind = t_ind + 1;
    links = load([folder 'contacts_time' num2str(ref_time,'%04.0f') '.dat']);
    links = links(:,1:2);
    tracks_t = tracks(tracks(:,time_col)==ref_time,:);
    index = 0;
    for i = 1:size(links,1)
        ref_id = links(i,2);
        ref_pos = tracks_t(tracks_t(:,id_col)==ref_id,1:3);
        ref_ori = tracks_t(tracks_t(:,id_col)==ref_id,9:11);
        part_id = links(i,1);
        part_pos = tracks_t(tracks_t(:,id_col)==part_id,1:3);
        rel_pos = part_pos - ref_pos;
        index = index + 1;
        cosines(t_ind,index) = dot(rel_pos,ref_ori)/norm(rel_pos);
        if abs(cosines(t_ind,index)) > 0.97 && ref_id < part_id
            count = count + 1;
            list(count,:,t_ind) = [ref_id part_id];
        end
    end
    c_dist(t_ind,:) = hist(cosines(t_ind,find(cosines(t_ind,:))),c_bin);
    c_dist(t_ind,:) = c_dist(t_ind,:)/trapz(c_bin,c_dist(t_ind,:));
    bias(t_ind) = sum(0.01*c_dist(t_ind,1:3))+sum(0.01*c_dist(t_ind,end-2:end))-6*0.01*0.5;
    flat(t_ind) = mean(c_dist(4:end-3));
end

clear i id_col index links part_id part_pos ref_id ref_pos ref_ori rel_pos t_ind time_col tracks_t