%p_dt = zeros(1,301);
%list_new = [];
%num_lost = 0;
%for t_ind = 1:298
for t_ind = 1
    %t_ind
    t_ref = t_ind+55;
    ref_links = list(list(:,1,t_ind)~=0,1:2,t_ind);
    persist_links = ref_links;
    num_ref = size(ref_links,1);
    list_new = [list_new; ref_links];
    for dt = 4:(302-t_ind)
        tracks_tdt = tracks(tracks(:,14)==t_ref+dt,:);
        current_links = list(list(:,1,t_ind+dt)~=0,1:2,t_ind+dt);
        persist_links = intersect(persist_links,current_links,'rows');
        same_links = intersect(ref_links,current_links,'rows');
        %for i = 1:size(ref_links,1)
        %    if ~ismember(ref_links(i,:),same_links,'rows') && ismember(ref_links(i,1),tracks_tdt(:,15)) && ismember(ref_links(i,2),tracks_tdt(:,15))
        %        num_lost = num_lost + 1;
        %    end
        %end
        %p_dt(dt) = p_dt(dt) + num_lost;
        %num_lost = 0;
    end
end

clear i tracks_tdt