% This function returns a list of contacts for a single contact cutoff
% length.  In order to be in contact, two particles must be within rCutoff
% of each other both at time t and t+dt.  Enforcing both frames will
% eliminate spurious contacts involving grains that 'fall out'.  

% This function is called several times during the course of
% probeContactLength; it should be run one more time using the chosen
% contact length

% INPUTS: tracks, dt, rCutoff

% OUTPUTS: contacts -- this is the contact list; this will be formatted as
% a Nx1 cell array, where N is the total number of tracks.  Each entry will
% be a Tx1 cell array, where T is the total number of frames.  Each entry
% in *that* array will be a list of Z particle id's that the particle i is
% in contact with at time t.  Z for a single particle can of course change
% over time.  An empty entry would correspond to a particle that is either
% not tracked at that time, or has no contacts found

% cosines -- same format as contacts, this is a corresponding list of the
% cos \alpha 's between the relative position vector and relative velocity
% vector between each pair of contacts; mainly used for probing contact
% lengths

function [contacts,cosines] = getContacts(tracks, rCutoff, dt, tStart, tEnd)

N = max(tracks(:,end));
T = tEnd-dt-tStart+1;
contacts = cell(N,1);
for i = 1:N
    contacts{i} = cell(T,1);
end

cosines = contacts;

for t = 1:T
    fullPos_t = zeros(N,3);
    fullPos_tdt = zeros(N,3);
    
    tracks_t = tracks(tracks(:,end-1)==t,:);
    tracks_tdt = tracks(tracks(:,end-1)==t+dt,:);
    
    fullPos_t(tracks_t(:,end),:) = tracks_t(:,1:3);
    fullPos_tdt(tracks_tdt(:,end),:) = tracks_tdt(:,1:3);
    
    % Use of fullPos is meant to account for particles that are
    % not found in certain frames; for particles that are not found in
    % frame t or t+dt, the position should remain (0,0,0) 
    
    dists_t = squareform(pdist(fullPos_t));
%      dists_t(fullPos_t(:,1) == 0,:) = zeros(size(dists_t,1),1); % set rows and columns pertaining to missing particles to zero
%      dists_t(:,fullPos_t(:,1) == 0) = zeros(1,size(dists_t,2));
    dists_tdt = squareform(pdist(fullPos_tdt));
%     dists_tdt(fullPos_tdt(:,1) == 0,:) = zeros(size(dists_tdt,1),1);
%     dists_tdt(:,fullPos_tdt(:,1) == 0) = zeros(1,size(dists_tdt,2));
    
    [row_t,col_t] = find(dists_t > 0 & dists_t <= rCutoff);
    [row_tdt,col_tdt] = find(dists_tdt > 0 & dists_tdt <= rCutoff);
    contact_t = intersect([row_t col_t],[row_tdt col_tdt],'rows');
    for i = 1:N
        disp('Contact Analysis')
        beadNo=rCutoff
        time=t
        progressmeter=[i, N]
        contacts_i = contact_t(contact_t(:,1) == i,2)';
        contacts{i}{t} = contacts_i;
        cosines_add = zeros(length(contacts_i),1);
        for j_ind = 1:length(contacts_i)
            j = contacts_i(j_ind);
            rij_t = fullPos_t(j,:)-fullPos_t(i,:);
            rij_tdt = fullPos_tdt(j,:)-fullPos_tdt(i,:);
            velij_t = rij_tdt - rij_t;
            cosij = dot(velij_t,rij_t)/(norm(rij_t)*norm(velij_t));
            cosines_add(j_ind) = cosij;
        end
        cosines{i}{t} = cosines_add;
    end
end