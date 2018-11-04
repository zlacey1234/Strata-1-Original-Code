function newtracks = TrackConverter(oldtracks,min_length)

% newtracks returns a tracks file with columns: x, y, z, u, v, w, t,
%   particle id#
% oldtracks should be in the format given by Ouellette's tracking code
% largeparts is the position list; necessary for grabbing data from other
% columns
% min_length is the lower bound for track length
% cyclic = 1 for cyclic shear, 0 for steady shear

clc
newtracks = zeros(1e8,11);
id = 0;
row = 1;
num_tracks = length(oldtracks);

for i = 1:num_tracks;
    disp(['i = ' num2str(i)]);
    track = oldtracks(i);
    if track.len < min_length
        continue;
    end
    id = id + 1;
    xlist = transpose(track.X);
    ylist = transpose(track.Y);
    zlist = transpose(track.Z);
    tlist = transpose(track.T)-ones(size(track.T'));
    idlist = id*ones(track.len,1);
    alist = transpose(track.A);
    blist = transpose(track.B);
    clist = transpose(track.C);
    dlist = transpose(track.D);
    elist = transpose(track.E);
    flist = transpose(track.F);
    newtracks(row:row+track.len-1,:) = [xlist ylist zlist alist blist clist dlist elist flist tlist idlist];
    row = row+track.len;
end

if row < 1e8
    newtracks(row:end,:) = [];
end