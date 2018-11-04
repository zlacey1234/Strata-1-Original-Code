
%% create id structure format for easier computation
clear strc
d=1;
for id=min(tracks(:,15)):max(tracks(:,15))
    if isempty(tracks(tracks(:,15)==id,15)==1)
    continue
    end
strc{d,1}=tracks(tracks(:,15)==id,:);
d=d+1;
end
%% flip vectors that flipped about 180deg between frames
for d=1:length(strc)
    for f=2:length(strc{d,1}(:,1))
        if norm(strc{d,1}(f,9:11)+strc{d,1}(f-1,9:11))<1
            strc{d,1}(f,9:11)=-strc{d,1}(f,9:11);
        end
    end
end
%% smoothing (cant smooth the whole collection of beads!)
for d=1:length(strc)
    if length(strc{d,1}(:,1))>=3
    strc{d,1}(:,9)=fastsmooth(strc{d,1}(:,9),3,1,1);
    strc{d,1}(:,10)=fastsmooth(strc{d,1}(:,10),3,1,1);
    strc{d,1}(:,11)=fastsmooth(strc{d,1}(:,11),3,1,1);
    end
end

%% find projections onto plane of rotation
for d=1:length(strc)
    strc{d,2}=zeros(length(strc{d,1}(:,1)),3);
    for f=2:length(strc{d,1}(:,1))-1
        strc{d,2}(f,1:3)=cross((strc{d,1}(f,9:11)-strc{d,1}(f-1,9:11)),(strc{d,1}(f+1,9:11)-strc{d,1}(f,9:11)))... %creates normalized vector along axis of rot
            /norm(cross((strc{d,1}(f,9:11)-strc{d,1}(f-1,9:11)),(strc{d,1}(f+1,9:11)-strc{d,1}(f,9:11))));%cross(2-1,3-2)/norm of cross where 1 is the frame before frame f(current frame)
%find magnitude of w
        proj3=strc{d,1}(f+1,9:11)-dot(strc{d,1}(f+1,9:11),strc{d,2}(f,1:3)).*strc{d,2}(f,1:3);%orientation vector - component along omega (=component in plane perp to w)
        proj1=strc{d,1}(f-1,9:11)-dot(strc{d,1}(f-1,9:11),strc{d,2}(f,1:3)).*strc{d,2}(f,1:3);
        strc{d,2}(f,1:3)=acos(dot(proj1,proj3)/(norm(proj1)*norm(proj3)))*strc{d,2}(f,1:3);%=the angle between the projections*omega. ie the angle subtended is ~the magnitude of omega
    end
end
%% find z component of rotation
for d=1:length(strc)%cycle through ID
   %strc{d,2}=zeros(length(strc{d,1}(:,1)),1);
    for f=2:length(strc{d,1}(:,1))   %cycle through frame  
%find magnitude of w
       cr=cross([strc{d,1}(f,9:10) 0] ,[(strc{d,1}(f,9:10)-strc{d,1}(f-1,9:10)) 0]);% cross between OR vector and dOR vector between frames
       if cr(3)>0
        strc{d,1}(f,16)=-2*(pi/180)+1*abs(acos(dot(strc{d,1}(f-1,9:10),strc{d,1}(f,9:10))/(norm(strc{d,1}(f-1,9:10))*norm(strc{d,1}(f,9:10)))));%rotation relative to solid body of disk
       else 
        strc{d,1}(f,16)=-2*(pi/180) +-1*abs(acos(dot(strc{d,1}(f-1,9:10),strc{d,1}(f,9:10))/(norm(strc{d,1}(f-1,9:10))*norm(strc{d,1}(f,9:10)))));%the z component of omega_bead
       end%could we potentially just use w=R x V in place of the acos??
    end
end
%% convert smoothed tracks in strc back to a long matrix
trackfile=[];
for d=1:length(strc)
    trackfile=[trackfile; strc{d,1}];
end
%% plots
cm=colormap(jet(200));
hold on
for ix=1:373%1:100:length(strd)
    co=rand(1,3);
plot3(strc{ix,1}(:,1),strc{ix,1}(:,2),strc{ix,1}(:,3),'-','Color',co)%plot tracks in 3d
 for f=2:length(strc{ix,1}(:,1)) 
    plot3([strc{ix,1}(f,1),strc{ix,1}(f,1)+strc{ix,2}(f,4)*7],[strc{ix,1}(f,2),strc{ix,1}(f,2)+strc{ix,2}(f,5)*7],...
        [strc{ix,1}(f,3),strc{ix,1}(f,3)+strc{ix,2}(f,6)*7],'k-');
 end
end