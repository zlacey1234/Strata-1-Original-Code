function plot_cycles_proposal(cycles,contacts)
%close(figure(2))

% cycles = c3, c4, c5, or c6; repeat function to plot multiple cycle types
% contacts = persistant network with coordinates

colormap(jet);
cmap = colormap;
cyc_list = cycles(:);
cyc_id = unique(cyc_list);
max_id = max(cyc_list);
mode_id = mode(cyc_list);
id_scale = round(64*hist(cyc_list,max_id) / sum(cyc_list == mode_id));

figure(2)
for i = 1:length(cyc_id)
    r = unique(contacts(contacts(:,1)==cyc_id(i),3:5),'rows');
    c_point = id_scale(cyc_id(i));
    if r(3) > 100 && r(3) < 115 && c_point > 0
        color_vec = cmap(c_point,:);
        plot3(r(1),r(2),r(3),'o','MarkerFaceColor',color_vec,'MarkerEdgeColor',[1,1,1])
        hold on
    end        
end

hold off
axis equal
box on