function plot_cycles(cycles,persist,A)
%close(figure(2))

% cycles = c3, c4, c5, or c6; repeat function to plot multiple cycle types
% persist = persistant network with coordinates
% A = adjacency matrix

cyc_length = length(cycles(1,:));
figure(1)
hold on
for i = 1:length(cycles(:,1))
    parts = cycles(i,:);
    for j = 1:cyc_length-1
        r1 = persist(parts(j),:);
        %if r1(1) > 0 && r1(2) > 0
            for k = j+1:cyc_length
                if A(parts(j),parts(k)) == 1
                    r2 = persist(parts(k),:);
                    if abs(r1(1) - r2(1)) < 48 && abs(r1(2) - r2(2)) < 48
                        plot([r1(1) r2(1)],[r1(2) r2(2)],'r.-','LineWidth',2);
                        %plot3([r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)],'r.-','MarkerSize',8);
                        hold on
                    end
                end
            end
        %end
    end
end


hold off

axis equal