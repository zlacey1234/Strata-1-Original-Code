function plot_contacts(A,coords)

if A ~= transpose(A)
    disp('Improper A');
end

figure(1)
plot(coords(:,1),coords(:,2),'.')
%plot3(coords(:,1),coords(:,2),coords(:,3),'.')
hold on

numparts = length(A(:,1));

for i = 1:numparts-1
    for j = (i+1):numparts
        if A(i,j) == 1 %&& (((coords(i,3) == -4 || coords(i,3) == -2.5) && coords(i,1) < -2 && coords(i,1) > -4 && coords(i,2) > -4.25 && coords(i,2) < -2.75 ...
                %&& (coords(j,3) == -4 || coords(j,3) == -2.5) && coords(j,1) < -2 && coords(j,1) > -4 && coords(j,2) > -4.25 && coords(j,2) < -2.75) ...
                %|| (coords(i,3) == -3.25 && coords(i,1) <= -2 && coords(i,1) >= -4 && coords(i,2) > -4.5 && coords(i,2) < -2.5 && ...
                %coords(j,3) == -3.25 && coords(j,1) <= -2 && coords(j,1) >= -4 && coords(j,2) > -4.5 && coords(j,2) < -2.5) ...
                %|| ((coords(i,3) == -4 || coords(i,3) == -2.5) && coords(i,1) < -2 && coords(i,1) > -4 && coords(i,2) > -4.25 && coords(i,2) < -2.75 ...
                %    && coords(j,3) == -3.25 && coords(j,1) <= -2 && coords(j,1) >= -4 && coords(j,2) > -4.5 && coords(j,2) < -2.5) ...
                % || (coords(i,3) == -3.25 && coords(i,1) <= -2 && coords(i,1) >= -4 && coords(i,2) > -4.5 && coords(i,2) < -2.5 && ...
                %    (coords(j,3) == -4 || coords(j,3) == -2.5) && coords(j,1) < -2 && coords(j,1) > -4 && coords(j,2) > -4.25 && coords(j,2) < -2.75))
            plot([coords(i,1) coords(j,1)],[coords(i,2) coords(j,2)],'b-')
            %if coords(i,3) == coords(j,3)
            %    plot3([coords(i,1) coords(j,1)],[coords(i,2) coords(j,2)],[coords(i,3) coords(j,3)],'b.-','MarkerSize', 30)
            %    hold on
            %else
            %    plot3([coords(i,1) coords(j,1)],[coords(i,2) coords(j,2)],[coords(i,3) coords(j,3)],'r-')
            %    hold on
            %    plot3([coords(i,1) coords(j,1)],[coords(i,2) coords(j,2)],[coords(i,3) coords(j,3)],'b.','MarkerSize',30)
            %    hold on
            %end
        end
    end
end

axis equal
hold off