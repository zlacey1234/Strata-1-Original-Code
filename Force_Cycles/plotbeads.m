function plotbeads(config,A)

figure(2)

for i = 1:length(config(:,1))
    x = config(i,1);
    y = config(i,2);
    r = config(i,3)*config(i,4)/2;
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp);
    hold on
end

numparts = length(A(:,1));

for i = 1:numparts-1
    for j = (i+1):numparts
        if A(i,j) == 1 && abs(config(i,1) - config(j,1)) < 48 && abs(config(i,2) - config(j,2)) < 48
            plot([config(i,1) config(j,1)],[config(i,2) config(j,2)],'b-','LineWidth',2)
            hold on
        end
    end
end

plot(config(:,1),config(:,2),'.')
%axis([-10 10 -10 10])
axis equal