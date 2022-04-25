clear all
close all
% points = [-1 -1 -1; 1 -1 -1; -1 1 -1; 1 1 -1; -1 -1 1; 1 -1 1; -1 1 1; 1 1 1];
% points_center = .3*points;
% points = [points; points_center];
% 
% x = -1:1:1;
% y = -1:1:1;
% z = -1:1:1;
points = [-1 -1 -1; 0 -1 -1; 1 -1 -1; -1 0 -1; 0 0 -1; 1 0 -1; -1 1 -1; 0 1 -1; 1 1 -1];
z = [0 0 1];
points1 = points + z;
points2 = points + 2*z;
points = [points; points1; points2];
points = [points; .3*points];
dif = 27;
points(41,:) = [];
for i=1:27
    if (i ~= 14)
    points(end+1,:) = halfway(points(i, :), points(i+dif, :));
    else
        dif = dif - 1;
    end
end
figure;
for i=1:length(points)
    if (i < 28)
    scatter3(points(i,1), points(i,2), points(i,3), 'r');
    elseif (i < 54)
    scatter3(points(i,1), points(i,2), points(i,3), 'b');
    else 
    scatter3(points(i,1), points(i,2), points(i,3), 'g');
    end
    text(points(i,1), points(i,2), points(i,3), num2str(i));
    hold on
end
hold off
xlabel('x');
ylabel('y');
zlabel('z');
figure;
for i=1:length(points)
    if (points(i,3) > 0)
       scatter3(points(i,1), points(i,2), points(i,3), 'r');
       text(points(i,1), points(i,2), points(i,3), num2str(i));
    end
    hold on
end
hold off
xlabel('x');
ylabel('y');
zlabel('z');