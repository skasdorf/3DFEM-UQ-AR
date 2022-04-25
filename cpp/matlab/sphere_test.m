close all
%elist = [9,8,3,10];
%elist = [8 10 3 15 374 45 9 649 650 651 652 653 654];
%elistorig = [9 8 10 3 15 374 45];
%elist = [9 8 10 3 15 374 45 14]
%elist = [9, 652];
elist =  [9 653];
figure;
for edex = 1:length(elist);
allpoints = dlmread('points_ref.txt');
element1 = dlmread('element_ref.txt');
element1 = element1(:,3:end);
element1 = element1(elist(edex),:);
points = [];
for i=1:length(element1)
    index = element1(i);
    points(end+1,:) = allpoints(index, 2:4);
end

% points = [points; .3*points];
for i=1:length(points)
    if (edex == 1)
    scatter3(points(i,1), points(i,2), points(i,3), 'r');
    text(points(i,1), points(i,2), points(i,3), num2str(i), 'FontSize', 14);
    else
     scatter3(points(i,1), points(i,2), points(i,3), 'g');
     text(points(i,1), points(i,2), points(i,3), "         " + num2str(i), 'FontSize', 10);
    end
    
    hold on
end
end
hold off