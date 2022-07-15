clear all;

global globalIxedge globalVal globalIid globalEdges globalNedge

x = [1.68556, 2.9367, 3.2575, 3.5056, 3.7024, 3.874, 4.001, 4.158, 4.294, 4.38, 4.482, 4.574, 4.71, 4.826, 4.95, 5.081, 5.25, 5.481, 5.717, 6.188, 7.676];
y = x;
z = [0.1475654 , 0.16199998, 0.17503536, 0.18319088, 0.18684941, 0.19110608, 0.19381687, 0.19729595, 0.20060083, 0.20317237, 0.20534391, 0.20682494, 0.20821091, 0.20966452, 0.2107259 ,0.21157703, 0.21215075, 0.21251113, 0.21285987, 0.2130213 , 0.25030509];
d = variogram(x', z', 'plotit', true, 'nrbins', 10)

for int i = 1:size(globalVal,1)
    for int j = 1:size(globalVal,1)
        varMat[i,j] = 


% scatter(globalIid(:,2), globalIid(:,1))
% x = rand(100,1)*4-2;  
% y = rand(100,1)*4-2;
% z = 3*sin(x*15)+ randn(size(x));
% d = variogram(x,z, 'plotit', true)