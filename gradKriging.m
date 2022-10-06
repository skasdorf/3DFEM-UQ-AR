clear all;

x = linspace(0, 10, 1000);
stdDev = 1.0;
mean = 5.0;

%create random variable with known pdf
y = pdf('Normal', x, mean, stdDev);

%random y derivative
dy = diff(y) ./ diff(x);

%xSamples
xSample = [x(1), x(250), x(500), x(750), x(1000)];
ySample = [y(1), y(250), y(500), y(750), y(1000)];

%setting up bins
maxDist = 10;
nBins = 10;
binTol = maxDist / nBins;




dMat = zeros(size(xSample,2)^2, 4);
index = 1;
for i = 1:size(xSample,2)
    for j = 1:size(xSample,2)
        
        dMat(index, 1) = i;
        dMat(index, 2) = j;
        dMat(index, 3) = abs(xSample(i) - xSample(j));
        dMat(index, 4) = (ySample(i) - ySample(j)).^2;
        index = index + 1;
        
    end
end


edges = linspace(0, maxDist, nBins+1);
edges(end) = maxDist;

binIdx = zeros(size(xSample, 2)^2,1);
binCounts = zeros(size(edges,2)-1,1);


