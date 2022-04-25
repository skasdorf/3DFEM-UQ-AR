
close all
clear all

plotType = 'Scattered'; %
%plotType = 'Total'; %

frequency = 3E8;
lambda = 3E8/frequency;
k = 2*pi/lambda;
a = 4;
Epsr = 2.25;

data2 = dlmread('FieldPlot_1_2182.dat');
%data2 = dlmread('PMLABC_WORKS.dat');
count = 0;
for i = 1:size(data2,1)
   
   if(abs(data2(i,3))  == 0) && (sqrt(data2(i,1).^2 + data2(i,2).^2) < 5)
      count = count + 1;
      data(count,:) = data2(i,:);
      if(strcmp(plotType,'Scattered'))
          %data(count,8) = sqrt(data(count,8).^2 + data(count,9).^2);
          data(count,8) = data(count,8);
      else
          data(count,8) = data(count,8) + real(exp(-1j*k*data(count,1))); %+ is usually correct
      end
   end
end

clear data2;

x = data(:,1);
y = data(:,2);
z = data(:,8);
n=300; %number of subdivisions for plotting purposes

[xi,yi]=meshgrid(linspace(min(x),max(x),n),linspace(min(y),max(y),n));
zi=griddata(x,y,z,xi,yi,'linear');
for i = 1:size(zi,1)
    for j = 1:size(zi,2)
        xtemp = xi(i,j);
        ytemp = yi(i,1);
        %if(strcmp(plotType,'Total'))
           %if(isnan(zi(i,j)))
               %zi(i,j) = 0;
              %zi(i,j) =  real(exp(-1j*k*xtemp));
           %end
        %end
    end
end
%clear x; clear y;
clear z; clear data;
%contourf(xi,yi,zi,25)






%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Plot from FEM output file
figure(1)
% for i = 1:size(zi,1)
%     for j = 1:size(zi,2)
%         if(norm([xi(1,i),yi(j,1)]) >= 1.3)
%            zi(j,i) = NaN; 
%         end
%     end
% end

imagesc(xi(1,:),yi(:,1)',zi)
%xlim([-1.5,1.5]);
%ylim([-1.5,1.5]);
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
%caxis([-1.8,.5]);
%caxis([-.3,.1])
%caxis([-2,2]);
%caxis([-2.5,4])
%Specific to almond
axis equal;
colormap(jet)

radians = linspace(0,2*pi,100);
for i = 1:length(radians)
    circleX(i) = cos(radians(i));
    circleY(i) = sin(radians(i));
end

hold on
plot(4*circleX,4*circleY,'k');
plot(4.3*circleX,4.3*circleY,'k');
plot(4.6*circleX,4.6*circleY,'k');


