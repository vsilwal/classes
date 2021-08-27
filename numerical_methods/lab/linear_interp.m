clear all
close all

N = 10;
xi = linspace(-1,1,N);
yi = [3 2 2 1 0 -3 -4 1 2 4]; % for N=10

%plot
plot(xi, yi, 'o','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k')

% %=========================
% y1 = yi(1); y2 = yi(2);
% x1 = xi(1); x2 = xi(2);
% 
% x = linspace(x1,x2,n);
% 
% y = y1 + (x - x1)*(y2 - y1)/(x2 - x1);
% hold on
% plot(x,y)
% %=========================
% y3 = yi(3);
% x3 = xi(3);
% 
% x = linspace(x2,x3,n);
% 
% y = y2 + (x - x2)*(y3 - y2)/(x3 - x2);
% plot(x,y)

hold on
n = 100;
for i=1:N-1
    y0 = yi(i); y1 = yi(i+1);
    x0 = xi(i); x1 = xi(i+1);
    
    x = linspace(x0,x1,n);
    
    y = y0 + (x - x0)*(y1 - y0)/(x1 - x0);
    
    plot(x,y)
end
