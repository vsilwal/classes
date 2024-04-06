clear all
close all

fx = @(x) log(1+x);

xvec = -1:0.1:1;

yvec = fx(xvec);

plot(xvec,yvec,'o-')
hold on
%===========

% 1. Create function handle for the polynomial
% obtained using taylor series expansion
fx1 = @(x) x;
fx2 = @(x) x - (x.^2)/2;
fx3 = @(x) x - (x.^2)/2 + (x.^3)/3;
fx4 = @(x) x - (x.^2)/2 + (x.^3)/3 -(x.^4)/4;

% 2. xnew and ynew values obtained using the new 
% function handle
y1 = fx1(xvec);
y2 = fx2(xvec);
y3 = fx3(xvec);
y4 = fx4(xvec);

% 3. plot (xnew, ynew)
plot(xvec,y1,'r-');
plot(xvec,y2,'g-');
plot(xvec,y3,'b-');
plot(xvec,y4,'k-');