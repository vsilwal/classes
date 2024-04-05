% clear everything
clc 
clear all
close all

% Set length of X and Y
xlen = 1000;
ylen = 100;

% create vector
xinc = 1; yinc = 1;
xvec = [1:xinc:xlen];
yvec = [1:yinc:ylen];

% create a 2D mesh
[X,Y] = meshgrid(xvec,yvec);

Z = zeros(size(X)); % zeros(100,1000)

% define the shape of object
xmin = 30; xmax = 50;
ymin = 40; ymax = 80;
Z(xmin:xmax,ymin:ymax) = 1;

figure(1)
imagesc(Z)

% number of simulation steps
nsteps = 100;
step_size = 10;

% loop over
for ii=1:nsteps
    Z(xmin:xmax,ymin+ii-step_size:ymax+ii-step_size) = 0;
    Z(xmin:xmax,ymin+ii:ymax+ii)=1;
    figure(2)
    imagesc(Z)
end