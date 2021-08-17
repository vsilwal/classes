clc
clear all; close all

% set length of X and Y
xlen = 1000; % larger because we want object to move in x-direction
ylen = 100;

% create vectors
xvec = [1:1:xlen];
yvec = [1:1:ylen];

% create 2D mesh of size xlen and ylen
[X,Y] = meshgrid(xvec,yvec);

% Initialize all zeros
Z = zeros(size(X));

% define shape of moving object
xmin=30; xmax=50;
ymin=30; ymax=60;
Z(xmin:xmax,ymin:ymax)=1;

% plot initial condition
figure(1); axis equal
imagesc(Z)

% number of simulation steps
nsteps = 100; 
step_size = 1;

% slide the object
for ii=1:nsteps
XXX % Write this block
    figure(1)
    imagesc(Z)
end

