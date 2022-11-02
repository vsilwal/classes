% 2D Advection Equation
%
% du/dt + v_x du/dx = 0
% Future: Could also modify to represent Advection-convection equation

clc,clear,close all

vx = 5; % m/s
vy = 3; % m/s

dt = .4;

xmin = -500; xmax = 500;
ymin = -500; ymax = 500;

Nx = 41;
Ny = 41;
Nt = 200;

%------------------------
% mesh domain
xvec = linspace(xmin,xmax,Nx);
yvec = linspace(ymin,ymax,Ny);
dx = (xmax-xmin)/(Nx-1);
dy = (ymax-ymin)/(Ny-1);

[X,Y]=meshgrid(xvec,yvec);

% Inital condition
% Gaussian source
r = 50; % variance
u = exp(-((X - xvec(nearest(Nx/2))).^2 + (Y - yvec(nearest(Ny/2))).^2)/(2*r^2)); 

% Use square box to see numerical artifacts
%u = zeros(Nx,Ny);
%u(nearest(Nx/2)-5:nearest(Nx/2)+5,nearest(Ny/2)-10:nearest(Ny/2)+10) = 5;

figure(1)
subplot(2,2,1)
surf(u)

% Initialize for next step
unew = zeros(Nx,Ny);

rx = 0.5*vx*dt/dx;
ry = 0.5*vy*dt/dy;

% loop over time
for n = 1:Nt
    
    for jj = 2:Ny-1
        for ii = 2:Nx-1
            unew(ii,jj) = ???
        end
    end
    
    % update for next step
    u = unew;
    
    % Boundary condition (Dirichlet)
    u(1,:)  = 0; u(Nx,:) = 0;
    u(:,Ny) = 0; u(:,1) = 0;
    
    % plot
    figure(1), 
    subplot(2,2,2)
    surf(u)
    
    subplot(2,2,3)
    imagesc(u)
    
end
colorbar
