% 
% From Pradip Roy
% PhD Candidate
% Department of Earth Science, IITR
%
% 29/10/2021

clear all
close all

% Physical parameters

L = 100;                    % Length of modeled domain [m]
Tmagma = 1200;              % Temperature of magma [C]
Trock = 300;                % Temperature of country rock [C]
kappa = 1e-6;               % Thermal diffusivity of rock [m2/s]
W = 5;                      % Width of dike [m]
day = 3600*24;              % # seconds per day
dt = 3*day;                 % Timestep [s]

% Numerical parameters

nx = 201;                   % Number of gridpoints in x-direction
nt = 100;                   % Number of timesteps to compute
dx = L/(nx-1);              % Spacing of grid
x = -L/2:dx:L/2;            % Grid

% Setup initial temperature profile

T = ones(size(x))*Trock;
T(find(abs(x)<=W/2)) = Tmagma;

% Construct A matrix

A = sparse(nx,nx);
A(1,1)=1;
c = kappa*dt/dx^2;
for i=2:nx-1
    XXX
end
A(nx,nx)=1;
B = XXX;

% construct C matrix
C = sparse(nx,nx);
C(1,1)=1;
for i=2:nx-1
    XXX
end
C(nx,nx)=1;

time = 0;
T = T';
for n=1:nt                  % Timestep loop
    % Compute new temperature
    Tnew = zeros(nx,1);
    Tnew = XXX;
    
    % Set boundary conditions
    Tnew(1) = T(1);
    Tnew(nx) = T(nx);
    
    % Update temperature and time
    T = Tnew;
    time = time+dt;
    
    % Plot solution
    figure(1), clf
    plot(x,Tnew);
    Ymat(n,:) = n*dt*(ones(1,nx));
    Xmat(n,:) = x;
    Tmat(n,:) = Tnew;
    axis([-50 50 200 1300])
    xlabel('x [m]')
    ylabel('Temperature [^oC]')
    title(['Temperature evolution after ',num2str(time/day),' days'])
    drawnow
end
figure
surf(Xmat,Ymat,Tmat)
%shading flat