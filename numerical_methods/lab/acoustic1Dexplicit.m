clear all
close all

% Parameter Configuration
nx   = 1000;         % number of grid points in x-direction
dx   = 0.5;          % grid point distance in x-direction
c0   = 333.;         % wave speed in medium (m/s)
isrc = 500;          % source location in grid in x-direction
ir   = 730;          % receiver location in grid in x-direction
nt   = 100;         % maximum number of time steps
dt   = 0.0010;       % time step

% CFL Stability Criterion
eps  = c0 * dt / dx; % epsilon value
% This should be less than 1
fprintf('Stability criterion = %f \n', eps)

%% Plot Source Time Function
f0   = 25.; % dominant frequency of the source (Hz)
t0   = 4. / f0; % source time shift
fprintf('Source frequency = %f Hz \n', f0);

% Source time function (Gaussian)
src  = zeros(nt + 1,1);
time = linspace(0 * dt, nt * dt, nt);
% 1st derivative of a Gaussian
src  = -2. * (time - t0) * (f0.^2) .* (exp(-1.0 * (f0.^2)* (time - t0).^2));

%%
% Plot Snapshot & Seismogram (PLEASE RERUN THIS CODE AGAIN AFTER SIMULATION!)

% Initialize empty pressure
p    = zeros(nx,1); % p at time n (now)
pold = zeros(nx,1); % p at time n-1 (past)
pnew = zeros(nx,1); % p at time n+1 (present)
d2px = zeros(nx,1); % 2nd space derivative of p

% Initialize model (assume homogeneous model)
c    = zeros(nx,1);
c    = c + c0;       % initialize wave velocity in model

% Initialize coordinate
x    = 1:nx;
x = x';
x    = x * dx;       % coordinate in x-direction

% Initialize empty seismogram
seis = zeros(nt,1);

%%
for n=1:1:nt
    for j=2:nx-1
        % compute double derivative of p w.r.t x (d2p/dx2)
        d2px(j) = ???;
    end
    
    % Time extrapolation
    pnew = ???;
    
    % Add Source Term at isrc
    % Absolute pressure w.r.t analytical solution
    pnew(isrc) = pnew(isrc) + (dt^2 * src(n))/(dx);
    
    % Remap Time Levels
    pold = p;
    p = pnew;
    
    % Output Seismogram
    seis(n) = p(ir);
    
    plot(x,p)
    ylim([-2e-3 2e-3])
    drawnow
end