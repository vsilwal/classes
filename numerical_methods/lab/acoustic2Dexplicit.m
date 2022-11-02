clear all
close all
% Parameter Configuration 
% -----------------------

nx   = 500;          % number of grid points in x-direction
nz   = nx;           % number of grid points in z-direction
% Note: regular 2D grid, dz = dx
dx   = 1.;           % grid point distance in x-direction
dz   = dx;           % grid point distance in z-direction
c0   = 580.;         % wave velocity in medium (m/s)
isx  = 250;          % source location in grid in x-direction
isz  = isx;          % source location in grid in z-direction
irx  = 330;          % receiver location in grid in x-direction
irz  = isz;          % receiver location in grid in z-direction
nt   = 750;          % maximum number of time steps
dt   = 0.0010;       % time step

% CFL Stability Criterion
% -----------------------
eps  = c0 * dt / dx; % epsilon value

fprintf('Stability criterion = %f', eps)

%% Plot Source Time Function 
% -------------------------

f0   = 40.; % dominant frequency of the source (Hz)
t0   = 4. / f0; % source time shift

fprintf('Source frequency = %f Hz ', f0)

% Source time function (Gaussian)
% -------------------------------
src  = zeros(nt + 1);
time = linspace(0 * dt, nt * dt, nt);
% 1st derivative of a Gaussian
src  = -2. * (time - t0) * (f0 ^ 2) .* (exp(-1.0 * (f0 ^ 2) * (time - t0).^ 2));

%% Initialize Empty Pressure
% -------------------------
p    = zeros(nz, nx); % p at time n (now)
pold = zeros(nz, nx); % p at time n-1 (past)
pnew = zeros(nz, nx); % p at time n+1 (present)
d2px = zeros(nz, nx); % 2nd space derivative of p in x-direction
d2pz = zeros(nz, nx); % 2nd space derivative of p in z-direction

% Initialize Velocity Model (assume homogeneous model)
% ----------------------------------------------------
c    = zeros(nz, nx);
c    = c + c0;             % initialize wave velocity in model

% Initialize Grid
x    = 1:nx;
x    = x * dx;             % coordinate in x-direction
z    = 1:nz;
z    = z * dz;             % coordinate in z-direction

% Initialize Empty Seismogram
% ---------------------------
seis = zeros(nt);

%%
% Calculate Partial Derivatives
% -----------------------------
for it = 1:1:nt
    
    for i = 2:nx - 1
        d2px(i, :) = XXX;
    end
    
    for j = 2:nz - 1
        d2pz(:, j) = XXX;
    end
    
    % Time Extrapolation
    % ------------------
    pnew = XXX;
    
    % Add Source Term at isz and isx
    % ------------------------------
    % Absolute pressure w.r.t analytical solution
    pnew(isz, isx) = pnew(isz, isx) + (src(it) * (dt ^ 2)) / (dx * dz) ;
    
    % Remap Time Levels
    % -----------------
    pold= p;
    p =  pnew;
    
    % Output Seismogram
    % -----------------
    seis(it) = p(irz, irx);
    
    % Plot
    %surf(p)
    %shading flat
    %colorbar
    imagesc(p)
    drawnow    
end
    hold on
    plot(isx,isz,'rp')
    plot(irx,irz,'rv')
    
figure
plot(time, seis(:,1))
