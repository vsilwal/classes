close all; clear all

% Earthquake location 
% ACTUAL - this is for checking whether are results are accurate or not
eqloc_actual = [0,0,10];
v0 = 5; % km/s % Homogeneous earth velocity 
t0 = 0; % origin time of earthquake

%Station locations (x y z)
sta_loc = [35     9     0;
   -44    10     0;
   -11   -25     0;
    23   -39     0;
    42   -27     0;
   -12    50     0;
   -45    16     0;
     5   -19     0;
    -1   -11     0;
    20    11     0];
% X coordiantes of stations
xrec = sta_loc(:,1);
% Y coordiantes of stations
yrec = sta_loc(:,2);
% Z coordiantes of stations
zrec = sta_loc(:,3);

% Number of data and model
N = length(sta_loc);    % 10 stations (= P arrival times)
M = 5;                  % xr,yr,zr,t0,v

% plot station-receiver configuration
figure(1)
plot(xrec,yrec, 'kv','MarkerSize',10,'LineWidth',2); hold on
plot(eqloc_actual(1),eqloc_actual(2),'rp','MarkerSize',10,'LineWidth',2);
xlabel('X-axis (km)'); ylabel('Y-axis (km)');
title('Station and Receiver Configuration');

% Compute travel-time
travel_time = @(xr,yr,zr,xs,ys,zs,t,v) t + (sqrt((xr-xs).^2 + (yr-ys).^2 + (zr-zs).^2 ))/v;

% Compute synthetic observed data
t_obs = travel_time(xrec,yrec,zrec,eqloc_actual(1),eqloc_actual(2),eqloc_actual(3),t0,v0);

% XXX add noise to data
sigma_d = 0.1;
% W = diag(ones(10,1)/sigma_d); % Weight Matrix
% Equivalently you can use Covariance matrix (more appropriate)
Cd = sigma_d*sigma_d * diag(ones(N,1)); 


%% Start Inversion

% Inital model (or starting model) - Prior
eqloc_initial = [3,4,10];   % assumed earthquake location
ti = 2;                     % assumed origin at 2 sec
vi = 4;                  % km/s % assumed velocity

plot(eqloc_initial(1),eqloc_initial(2), 'c*','MarkerSize',10,'LineWidth',2); hold on

% model parameters
m1 = eqloc_initial(1);
m2 = eqloc_initial(2);
m3 = eqloc_initial(3);
m4 = ti;
m5 = vi;

% Create functions for estimating partial derivative matrix (G)
dgdx = @(xr,yr,zr,xs,ys,zs,v) (xs-xr)./(v*sqrt((xr-xs).^2 + (yr-ys).^2 + (zr-zs).^2 ));
dgdy = XXX
dgdz = XXX
dgdt = XXX
dgdv = XXX

n_iterations = 4;

for ii=1:n_iterations
    % Compute partial derivative matrix (G)
    G(:,1) = dgdx(xrec,yrec,zrec,m1,m2,m3,m5);
    G(:,2) = XXX
    G(:,3) = XXX
    G(:,4) = XXX
    G(:,5) = XXX
   
    ts_ii = travel_time(xrec,yrec,zrec,m1,m2,m3,m4,m5);
    delta_d = (t_obs - ts_ii);
    
    % KEY: Generalized Inverse
    dm = XXX

    % Residual
    residual(ii) = sum((t_obs - ts_ii).^2)/length(t_obs);
    
    % Update model
    m1 = m1 + dm(1);
    m2 = m2 + dm(2);
    m3 = m3 + dm(3);
    m4 = m4 + dm(4);
    m5 = m5 + dm(5);
    
    cond(G)
    
    % Plot updated location
    plot(m1,m2, 'm*','MarkerSize',10,'LineWidth',2); hold on
end

% model covariance
Cm = XXX

figure(2)
plot(residual,'o-','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
xlabel('Number of Iteration'); ylabel('Residual')


