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
%t_obs = t_obs + rand(10,1);
% W = diag(ones(10,1)/sigma_d); % Weight Matrix
% Equivalently you can use Covariance matrix (more appropriate)
Cd = sigma_d*sigma_d * diag(ones(N,1)); 

%% Perform grid search
% Homogeneously sample X and Y (keep Z,t,v same as true model -  this
% is done only for plotting understanding a simple problem; in reality we
% will perform grid search over all components

Nsamples = 1e5;
xsamp = -50 + 100*rand(Nsamples,1);
ysamp = -40 + 90*rand(Nsamples,1);
zsamp = eqloc_actual(3)*ones(Nsamples,1);
tsamp = t0*ones(Nsamples,1);
vsamp = v0*ones(Nsamples,1);

for ii=1:Nsamples
    t_samp(:,ii) = travel_time(xrec,yrec,zrec,xsamp(ii),ysamp(ii),zsamp(ii),tsamp(ii),vsamp(ii));
    residual_samp(:,ii) = t_obs - t_samp(:,ii);
    residual_samp_all(ii,1) = sum(residual_samp(:,ii).^2)/length(t_obs);
end

figure
scatter(xsamp,ysamp,100,residual_samp_all,'filled'); hold on; colorbar; colormap('jet')
plot(xrec,yrec, 'kv','MarkerSize',10,'LineWidth',2); hold on
plot(eqloc_actual(1),eqloc_actual(2),'rp','MarkerSize',10,'LineWidth',2);
xlabel('X-axis (km)'); ylabel('Y-axis (km)');
title(sprintf('Residuals for %.0f random samples;\n Only forward simulation - No inversion',Nsamples));

% Residuals for one station only
figure
sta_idx = 1;
scatter(xsamp,ysamp,100,abs(residual_samp(sta_idx,:)),'filled'); hold on; colorbar; colormap('jet')
plot(xrec,yrec, 'kv','MarkerSize',10,'LineWidth',2); hold on
plot(eqloc_actual(1),eqloc_actual(2),'rp','MarkerSize',10,'LineWidth',2);
xlabel('X-axis (km)'); ylabel('Y-axis (km)');
title(sprintf('Residuals at Station %d for %.0f random samples;\n Only forward simulation - No inversion',sta_idx,Nsamples));


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
    cond(G)
    
    % Residual
    residual(ii) = sum((t_obs - ts_ii).^2)/length(t_obs);
    
    % Update model
    m1 = m1 + dm(1);
    m2 = m2 + dm(2);
    m3 = m3 + dm(3);
    m4 = m4 + dm(4);
    m5 = m5 + dm(5);
    
    % Plot updated location
    plot(m1,m2, 'm*','MarkerSize',10,'LineWidth',2); hold on
end

% model covariance
Cm = XXX

% 95% confidence interval (Aster (Eq 2.32))
m_err = 1.96*sqrt(diag(Cm));

sprintf('Final Solution:\n x = %.2f +- %.2f,\n y = %.2f +- %.2f,\n z = %.2f +- %.2f',m1,m_err(1),m2,m_err(2),m3,m_err(3))

figure
plot(residual,'o-','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10)
xlabel('Number of Iteration'); ylabel('Residual')

%% Plot location with error ellipsoid
% This is when we ALSO include off-diagonal elements of model covariance matrix (Cm)
% 
% Error estimates in this case are larger than those estimated in last
% section using the diagonal elements only

% get 95% conf intervals (See Aster (Eq 2.33))
PCONF = 0.95; % 95% percentile of the chi-square distribution with 2 degrees of freedom
dof = 2; % (m1,m2) (m2,m3) etc... projecting on 2 model parameter space at a time 
deltachisq = chi2inv(PCONF,2);
DELTA = sqrt(deltachisq);

%-------------(m1,m2)--------------------
% Get (x,y) subset of Cm (NOT the z,t and v components)
C = Cm((1:2),(1:2));

% get eigenalues and eigenvectors of model covariance vector
% The eigen vectors will be the semi-major axis of the error ellipses
[u,lam] = eig(inv(C));

%generate a vector of angles from 0 to 2*pi
theta = (0:.01:2*pi)';
%calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
figure; 
subplot(1,3,1)
hold on
plot(m1+r(:,1),m2+r(:,2),'k');
fill(m1+r(:,1),m2+r(:,2),'r');
plot(m1,m2, 'mp','MarkerSize',10,'LineWidth',2);
xlabel('m_1 (x)'); ylabel('m_2 (y)');

% Similarly (m2,m3) and (m1,m3)

%------------(m2,m3)---------------------
% Get (x,y) subset of Cm (NOT the z,t and v components)
C = Cm((2:3),(2:3));

% get eigenalues and eigenvectors of model covariance vector
% The eigen vectors will be the semi-major axis of the error ellipses
[u,lam] = eig(inv(C));

%generate a vector of angles from 0 to 2*pi
theta = (0:.01:2*pi)';
%calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
subplot(1,3,2)
hold on
plot(m2+r(:,1),m3+r(:,2),'k');
fill(m2+r(:,1),m3+r(:,2),'r');
plot(m2,m3, 'mp','MarkerSize',10,'LineWidth',2);
xlabel('m_2 (y)'); ylabel('m_3 (z)');

%------------(m1,m3)---------------------
% Get (x,y) subset of Cm (NOT the z,t and v components)
C = Cm([1,3],[1,3]);

% get eigenalues and eigenvectors of model covariance vector
% The eigen vectors will be the semi-major axis of the error ellipses
[u,lam] = eig(inv(C));

%generate a vector of angles from 0 to 2*pi
theta = (0:.01:2*pi)';
%calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
subplot(1,3,3)
hold on
plot(m1+r(:,1),m3+r(:,2),'k');
fill(m1+r(:,1),m3+r(:,2),'r');
plot(m1,m3, 'mp','MarkerSize',10,'LineWidth',2);
xlabel('m_1 (x)'); ylabel('m_3 (z)');

