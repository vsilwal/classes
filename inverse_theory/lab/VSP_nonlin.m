% Vertical Seismic Profile Example
% Using non-linear approach
%
% Vipul Silwal
%
% For liner solution:
% Parameter Estimation and Inverse Problems
% Author: Aster, Brocher, and Thurber
% Third Edition (2019)
%
% Chapter 1
% Section 1.6: Question 3

clear all
close all

%==== Input Parameters =======
% Number of sensors (Unitless)
%N = 10;
N = 100;
%N = 1000;

% Depth of first and last sensor (in meters)
z_inital = 0;
z_final = 20;

% velocity gradient (v = v0 + kz)
k = 100;    % in m/s/m
v0 = 1000; % in m/s
%=============================

% discretizing the depth vector
z_points = linspace(z_inital,z_final,N+1);
sensor_depth = z_points(2:length(z_points));
sensor_depth = sensor_depth';

% sensor spacing in depth
dz = z_points(2) - z_points(1);

% creating mid-points for each layer
for i = 1:N
    z_mid(i) = z_points(i) + dz/2;
end

% making design matrix
G_base = tril(ones(N,N));
G = G_base * dz;

% making model_true space
% When our sensor spacing (and hence dz tends to 0) this should be close to
% the d_exact
m_true = 1./(v0 + k*z_mid);
m_true = m_true';

% calculating data_true
d_true = G * m_true;

%-----------------------
% g(m) non-linear forward model
% k - scaler; v0 - scaler; z - vector
gm = @(k,v0,z) (1/k)*(log(v0+k*z) - log(v0));

% Exact solution
d_exact = gm(k,v0,sensor_depth);

% Let observed be equal to exact model
d_obs = d_exact;
% Or add noise to it
std = .00005; % in seconds
error = (std)*randn(N,1);
d_obs = d_exact + error;

%-----------------------
% partial derivatives
% for model update using approximate Hessian
dgdk = XXX;
dgdv0 = XXX;

% Funcitons for second order partial derivative matrix
% for model update using true Hessian
% H(m) = hat(G)^T *hat(G) + sum_over_Ndata(g(m)-d)Ki
d2gdk2 = XXX;
d2gdv02 = XXX;
d2gdkdv0 = XXX;
d2gdv0dk = XXX; % same as d2Hdkdv0

%-----------------------
% Initial Guess
k_initial = 70;
v0_initial = 700;

% Initial synthetic data points
d_initial = gm(k_initial,v0_initial,sensor_depth);

% Check (and compare)
figure
subplot(1,2,1)
plot(sensor_depth,d_initial,'xb'); hold on
plot(sensor_depth,d_true,'or');
xlabel('sensor deopth')
ylabel('travel-time (sec)')
%-----------------------
% iteration block for 'approximate Hessian' non-linear case
n_iterataion = 5;

m1 = k_initial;
m2 = v0_initial;

for ii=1:n_iterataion
    % prepare partial derivative matrix (N X M)
    Gh(:,1) = dgdk(m1(ii),m2(ii),sensor_depth);
    Gh(:,2) = dgdv0(m1(ii),m2(ii),sensor_depth);
    
    % synthetic model
    d_ii = gm(m1(ii),m2(ii),sensor_depth);
    % data residual
    delta_d(:,ii) = d_obs - d_ii;
    d_residual(ii) = norm(delta_d(:,ii));
    
    % Generalized inverse
    % using approximate hessian (ignore second order partial derivate
    % matrix K
    dm(:,ii) =  XXX
    
    % update model
    m1(ii+1) = m1(ii) + dm(1,ii);
    m2(ii+1) = m2(ii) + dm(2,ii);
    plot(sensor_depth,d_ii)
    text(sensor_depth(end)-5,d_ii(end),sprintf('iteration = %d',ii))
end

% plot Residual with iteration
subplot(1,2,2)
plot(d_residual,'-o')
xlabel('iteration')
ylabel('misfit norm')

%---------------------------------
% iteration block for 'true Hessian' non-linear case
m1 = k_initial;
m2 = v0_initial;

figure; 
subplot(1,2,1)
hold on
plot(sensor_depth,d_initial,'xb'); hold on
plot(sensor_depth,d_true,'or');
xlabel('sensor deopth')
ylabel('travel-time (sec)')
for ii=1:n_iterataion
    % synthetic model
    d_ii = gm(m1(ii),m2(ii),sensor_depth);
    % data residual
    delta_d(:,ii) = d_obs - d_ii;
    d_residual(ii) = norm(delta_d(:,ii));
    
    % prepare partial derivative matrix (N X M)
    Gh(:,1) = dgdk(m1(ii),m2(ii),sensor_depth);
    Gh(:,2) = dgdv0(m1(ii),m2(ii),sensor_depth);
    
    % Using true hessian
    % First calculate Ki (second order partial derivative matrix)
    % Think of it as a 3D block of (M X M X N) dimension
    Ki(1,1,:) = XXX;
    Ki(1,2,:) = XXX;
    Ki(2,1,:) = XXX;
    Ki(2,2,:) = XXX;
    
    % Compute second order changes (ignored in approximate hessian)
    % Think each layer of K (M X M X N) is multiplied by data residual
    for jj=1:N
        K_temp(:,:,jj) = XXX;
    end
    % Now sum all data layer 
    % sum_over_Ndata(g(m)-d)Ki
    sum_K = sum(K_temp,3);
    
    % Compute Hessian
    % H(m) = hat(G)^T *hat(G) + sum_over_Ndata(g(m)-d)Ki
    H = XXX;
    
    % Generalized inverse using true hessian 
    dm2(:,ii) =  XXX;
    
    % update model
    m1(ii+1) = m1(ii) + dm2(1,ii);
    m2(ii+1) = m2(ii) + dm2(2,ii);
    plot(sensor_depth,d_ii)
    text(sensor_depth(end)-5,d_ii(end),sprintf('iteration = %d',ii))
end

% plot Residual with iteration
subplot(1,2,2)
plot(d_residual,'-o')
xlabel('iteration')
ylabel('misfit norm')
%---------------------------------

% plot surface plot - (k,v0,misfit)
k_vec = [50:10:150];
v0_vec = [00:100:7000];
[K,V0] = meshgrid(k_vec,v0_vec);

for ii = 1:N
    for jj = 1:length(k_vec)
        for kk = 1:length(v0_vec)
           d_syn(kk,jj,ii) = gm(k_vec(jj),v0_vec(kk),sensor_depth(ii));
           d_res(kk,jj,ii) = norm(d_syn(kk,jj,ii) - d_obs(ii));
        end
    end
end
%Total residual for all combinations of k and v0
d_res_total = sum(d_res,3);
figure
imagesc(d_res_total)
