% Problem Set1
% Parameter Estimation and Inverse Problems
% Author: Aster, Brocher, and Thurber
% Third Edition (2019)
% 
% Chapter 1
% Section 1.6: Question 3
% 
% Vertical Seismic Profile Example
% 
% Vipul Silwal

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

% evaluating data_exact
% exact data can be obtained by analytical solving the integral equation
% Solve the following definate integral
% d(z)  = integral_0_to_z (1/(v0 + ke)) de
d_exact = XXX;
d_exact = d_exact';

% TASKS: 
% 1. Check the difference between d_true and d_exact
% 2. Confirm that the difference decreases as number of sensors (N)
% increases

% plotting
figure;
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_exact,'or');
legend('d\_true','d\_exact');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data space');
grid on;

%======================Inversion=====================
% solving for m
m_exact = G\d_exact;

% TASKS: 
% 1. Check the difference between m_true and m_exact

% plotting
subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m_exact,'or');
legend('m\_true','m\_exact');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Model space');
grid on;

%===================Adding noise=====================
std = .000005; % in seconds
error = (std)*randn(N,1);
d_obs = d_exact + error;
m_inv = G\d_obs;
figure

% TASKS: 
% 1. Check the difference between m_true and m_inv after varying the number
% of sensors
% 2. Check the condition number of matrix G and confirm that it is the
% culprit

% plotting
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_obs,'or');
legend('d\_true','d\_exact+noise');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data space + noise');
grid on;

subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m_inv,'or'); 
legend('m\_true','m\_inv');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Inversion after adding noise to data');
grid on;

% ================ additional figures ==================
% Propagating wavefront
if (1)
    figure
    hold on
    
    for i = 1:10:N
    x = v0 * d_exact(i);
    y = sensor_depth(i);
    
    [x1,y1] = ellipse(x,y,pi/1,0,0,'r');
    idx = find (y1>0);
    y1(idx) = 0;
    plot(x1,y1,'r-', 'LineWidth', 3)
    end
    plot([0 0],[0 -20], 'k-', 'LineWidth', 10);  
    axis equal
end
