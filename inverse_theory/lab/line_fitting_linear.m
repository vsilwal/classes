clear all
close all

% Create points
N = 500;
x = linspace(1,10,N);
x = x';

% Select line
m = 2;
c = 3;
error_range = 2;

y = @(m,c) m*x+c;

% Actual line
d_actual = y(m,c);

% Observed data (with error)
errors = error_range*rand(N,1) - (error_range/2); % since rand is [0,1]
errors = error_range*randn(N,1); % since randn return [-1,1] 
d_obs = d_actual + errors;

% forward matrix
% For linear problem this will be equal to the partial derivative matrix
G = [ones(N,1) x];

% KEY: inversion
M = XXX
y_inv = y(M(2), M(1));

% plotting
plot(x,d_obs,'k.', 'MarkerSize',10); hold on
plot(x,d_actual,'-b', 'LineWidth',3)
plot(x,y_inv,'--r', 'LineWidth',2)
legend('data points','actual line','inverted line')

title(sprintf('Line Fitting Example \n m_{actual} = %.2f,c_{actual}= %.2f\n m_{inv}= %.2f, c_{inv}= %.2f',m,c,M(2),M(1)))
