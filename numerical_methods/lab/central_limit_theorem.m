close all
clear all

N = 1e6;     % Total Number of samples in Set
N_sub = 3;   % Number of samples in subset (from Set)
N_times1 = 50; % Number of subsets (small)
N_times2 = 1e5; % Number of subsets (large)

% Select distribution
iex = 3;

if iex==1
% Option 1: Homogeneous Distribution
range = 5;   % distribution range [0 5] 
x = range * rand(N,1);

elseif iex==2
% Option 2: Skewed Gaussian Distribution
% XXX

elseif iex==3
% Option 3: Triangular distribution
pd = makedist('Triangular','A',3,'B',200,'C',1000);
x = random(pd,N,1);
end

x_mean = sum(x)/N;
% calculate and check variancee
var1 =  sum((x_mean-x).^2);
var2 = var(x); 

% plot distribution
figure
subplot(2,2,1)
hist(x,20); hold on
ax = axis;
plot([x_mean x_mean], [ax(3) ax(4)], '--r', 'LineWidth',2)
title(sprintf('N = %d, mean = %.2f', N, x_mean))

% grab subset of x 
idx = randi([1 N],N_sub,1);
xsub_mean = sum(x(idx))/N_sub;
subplot(2,2,2)
hist(x(idx),20); hold on
ax2 = axis;
plot([xsub_mean xsub_mean], [ax2(3) ax2(4)], '--r', 'LineWidth',2)
xlim([ax(1) ax(2)])
title(sprintf('N = %d, mean = %.2f', N_sub, xsub_mean))

% Repeat above process multiple times
for ii=1:N_times1
    idx = randi([1 N],N_sub,1);
    xsub_mean(ii)= sum(x(idx))/N_sub;
end
xsub_dist_mean = sum(xsub_mean)/N_times1;
subplot(2,2,3)
hist(xsub_mean,20); hold on
ax3 = axis;
plot([xsub_dist_mean xsub_dist_mean], [ax3(3) ax3(4)], '--r', 'LineWidth',2);
xlim([ax(1) ax(2)])
title(sprintf('%d selections, N = %d, mean = %.2f', N_times1, N_sub, xsub_dist_mean))

% Repeat above process multiple times
for ii=1:N_times2
    idx = randi([1 N],N_sub,1);
    xsub_mean(ii)= sum(x(idx))/N_sub;
end
xsub_dist_mean = sum(xsub_mean)/N_times2;
var_sub = var(xsub_mean); % THis should be close to sigma/sqrt(N_times)
subplot(2,2,4)
hist(xsub_mean,20); hold on
ax3 = axis;
plot([xsub_dist_mean xsub_dist_mean], [ax3(3) ax3(4)], '--r', 'LineWidth',2);
xlim([ax(1) ax(2)])
title(sprintf('%d selections, N = %d, mean = %.2f', N_times2, N_sub, xsub_dist_mean))
