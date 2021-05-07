% Cross-well tomography
% Problem 4.10.3
% Parameter Estimation and Inverse Problems
% Author: Aster, Brocher, and Thurber
% Third Edition (2019)

% Vipul Silwal
% 29/4/2021

close all
clear all

% initialize config
dz = 100;       % vertical spacing between sources (and receivers)
z_max = -50;    % position of top most source
z_min = -1550;  % position of bottom most source

% SOURCE
zs = [z_min:dz:z_max]';
xs = zeros(length(zs),1);
Ns = length(zs); % Number of sources
% RECEIVER
zr = zs;
xr = zeros(length(zr),1)+1600; % receiver borehole is 1600 mt away

% plot configuration
x_vec = 0:dz:1600;
z_vec = -1600:dz:0;

hold on
for ii = 1:length(x_vec) %// loop over vertical lines
    plot([x_vec(ii) x_vec(ii)], [z_vec(1) z_vec(end)], 'k-');
end
for ii = 1:length(z_vec) %// loop over horizontal lines
    plot([x_vec(1) x_vec(end)], [z_vec(ii) z_vec(ii)], 'k-');
end

% plot source borehole
plot([x_vec(1) x_vec(1)], [z_vec(1) z_vec(end)], 'k-','Linewidth',5)
% plot receiver borehole
plot([x_vec(end) x_vec(end)], [z_vec(1) z_vec(end)], 'k-','Linewidth',5)
yyaxis left; ylabel('Source Borehole');
yyaxis right; ylabel('Receiver Borehole');

% plot sources
plot(xs, zs, 'k+','MarkerSize',10)
% plot receiver
plot(xr, zr, 'k+','MarkerSize',10)

% plot an example ray path
sid = 7; % source id

for ii = 1:length(zs) %// loop over horizontal lines
    plot([x_vec(1) x_vec(end)], [zs(sid) zs(ii)], 'r-');
end

%% Inversion begins

% Solve using SVD

load crosswell.mat % load 'G' matrix and noisy data 'dn'

mu = 0.1;

rank(G) % Check the rank of G, Does it have non-trivial null space?

% singular value decomposition
% Check out matlab functions: svd, pinv
[U,S,V] = svd(G);
Gcheck = U*S*V'; % because G=USV'

% Gdagger is pseudoinverse of G
% Gdagger = V * inv(S) * U'

Gdagger = XXX; % Moore-Penrose Pseudoinverse

%Gdagger = XXX; % use pinv and set the threshold

m_svd = Gdagger*dn;
M_svd = reshape(m_svd,16,16);
figure; pcolor(M_svd)

%% Solve using damped least-square

% XXX
%% Comapre the results