%
% Examples 2.1 and 2.2 
% from Parameter Estimation and Inverse Problems, 2nd edition, 2013
% by R. Aster, B. Borchers, C. Thurber
%

clear, close all, clc

% USER CHOICE
PCONF = 0.95;

% Load precomputed data
ddir = './';
load([ddir 'data1.mat']);
t = data1(:,1);
y = data1(:,2);
sigvec = data1(:,3);
ndata = length(t);
nparm = 3;

disp('displaying t, y, sigma')
[ t y sigvec ]

% build the parabolic system matrix
G = [ ones(ndata,1) t -0.5*t.^2 ];

% apply the weighting
yw = y ./ sigvec;
Gw = G ./ [sigvec sigvec sigvec];  % divide each row by an element in sigvec

% solve for the least-squares solution
disp('Least-squares solution')
%m = inv(Gw'*Gw)*Gw'*yw
%m = G\y
m = Gw\yw

% get the covariance matrix
ginv = inv(Gw'*Gw)*Gw';     % generalized inverse, pinv(Gw)
disp('Covariance matrix')
covm = ginv*ginv'

% get the 1.96-sigma (95%) conf intervals
disp('95% parameter confidence intervals (m-, mest, m+)')
%DELTA = 1.96;   % Eq. 2.30
DELTA = erfinv(PCONF)*sqrt(2);
dm = DELTA*sqrt(diag(covm));
[m-dm  m  m+dm]

% N-M degrees of freedom
dof = ndata-nparm;
disp(sprintf('Chi-square misfit for %i degrees of freedom (dof)',dof));
chi2 = norm((y - G*m)./sigvec)^2
%chi2 = norm(yw - Gw*m)^2

% Find the p-value for this data set
% see Eq. 2.21: p value is the integral from chi2 to infinity
disp('chi-square p-value')
p = 1-chi2cdf(chi2,dof)

% Find the parameter correlations
s = sqrt(diag(covm))
disp('correlation matrix')
r = covm./(s*s');

% Plot the data and model predicted data
xx = min(t)-1 : 0.05 : max(t)+1;        % denser sampling of x values
mm = m(1) + m(2)*xx - 0.5*m(3)*xx.^2;

figure(1); hold on
plot(xx,mm,'k');
errorbar(t,y,sigvec,'ko');
xlabel('Time (s)');
ylabel('Elevation (m)');
%print -deps2 c2fparabfig.eps

% Output covm and the eigenvalues/eigenvectors of covm.
disp('Covariance matrix for fitted parameters:')
covm
disp('Eigenvalues/eigenvectors of the inverse covariance matrix:');
[u,lam] = eig(inv(covm))
disp('95% confidence ellipsoid semiaxis lengths:');
semi_axes = [sqrt(chi2inv(PCONF,3)*(1./diag(lam)))]'

disp('95% confidence ellipsoid semiaxes:')
[semi_axes(1)*u(:,1), semi_axes(2)*u(:,2), semi_axes(3)*u(:,3)]

%-------------------------------
% Monte Carlo Section

% predictions for baseline model
y0 = G*m;

nsamp = 1e3;    % number of samples
mmc = zeros(3,nsamp);
chimc = zeros(nsamp,1);
for ii = 1:nsamp
  % Generate a trial data set of perturbed, weighted data
  ytrial = y0 + sigvec.*randn(ndata,1);
  ywtrial = ytrial./sigvec;
  mtrial = Gw\ywtrial;
  % Eq. 2.19
  chimc(ii) = norm( (G*mtrial - ytrial) ./ sigvec )^2;
  %chimc(ii) = norm( ywtrial - Gw*mtrial )^2;
  mmc(:,ii) = mtrial;
end

% Plot the histogram of chi squared values
figure(2); nr=3; nc=1; xmax = 30; ymax = 0.14;

dbin = 0.5;
for ii=1:3
    subplot(nr,nc,ii);
    [B,Bplot] = plot_histo(chimc,[0:dbin:xmax],ii);
    sum(Bplot)*dbin  % check
    xlabel('chi_{obs}^2');
    xlim([0 xmax]); if ii==3, ylim([0 ymax]); end
end

hold on;
xx = linspace(0,xmax,100);
dx = xx(2)-xx(1);
chitheo = chi2pdf(xx,dof);
sum(chitheo)*dx
plot(xx,chitheo,'r','linewidth',2);
title(sprintf('PDF for chi^2(nu=%i, x)',dof))
%ylabel(sprintf('PDF for chi^2(nu=%i, x)',dof))
%xlabel('x')
axis([0 xmax 0 ymax]);

% Plot the histograms of the model parameters
figure(3)
subplot(1,3,1); hist(mmc(1,:)); title('m_1 (m)')
subplot(1,3,2); hist(mmc(2,:)); title('m_2 (m/s)')
subplot(1,3,3); hist(mmc(3,:)); title('m_3 (m/s^2)')

% Plot the realizations of each pair of model parameters with the other
figure; nr=3; nc=2;
ax1 = [-50 50 85 110];
ax2 = [-50 50 7 12];
ax3 = [80 120 7 12];

subplot(nr,nc,1)
plot(mmc(1,:),mmc(2,:),'k.')
xlabel('m_1 (m)'); ylabel('m_2 (m/s)')
axis(ax1);

subplot(nr,nc,3)
plot(mmc(1,:),mmc(3,:),'k.')
xlabel('m_1 (m)'); ylabel('m_3 (m/s^2)')
axis(ax2);

subplot(nr,nc,5)
plot(mmc(2,:),mmc(3,:),'k.')
xlabel('m_2 (m/s)'); ylabel('m_3 (m/s^2)')
disp('Displaying Projections of 1000 Monte-Carlo models (fig 4)')
axis(ax3);

% Plot the 95% error ellipses for each pair of parameters

%generate a vector of angles from 0 to 2*pi
theta = (0:.01:2*pi)';
DELTA = sqrt(chi2inv(PCONF,2));  % 2 is the subspace dimension (p. 34)
%the radii in each direction from the center
r = zeros(length(theta),2);

%---------------

% compute the data for the m1, m2 ellipsoid
C = covm((1:2),(1:2));
[u,lam] = eig(inv(C));
%calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
subplot(nr,nc,2)
plot(m(1)+r(:,1),m(2)+r(:,2),'k');
fill(m(1)+r(:,1),m(2)+r(:,2),'r');
axis(ax1);
xlabel('m_1 (m)'); ylabel('m_2 (m/s)');

% compute the data for the m1, m3 ellipsoid
C = covm([1,3],[1,3]);
[u,lam] = eig(inv(C));
deltachisq = chi2inv(PCONF,2);
DELTA = sqrt(deltachisq);
% calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
% calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m3 ellipsoid
subplot(nr,nc,4)
plot(m(1)+r(:,1),m(3)+r(:,2),'k');
fill(m(1)+r(:,1),m(3)+r(:,2),'r');
axis(ax2);
xlabel('m_1 (m)'); ylabel('m_3 (m/s^2)');

% compute the data for the m2, m3 ellipsoid
C = covm([2,3],[2,3]);
[u,lam] = eig(inv(C));
deltachisq = chi2inv(PCONF,2);
DELTA = sqrt(deltachisq);
% calculate the x component of the ellipsoid for all angles
r(:,1) = (DELTA/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(1,2)*sin(theta);
% calculate the y component of the ellipsoid for all angles
r(:,2) = (DELTA/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(DELTA/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m2, m3 ellipsoid
subplot(nr,nc,6)
plot(m(2)+r(:,1),m(3)+r(:,2),'k');
fill(m(2)+r(:,1),m(3)+r(:,2),'r');
axis(ax3);
xlabel('m_2 (m/s)'); ylabel('m_3 (m/s^2)');
%print -deps2 c2fellipseproj.eps

%==========================================================================