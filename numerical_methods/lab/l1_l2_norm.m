close all
Nsamples=1000;
X = linspace(1,10,Nsamples);

% line
y = @(m,c,x) m*x + c;
Y = y(1,3,X) + 2*rand(1,Nsamples);

plot(X, Y,'o');
xlabel('X')
ylabel('Y')
title('Points distribution')

% error functions
%err1 = @(y1,y2) XXX;
%err2 = @(y1,y2) XXX;


% approximate range
N = 100;
mmin = -5; mmax = 15;
cmin = -10; cmax = 25;
%mrange = linspace(mmin,mmax,N)';
%crange = linspace(mmin,mmax,N)';
mrange = mmin + mmax*rand(N,1);
crange = cmin + cmax*rand(N,1);
mrange = sort(mrange);
crange = sort(crange);

err = zeros(length(mrange),length(crange));

% grid search
for ii = 1:length(mrange)
    for jj = 1:length(crange)
        % XXX Write something here
    end
end
figure
surf(err_l1)
title('L1 norm')
figure
surf(err_l2)
title('L2 norm')

% find minimum
% l1
% XXX Find minimum of error
m1 = mrange(row);
c1 = crange(col);
x = linspace(1,10);
figure(1); hold on
plot(x,y(m1,c1,x),'r')
% l2
% XXX Find minimum of error
m1 = mrange(row);
c1 = crange(col);
x = linspace(1,10);
figure(1); hold on
plot(x,y(m1,c1,x),'k')