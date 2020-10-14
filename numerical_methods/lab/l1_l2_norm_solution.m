close all
clear all
clc
X = linspace(1,10,10);
Y = [2 4 5 8 7 8 5 9 12 16];

plot(X, Y,'o','MarkerFaceColor','r','MarkerEdgeColor','k');
xlabel('X')
ylabel('Y')
title('Points distribution')

% error functions
err1 = @(y1,y2) sum(abs(y1 - y2));
err2 = @(y1,y2) sum((y1 - y2).^2);

% line
y = @(m,c,x) m*x + c;

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

% grid search
for ii = 1:length(mrange)
    for jj = 1:length(crange)
        ytest = y(mrange(ii),crange(jj),X);
        err_l1(ii,jj) = err1(Y, ytest);
        err_l2(ii,jj) = err2(Y, ytest);
    end
end
M = repmat(mrange,1,N);
C = repmat(crange',N,1);

figure
surf(C',M',err_l1, 'EdgeColor','None')
title('L1 norm')
figure
surf(err_l2,'EdgeColor','None')
title('L2 norm')

%============================================
% find probability
P = exp(-err_l1/1000);
k = sum(P);
p_L1 = P./k;
figure
surf(p_L1,'EdgeColor','None');
P = exp(-err_l2/1000);
k = sum(P);
p_L2 = P./k;
figure
surf(p_L2,'EdgeColor','None');

% Rejection Method
chance = rand(N*N,1);
ind = find(p_L1(:) > chance);
%=============================================

% find minimum
min_l1 = min(err_l1(:));
[row,col] = find(err_l1 == min_l1);
m1 = mrange(row);
c1 = crange(col);
x = linspace(1,10);
figure(1); hold on
plot(x,y(m1,c1,x),'r')
% l2
min_l2 = min(err_l2(:));
[row,col] = find(err_l2 == min_l2);
m2 = mrange(row);
c2 = crange(col);
x = linspace(1,10);
figure(1); hold on
plot(x,y(m2,c2,x),'--k')
