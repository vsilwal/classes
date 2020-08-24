close all
clear all
clc
X = [1:10];
Y = [2 4 5 8 7 8 5 9 12 30];

plot(X, Y,'o');
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
% mrange = -5:10;
% crange = -10:15;
mrange = -5 + 15*rand(N,1);
crange = -10 + 25*rand(N,1);

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
surf(C',M',err_l1)
title('L1 norm')
figure
surf(err_l2)
title('L2 norm')

% find probability
P = exp(-err_l1/1000);
k = sum(P);
p_L1 = P./k;
figure
surf(p_L1);
P = exp(-err_l2/1000);
k = sum(P);
p_L2 = P./k;
figure
surf(p_L2);

% Rejection Method
chance = rand(N*N,1);
ind = find(p_L1(:) > chance);


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
