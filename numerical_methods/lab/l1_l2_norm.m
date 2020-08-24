close all
X = [1:10];
Y = [2 4 5 8 7 8 5 9 12 16];

plot(X, Y,'o');
xlabel('X')
ylabel('Y')
title('Points distribution')

% error functions
%err1 = @(y1,y2) XXX;
%err2 = @(y1,y2) XXX;

% line
y = @(m,c,x) m*x + c;

% approximate range
mrange = -5:10;
%crange = Range of intercept;

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
% XXX Find minimum of error
m1 = mrange(row);
c1 = crange(col);
x = linspace(1,10);
figure(1); hold on
plot(x,y(m1,c1,x),'r')
% l2
