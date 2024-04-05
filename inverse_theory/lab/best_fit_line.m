clear all
close all

% initial points
x = [1 2 3 4 5 6 7 8];
yobs = [ 1 5 3 7 8 11 10 11];

% plot initial dataset (in red)
plot(x,yobs,'o','MarkerSize',10,'MarkerFaceColor','r')

% draw a test line (y = mx + c)
xvec = 1:0.01:8;
m = (yobs(end)-yobs(1))/(x(end)-x(1));
c = yobs(1);

for m = 0:0.5:2
    for c = -10:2:10
        yvec = m*xvec + c;

        % plot line
        hold on
        plot(xvec, yvec)

        % plot projection of initial dataset on line (in blue)
        ysyn = m*x + c;
        plot(x,ysyn,'o','MarkerSize',10,'MarkerFaceColor','b')

        % calculate error
        err_l1  = sum(abs(yobs-ysyn));
        %err_l2  = yobs  - ysyn
    end
end

plot(x,yobs,'o','MarkerSize',10,'MarkerFaceColor','r')