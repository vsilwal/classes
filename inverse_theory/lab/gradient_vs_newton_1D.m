clear all
close all

N = 10;
xi = linspace(-1,1,N);
x = linspace(-1,1,1000);
dx = x(2) - x(1);

% Test Example
yi = [3 2 2 1 0 -3 -4 1 2 4];

% Interpolate
y = interp1(xi,yi,x,'spline');

% plot
plot(xi, yi, 'o','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
plot(x,y,'Linewidth',2,'Displayname','in-built cubic spline');

%---------------------------------------------------
% Find minimum using gradient descent

% Initial guess
xm = -0.5; 

step_size = .01;


for ii=1:100
    % Numerical estimation of gradient
    [dfdx,ym(ii),ym(ii)] = gradx(xm(ii), x, y);

    % Model update
    xm(ii+1) = xm(ii) - step_size*dfdx;
end

% plot steepest descent points
plot(xm(1:end-1),ym,'o-k','LineWidth',2,'Displayname','gradient descent')

%---------------------------------------------------
% Find minimum using Newton method
% Additional hessian estimation will be required

% Initial guess
xm = 0.4; 

% No step-size required

for ii=1:5
    % Numerical estimation of gradient
    [dfdx,xm(ii),ym2(ii)] = gradx(xm(ii), x, y);
    
    [d2fdx2] = hessx(xm(ii), x, y);
    
    % Model update
    dm = -(dfdx/d2fdx2);
    xm(ii+1) = xm(ii) + dm;
    
    % plot parabola
    % y(dx) = f(x0) + f'(x0)dx + f''(x0)/2 dx^2
    dxvec = linspace(-0.25, 0.25,100);
    yvec = ym2(ii) + dfdx*dxvec + (d2fdx2/2)*(dxvec.*dxvec);
    plot(xm(ii)+ dxvec,yvec,'b')
end

% plot steepest descent points
plot(xm(1:end-1),ym2,'o-g','LineWidth',2,'Displayname','gradient descent')

%===================================================
% Local function for gradient estimation
function [dfdx,xm,ym] = gradx(xm, x, y)

[val,idx]=min(abs(x-xm));
ym = y(idx);
xm = x(idx);

dfdx = (y(idx+1)-y(idx-1))/(x(idx+1)-x(idx-1));
end

%===================================================
% Local function for hessian estimation
function [d2fdx2] = hessx(xm, x, y)

[val,idx]=min(abs(x-xm));

% Numerical second derivative
d2fdx2 = (y(idx+1) - 2*y(idx) + y(idx-1))/((x(idx+1)-x(idx))^2);
end

