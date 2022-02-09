clear all
close all

% define range
xvec = 0:5:180;
yvec = 0:5:180;

[X,Y] = meshgrid(xvec,yvec);

% define function
% NOTE: this is the misfit function (not the synthetic g(m))
fx = @(X,Y) (cosd(X).^2 + cosd(Y).^2).^2;

Z = fx(X,Y);

% plot function
figure(1)
surf(X,Y,Z); 
title(sprintf('f(x,y) = (cos(x)^2 + cos(y)^2)^2'));
xlabel('x'); ylabel('y'); zlabel('f(x,y)')

%---------------------------------------------------------------
% calculate gradient
% ANALYTICALLY
% This can also be used in partial derivative matrix
dfdx = @(X,Y) -4*(cosd(X).*sind(X)).*(cosd(X).^2 + cosd(Y).^2);
dfdy = @(X,Y) -4*(cosd(Y).*sind(Y)).*(cosd(X).^2 + cosd(Y).^2);

% Plot gradient
figure
subplot(2,2,1)
quiver(X,Y,dfdx(X,Y),dfdy(X,Y),'k')
title('Analytical estimation of gradient');
xlabel('x'); ylabel('y');

%---------------------------------------------------------------
% Numerical estimation of gradient (without analytical derivatives)
delx = 5; dely = 5;
% remove boundary points 
xvec_grad = xvec(2):delx:xvec(end-1);
yvec_grad = yvec(2):delx:yvec(end-1);
[numX,numY] = meshgrid(xvec_grad,yvec_grad);

% Numerical gradient using central difference
for ii=2:length(xvec)-1
    for jj=2:length(yvec)-1
        delzdelx(jj-1,ii-1) = XXX;
        delzdely(jj-1,ii-1) = XXX;
    end
end

subplot(2,2,2)
quiver(numX,numY,delzdelx,delzdely,'r')
title('Numerical estimation of gradient');
xlabel('x'); ylabel('y');

%---------------------------------------------------------------
% Compare with MATLAB inbuilt
[px,py] = gradient(Z);
subplot(2,2,3)
quiver(X,Y,px,py)
title('MATLAB in-built gradient');
xlabel('x'); ylabel('y');

%=================================================================
% Find minimum
% Gradient descent
% Hessian estimation not required

% Initial guess
x = 80;
y = 80;
delx = 1; dely = 1;
step_size = 10;

for ii=1:100
    % compute gradient analytically
    dgdx = dfdx(x(ii),y(ii)); dgdy = dfdy(x(ii),y(ii));
    % You can also compute gradient Numerically 
    %dgdx = (fx(x(ii)+delx,y(ii))- fx(x(ii)-delx,y(ii)))/(2*delx);
    %dgdy = (fx(x(ii),y(ii)+dely)- fx(x(ii),y(ii)-dely))/(2*dely);

    % Model update
    x(ii+1) = x(ii) - step_size*dgdx;
    y(ii+1) = y(ii) - step_size*dgdy;
end

figure(1)
hold on
plot3(x,y,fx(x,y),'o-r','LineWidth',2)
quiver(X,Y,px,py)

%=================================================================
% Find minimum
% Newtons method

% Hessian 
d2fdx2 = @(X,Y) -(4*cosd(2*X).*(cosd(X).^2 + cosd(Y).^2) - 2*(sind(2*X)).^2);
d2fdy2 = @(X,Y) -(4*cosd(2*Y).*(cosd(X).^2 + cosd(Y).^2) - 2*(sind(2*Y)).^2);
d2fdxdy = @(X,Y) -2*sind(2*X).*cosd(2*Y);
d2fdydx = @(X,Y) -2*sind(2*X).*cosd(2*Y);

% Initial guess
x = 100;
y = 100;

for ii = 1:100
    Gh = [dfdx(x(ii),y(ii)) dfdy(x(ii),y(ii))];
    
    % second order partial derivative matrix 
    K = [d2fdx2(x(ii),y(ii)) d2fdxdy(x(ii),y(ii)); 
        d2fdydx(x(ii),y(ii)) d2fdy2(x(ii),y(ii))];
        
    % Hessian gives the direction of descent
    H = XXX; 
    
    dm = -inv(H) * Gh' * fx(x(ii),y(ii));
    
    % Model update
    x(ii+1) = x(ii) + dm(1);
    y(ii+1) = y(ii) + dm(2);
end

figure(1)
hold on
plot3(x,y,fx(x,y),'o-k','LineWidth',2)
quiver(X,Y,px,py)
