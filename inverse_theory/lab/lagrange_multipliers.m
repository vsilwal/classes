clear all
close all

% get function values
% function to be minimized
N = 100;
[X,Y,Z] = peaks(N);

% get gradient
[dfdx,dfdy] = gradient(Z,X(1,2)-X(1,1),Y(2,1)-Y(1,1));

% plot function
s=surf(X,Y,Z); 
%s.EdgeColor = 'none';
hold on
xlabel('x'); ylabel('y');
% plot gradient
quiver(X,Y,dfdx,dfdy,'k')
title('3(1-x)^2 e^{(-(x^2) - (y+1)^2)} - 10(x/5 - x^3 - y^5) e^{-x^2-y^2} - 1/3 e^{-(x+1)^2 - y^2}')

% Now see in contour
figure
contour(X,Y,Z); hold on
quiver(X,Y,dfdx,dfdy,'k')
xlabel('x'); ylabel('y');

%==================================================
if 1
% gradient descent
% initial guess

xm = 1.5;
ym = -0.5;

% index corresponding to this location
[idx,idy] = return_index(X,Y,xm,ym);
% check
plot(X(idy,idx), Y(idy,idx),'og','MarkerSize',5)

% perform iteration
step_size = .1;

for ii = 1:40    
    [idx(ii),idy(ii)] = return_index(X,Y,xm(ii),ym(ii));
    
    dfdx_xm = dfdx(idy(ii),idx(ii));
    dfdy_ym = dfdy(idy(ii),idx(ii));
    
    % Model update
    xm(ii+1) = xm(ii) - step_size*dfdx_xm;
    ym(ii+1) = ym(ii) - step_size*dfdy_ym;
end

plot(xm,ym,'o-r','LineWidth',2)
title('gradient descent');
figure(1)
plot3(xm(1:end-1),ym(1:end-1),Z(sub2ind([N N],idy,idx)),'o-r','LineWidth',2)
end

%==================================================
% constraint function
% Let circle x^2 + y^2 = 1 be the constraint function
theta = [0:.05:2*pi];
radius = 1;
x = radius*cos(theta);
y = radius*sin(theta);

% Get location of constarin function in mesh
for ii = 1:length(x)
    [idx(ii),idy(ii)] = return_index(X,Y,x(ii),y(ii));
    
    dfdx_xm(ii) = dfdx(idy(ii),idx(ii));
    dfdy_ym(ii) = dfdy(idy(ii),idx(ii));

end
figure
% plot minimizing function
s=surf(X,Y,Z); 
%s.EdgeColor = 'none';
hold on
xlabel('x'); ylabel('y');
% plot constrain function
plot3(x,y,Z(sub2ind([N N],idy,idx)),'-r','LineWidth',5)
title('Peaks')

% plot gradient along circle
figure
contour(X,Y,Z); hold on
plot(x,y,'k');
quiver(x,y,dfdx_xm,dfdy_ym,'b');
quiver(x,y,2*x,2*y,'r')
%quiver(X,Y,dfdx,dfdy,'k')

% Check if minimizing function and contraint function gradient vectors
% align
% Use cross product
v1 = [dfdx_xm' dfdy_ym'];
v1_norm = (v1(:,1).^2 + v1(:,2).^2).^(1/2);
v2 = [2*x' 2*y'];
v2_norm = (v2(:,1).^2 + v2(:,2).^2).^(1/2);


%==================================================
% function for finding index specific to a location in mesh
function [idx,idy] = return_index(X,Y,xm,ym)
[valx,idx] = min(min(abs(X-xm)));
[valy,idy] = min(min(abs(Y'-ym)));
end