function [y,m,G] = quad_spline(xi,yi,x)
% Function for performing quadratic spline interpolation
% 
%-------------------
% Total number of points = N
% Total number of connecting functions = N-1
% Total number of unknowns = 3(n-1) % for quadratic
% Total number of unknowns = 4(n-1) % for cubic (BUT NOT SOLVED HERE)
% 
% Dimnesions: 
% G = 3(N-1) x 3(N-1)
% m = 3(N-1) x 1
% d = 3(N-1) x 1
%-------------------

% check length(xi) == length(yi) 
if length(xi)~=length(yi)
    error('Length of input X and input Y must be equal');
end

n = length(xi);

% pre-allocation
G = zeros(3*(n-1));
y = zeros(length(x),1);

%------------------------------------
% final data vector of length 3(n-1)
% constraint 1
temp = [];
for ii = 2:n-1
    temp = XXX;
end
temp = [yi(1); temp; yi(end)];
% constraint 2
temp = [temp; zeros(n-2,1)];
% contraint 3
temp = [temp; 0];

% Final data vector
d = temp;

%----------------------------------
% Compute G

for ii = 2:n-1
    XXX
    XXX
end
XXX
XXX

% Constraint 2
for ii=2:n-1
    XXX
end

% Constraint 3
G(3*(n-1),1) = 1;
%-----------------------------

% compute coefficients
m = XXX;

% Compute polynomials 
count = 1; 
for ii = 1:length(x)
    if x(ii) <= xi(count+1)
        XXX
    else
        XXX
    end
end

%% =================================================
% EXAMPLE
if 0
    clear all
    close all
    n = 10
    xi = linspace(-1,1,n);
    x = linspace(-1,1,1000);
    % Test Example
    yi = [3 2 2 1 0 -3 -4 1 2 4];
    %yi = [3 2 2 1 0 -3 -4 1 2 4 3 2 2 1 0 -3 -4 1 2 4]; % for n=10
    
    % Function Call
    [y2,m,G] = quad_spline(xi,yi,x);
    
    y3 = interp1(xi,yi,x,'spline');
    
     %plot
    plot(xi, yi, 'o','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k')
    hold on
    
    % plot
    plot(x,y2,'Linewidth',2,'Displayname','quadratic spline');
    % plot
    plot(x,y3,'Linewidth',2,'Displayname','in-built cubic spline');
    legend
end
