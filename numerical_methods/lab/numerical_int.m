clc
close all
clear all

iex = 2;

%% EXAMPLE 1:  Trapezoidal rule
if iex ==1
    % Make an example Parabola
    % y = ax^2 + bx + c
    
    % Input coefficients
    a = -1;
    b = 4;
    c = 0;
    
    % Intialize x
    x = [-1:.01:5];
    
    % Initialize function
    f = @(a,b,c,x) a*x.^2 + b*x + c;
    
    % Call the function and find value of Y at X
    y = f(a,b,c,x);
    
    % Discretize X for summation
    xmin = 0; xmax = 4;
    Nsamples = 6; % Number of points for dividing x (N-1 trapezoids)
    x_sum = linspace(xmin,xmax,Nsamples);
    h = x_sum(2) - x_sum(1);
    
    % Y- points at x for summation
    y_sum = f(a,b,c,x_sum);
    
    % Area using TRAPEZOIDAL RULE
    area1 = XXX; % KEY
    
    % plot
    figure
    plot(x,y,'k','LineWidth',2);
    grid on, hold on
    
    % color area under the trapezoid
    for ii = 1:length(x_sum)-1
        area([x_sum(ii) x_sum(ii+1)],[y_sum(ii) y_sum(ii+1)]);
    end
    
    % plot circluar markers
    for ii = 1:length(x_sum)
        plot([x_sum(ii)],[y_sum(ii)],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
    end
    title(sprintf('Trapezoidal Rule \n %dx^2 + %dx + %d; Area = %0.2f; N = %d',a,b,c,area1,Nsamples));
end

%% EXAMPLE 2:  Midpoint rule
if iex ==2
    % Make an example Parabola
    % y = ax^2 + bx + c
    
    % Input coefficients
    a = -1;
    b = 4;
    c = 0;
    
    % Intialize x
    x = [-1:.01:5];
    
    % Initialize function
    f = @(a,b,c,x) a*x.^2 + b*x + c;
    
    % Call the function and find value of Y at X
    y = f(a,b,c,x);
    
    % Discretize X for summation
    xmin = 0; xmax = 4;
    Nsamples = 6; % Number of points for dividing x (N-1 rectangles)
    x_sum = linspace(xmin,xmax,Nsamples);
    h = x_sum(2) - x_sum(1);
    
    % Find mid-points of x
    x_mid = conv(x_sum,[0.5 0.5],'valid');
    
    % Y- points at mid-points for summation
    y_mid = f(a,b,c,x_mid);
    
    % Area using MID-POINT RULE
    area1 = XXX; % KEY
    
    % plot
    figure
    plot(x,y,'k','LineWidth',2);
    grid on, hold on  
    
    % color area under the rectangle
    for ii = 1:length(x_sum)-1
        area([x_sum(ii) x_sum(ii+1)],[y_mid(ii) y_mid(ii)]);
    end
    
    % plot circluar markers
    for ii = 1:length(x_mid)
        plot([x_mid(ii)],[y_mid(ii)],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
    end
    title(sprintf('Mid-point Rule \n %dx^2 + %dx + %d; Area = %0.2f; N = %d',a,b,c,area1,Nsamples))
end

%% EXAMPLE 2:  Simpson's Method
if iex == 3
    % Input coefficients
    a = -1;
    b = 4;
    c = 0;
    
    % Intialize x
    x = [-1:.01:5];
    
    % Initialize function
    f = @(a,b,c,x) a*x.^2 + b*x + c;
    
    % Find y
    y = f(a,b,c,x);
    
    % Discretize X for summation
    xmin = 0; xmax = 4;
    
    % Number of points for dividing x (N-1 rectangles)
    % even N means odd number of subintervals (simpsons 3/8 rule)
    % odd N means even number of subintervals (simpsons even rule)
    Nsamples = 7; 
    
    x_sum = linspace(xmin,xmax,Nsamples);
    h = x_sum(2) - x_sum(1);
    
    % Y- points at x for summation
    y_sum = f(a,b,c,x_sum);
    
    % Call simpson function for integration
    area1 = XXX;
    
    % plot
    figure
    plot(x,y,'k','LineWidth',2);
    grid on, hold on  

    % color area under the curve
    for ii = 1:length(x_sum)-1
        idx = x>x_sum(ii)&x<x_sum(ii+1);
        area(x(idx),y(idx));
    end
    
    % plot circluar markers
    for ii = 1:length(x_sum)
        plot([x_sum(ii)],[y_sum(ii)],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
    end
    if mod(length(x_sum),2)==1
        title(sprintf('Simpsons Rule \n %dx^2 + %dx + %d; Area = %0.2f; N = %d',a,b,c,area1,Nsamples));
    else
        title(sprintf('Simpsons 3/8 Rule \n %dx^2 + %dx + %d; Area = %0.2f; N = %d',a,b,c,area1,Nsamples));
    end
end

% Function for simpson's integration
function area1 = simpson(x,y)
h = x(2)-x(1);
N = length(x);
if mod(N,2)==1 
    % odd number of points (even number of intervals)
    % See class notes for derivation
    area1 = XXX;
else
    % even number of points
    % odd number of intervals (mod(N,2) == 0); % 3/8 rule
    % https://en.wikipedia.org/wiki/Simpson's_rule#Simpson.27s_3.2F8_rule
    area1 = (3*h/8)*(y(1) + sum(3*y(2:3:end-2)) + sum(3*y(3:3:end-1)) + sum(2*y(4:3:end-4)) + y(end));
end
end
