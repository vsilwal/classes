clc, clear all, close all

iex = 1; % Example number

%% EXAMPLE 1: Homogeneous probability density
if iex == 1
    % 2D box
    N = 1000; % Number of sample points
    x = rand(N,1);
    y = rand(N,1);
    figure
    plot(x,y,'.')
    
    %------------------------------------------------
    % sphere
    % x = r cos(theta) sin (phi)
    % y = r sin(theta) sin (phi)
    % z = r cos(phi)
    % theta = [0 2pi) from x-axis on x-y plane
    % phi = [0 pi] from +z axis to -z axis
    %
    % homogeneous distribution
    N = 5000;
    r = 5; % radius
    theta = 2 * pi * rand(N,1);
    %phi = pi * rand(N,1);           % INCORRECT  (Even though phi is homogeneous the samples generated on the surface of the sphere are NOT)
    phi = acos(2 * rand(N,1) - 1);   % CORRECT
    
    x = r .* cos(theta) .* sin(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(phi);
    
    figure
    subplot(2,2,1)
    scatter3(x,y,z,'.k');  % Check out x-y plane view
    subplot(2,2,2)
    scatter3(x,y,z,'.k');  % Check out x-y plane view
    subplot(2,2,3)
    hist(theta)
    subplot(2,2,4)
    hist(phi)
    xlim([0 4])
    axis equal
end
%------------------------------------------------

%% EXAMPLE 2: Normal Distribution
if iex == 2
    x = [-10:.1:10];
    
    % create inline function
    %XXX f = @(x,mu,sigma)(1/sqrt(2*pi*sigma^2)) * exp(-((x-mu).^2)/2*sigma^2);
    f = @(x,mu,sigma)(1/sqrt(2*pi*sigma^2)) * exp(-((x-mu).^2)/(2*sigma^2));
    
    % standard normal distribution
    mu = 0;
    sigma = 1;
    f1 = f(x,mu,sigma);
    
    mu = 0;
    sigma = 2;
    f2 = f(x,mu,sigma);
    
    mu = 5;
    sigma = 2;
    f3 = f(x,mu,sigma);
    
    % plotting
    figure
    hold on
    plot(x,f1)
    plot(x,f2)
    plot(x,f3)
    legend('mu = 0, sigma =1','mu = 0, sigma =2','mu = 5, sigma =2')
    grid on
    xlabel('x')
    ylabel('f(x)')
    
    % Problem
    % 1. Find the area under each of these curves
end