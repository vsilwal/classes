clc, clear all, %close all

% Select example
iex = 5; 

%% Example 1: dy = cos(x)
if iex == 1
    c = [-15:0.5:50]; % a subsetted range of C (in Reality this could be any real value) - constant of integral
    x = [-5:0.1:5]; % pick x range (function actually extends from -inf to + inf)
    
    dy = @(x) cos(x);
    y = @(x) sin(x);  % KEY
    
    figure; hold on; grid on;
    for ii=1:length(c)
        f = y(x) + c(ii);
        plot(x,f)
    end
    xlabel('x');
    ylabel('y');
    title('dy = cos x');
end

%% Example 2: dy = 2y
if iex == 2
    c = [-10:1:10]; % a subsetted range of C (in Reality this could be any real value)
    x = [-0.5:0.1:1.5]; % pick x range (function actually extends from -inf to + inf)
    
    y = @(x) exp(2*x); % KEY
    
    figure; hold on; grid on;
    for ii=1:length(c)
        f = y(x) * c(ii);
        %if c(ii)==1
        plot(x,f,'LineWidth',2)
        %end
    end
    xlabel('x');
    ylabel('y');
    title('dy = 2y');
end

%% Example 3: Direction field (falling object)
if iex == 3
    t = [0:0.5:5];   % sec
    v0 = [0:5:50];   % m/s
    gravity = 9.8;   % m/s^2
    air_drag = 0.47; % gamma
    arrow_len = 0.5;
    
    figure; hold on
    for ii=1:length(t)
        for jj = 1:length(v0)
            dvdt = gravity - (air_drag * v0(jj)); % KEY: tangent (derivative)
            % plot
            u = cos(atan(dvdt)) * arrow_len; % x-component of arrow
            v = sin(atan(dvdt)) * arrow_len; % y-component of arrow
            quiver(t(ii),v0(jj),u,v,'k'); % plot arrow (XXX: something not right)
        end
    end
    xlabel('t');
    ylabel('v');
    title('Direction field for a falling object; $$\frac{dv}{dt} = g - \frac{\gamma v}{m}$$','interpreter','latex')
end

%% Example 4: Direction field [ y' = cos(x) cos(y) ]
if iex == 4
    x = [0:0.4:2*pi];
    y = [-pi:0.4:pi];
    arrow_len = 0.4;
    
    figure; hold on
    for ii=1:length(x)
        for jj = 1:length(y)
            dydx = 2 * cos(x(ii)) * cos(y(jj)); % KEY
            % plot
            u = cos(atan(dydx)) * arrow_len; % x-component of arrow
            v = sin(atan(dydx)) * arrow_len; % y-component of arrow
            quiver(x(ii),y(jj),u,v,'k'); % plot arrow
        end
    end
    xlabel('x');
    ylabel('y');
    title('Direction field for $$y'' = 2 cos(x) cos(y)$$','interpreter','latex')
end

%% Example 5: Euler's method [ dy/dx = 2x ]
%% Go through the class notes for modified Euler method
if iex == 5
    % initial condition
    x0 = 0;
    y0 = 4;
    
    % initialize
    x = [0:0.1:2];
    
    % Differential equaiton and actual solution
    idiff = 1;
    if idiff == 1
        dydx = @(x,y) 2*x;
        y_in = @(x) x.^2; % analytical solution  [ dy/dx = 2x ]
        c = y0 - y_in(x0);
        y = @(x) x.^2 + c;
    elseif idiff == 2
        dydx = @(x,y) -y;
        y_in = @(x) exp(-x);
        c = y0/y_in(x0);
        y = @(x) c*exp(-x) ; % analytical solution  [ dy/dx = -y ]
    end
    
    % Plot actual solution
    %close all; 
    figure
    plot(x,y(x),'k','LineWidth',2); % plot analytical solution
    hold on
    
    % Compute numerical solution
    deltax = 0.4; % discretization in x
    xvec = [x0:deltax:x(end)];
    yvec(1) = y0;
    for ii=1:length(xvec)
        % Apply Euler's method
        yvec(ii+1) = yvec(ii) + (dydx(xvec(ii),yvec(ii)))*deltax; % GENERAL
        %yvec(ii+1) = yvec(ii) + (2*xvec(ii))*deltax; % KEY
        %yvec(ii+1) = yvec(ii) + (-yvec(ii))*deltax; % KEY
        plot(xvec(ii),yvec(ii),'o','MarkerFaceColor','r','MarkerEdgeColor','k')
    end
    
    deltax = 0.1; % smaller discretization (more accurate results)
    xvec = [x0:deltax:x(end)];
    yvec(1) = y0;
    for ii=1:length(xvec)
        yvec(ii+1) = yvec(ii) + (dydx(xvec(ii),yvec(ii)))*deltax; % GENERAL
        %yvec(ii+1) = yvec(ii) + (2*xvec(ii))*deltax; % KEY
        %yvec(ii+1) = yvec(ii) + (-yvec(ii))*deltax; % KEY
        plot(xvec(ii),yvec(ii),'o','MarkerFaceColor','b','MarkerEdgeColor','k')
    end
    
    %--------------------
    % Runge Kutte 4 (RK4)
    
    % Initailize x
    deltax = 0.2; % smaller discretization (more accurate results)
    xvec = [x0:deltax:x(end)];
    yvec(1) = y0;
    
    % Inline functions for k1,k2,k3,k4
    f = dydx;
    k1 = @(x,y) deltax * f(x, y);
    k2 = @(x,y) deltax * f(x + deltax/2, y + k1(x,y)/2);
    k3 = @(x,y) deltax * f(x + deltax/2, y + k2(x,y)/2);
    k4 = @(x,y) deltax * f(x + deltax, y + k3(x,y));
    
    % loop over x
    % compute y(ii+1)
    for ii=1:length(xvec)
        yvec(ii+1) = yvec(ii) + (1/6*(k1(xvec(ii),yvec(ii))...
            + 2*k2(xvec(ii),yvec(ii)) + 2*k3(xvec(ii),yvec(ii)) + k4(xvec(ii),yvec(ii))));
        plot(xvec(ii),yvec(ii),'o','MarkerFaceColor','c','MarkerEdgeColor','k')
    end
    %--------------------
    
    
    xlabel('x'); ylabel('y');
    title('Euler solution to $\frac{dy}{dx} = -y$','interpreter','latex')
end
