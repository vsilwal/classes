clc, clear all, close all

iex = 2; % Example number

%% Fourier Example 1 (f(x) = x)
if iex == 1
    x = [-pi:.01:pi]; % range of X
    N = 3; % number of fourier coefficients
    
    % actual function
    fx = @(x) x;
    
    % plot
    figure;
    plot(x,fx(x),'k','LineWidth',3,'DisplayName','f(x)'); grid on
    xlabel('x')
    ylabel('f(x)')
    hold on;
    
    % Compute coefficients (See Lecture Notes)
    a0 = 0;
    
    for ii = 1:N
        an(ii) = 0;
        bn(ii) = (2/ii) * (-1)^(ii+1);
    end
    
    fx_ii = zeros(1,length(x)); % Initialize
    for ii = 1:N
        fourier_coeff = (an(ii)*cos(ii*x) + bn(ii)*sin(ii*x));
        fx_ii = fx_ii + fourier_coeff;
        plot(x,fx_ii,'LineWidth',2,'DisplayName',['N = ',num2str(ii)]) % Uncomment this line to see ALL of them
    end
    legend show
end

%% Newton's method
if iex == 2
    x = [-1.5:.01:2.5]; % range of X
    N = 6; % number of iterations for finding the root
    
    % actual function
    fx = @(x) 2*x.^3 - 3*x.^2 - 3*x + 2;
    dfx = @(x) 6*x.^2 - 6*x -3;
        
    % plot
    figure;
    plot(x,fx(x),'k','LineWidth',3,'DisplayName','f(x)'); grid on
    xlabel('x')
    ylabel('f(x)')
    hold on;
    
    % Iterate over to find roots
    x_ii = -0.16;
    for ii = 1:N
        x_next = x_ii - fx(x_ii)/dfx(x_ii);
        y = XXX %(y - y0) = m(x - x0)

        % plot
        plot(x,y,'LineWidth',2,'DisplayName',['N = ',num2str(ii)])
        % Update x
        x_ii = x_next;
        % print
        disp(sprintf('Iteration = %d, x = %0.3f',ii, x_ii));
    end
    legend show
    plot([x(1) x(end)],[0 0],'--k','DisplayName','x-axis'); 
    plot([0 0],[-10 10],'--k','DisplayName','y-axis');
end
    

