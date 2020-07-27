% This program describes a moving 1-D wave
% using the finite difference method
clc
close all;
clear all;
%-------------------------------------------------------------------------%
% Initialization
Nx = 101;       % x-Grids
dx = 1;         % Step size
x(:,1) = (0:Nx-1)*dx;
mpx = (Nx+1)/2; % Mid point of x axis
                % ( Mid pt of 1 to 101 = 51 here )
                
T = 1001;       % Total number of time steps
f = 10;         % frequency of source
dt = 0.001;     % Time-Step
t(:,1)= (0:T-1)*dt;
v = 500;        % Wave velocity
c = v*(dt/dx);  % CFL condition
U = zeros(T,Nx);  % U(x,t) = U(space,time)
s1 = floor(T/f);  
%-------------------------------------------------------------------------%
% Initial condition
U((1:s1),1) = sin(2*pi*f.*t(1:s1));
U((1:s1),2) = sin(2*pi*f.*t(1:s1));
%-------------------------------------------------------------------------%
% Finite Difference Scheme
for j = 3:T
    for i = 2:Nx-1
        U1 = 2*U(j-1,i)-U(j-2,i);
        U2 = U(j-1,i-1)-2*U(j-1,i)+U(j-1,i+1);
        U(j,i) = U1 + c*c.*U2;    
    end                   
end
%-------------------------------------------------------------------------%
% Plot the results
plot_times = [50 150 250 310];
for i = 1:4
    
  figure(i)
  k = plot_times(i);
  plot(x,U(k,:),'linewidth',2);
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(t(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end
%-------------------------------------------------------------------------%
% Movie for the travelling wave
for j = 1:T              
  plot(x,U(j,:),'linewidth',2);
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(5);
  set(fh, 'color', 'white'); 
  F=getframe;
            
end
movie(F,T,1)
%-------------------------------------------------------------------------%