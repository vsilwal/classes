% Problem Set1
% Parameter Estimation and Inverse Problems
% Author: Aster, Brocher, and Thurber
% Third Edition (2019)
% 
% Chapter 1
% Section 1.6: Question 3
% 
% Vertical Seismic Profile Example
% 
% Vipul Silwal

clear all
close all

%==== Input Parameters =======
% Number of sensors (Unitless)
%N = 10;
N = 100;
%N = 1000;

% Depth of first and last sensor (in meters)
z_inital = 0;
z_final = 20;

% velocity gradient (v = v0 + kz)
k = 100;    % in m/s/m 
v0 = 1000; % in m/s
%=============================

% discretizing the depth vector
z_points = linspace(z_inital,z_final,N+1);
sensor_depth = z_points(2:length(z_points));

% sensor spacing in depth
dz = z_points(2) - z_points(1);

% creating mid-points for each layer
for i = 1:N
    z_mid(i) = z_points(i) + dz/2;
end

% making design matrix
G_base = tril(ones(N,N));
G = G_base * dz;

% making model_true space
% When our sensor spacing (and hence dz tends to 0) this should be close to
% the d_exact
m_true = 1./(v0 + k*z_mid);
m_true = m_true'; 

% calculating data_true
d_true = G * m_true;

% evaluating data_exact
% exact data can be obtained by analytical solving the integral equation
% Solve the following definate integral
% d(z)  = integral_0_to_z (1/(v0 + ke)) de
d_exact = XXX;
d_exact = d_exact';

% TASKS: 
% 1. Check the difference between d_true and d_exact
% 2. Confirm that the difference decreases as number of sensors (N)
% increases

% plotting
figure;
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_exact,'or');
legend('d\_true','d\_exact');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data space');
grid on;

%======================Inversion=====================
% solving for m
m_exact = XXX;

% TASKS: 
% 1. Check the difference between m_true and m_exact

% plotting
subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m_exact,'or');
legend('m\_true','m\_exact');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Model space');
grid on;

%===================Adding noise=====================
std = .000005; % in seconds
error = (std)*randn(N,1);
d_obs = d_exact + error;
m_inv = XXX;
figure

% TASKS: 
% 1. Check the difference between m_true and m_inv after varying the number
% of sensors
% 2. Check the condition number of matrix G and confirm that it is the
% culprit

% plotting
subplot(2,1,1);
hold on
plot(z_points(2:length(z_points)),d_true,'xb');
plot(z_points(2:length(z_points)),d_obs,'or');
legend('d\_true','d\_exact+noise');
xlabel('z (sensor depth)');
ylabel('data\_space (travel time)');
title('Data space + noise');
grid on;

subplot(2,1,2);
hold on
plot(z_mid,m_true,'xb');
plot(z_mid,m_inv,'or'); 
legend('m\_true','m\_inv');
ylabel('model\_parameter (slowness)');
xlabel('z (interval midpoint)');
title('Inversion after adding noise to data');
grid on;

% ================ additional figures ==================
% Propagating wavefront
if (1)
    figure
    hold on
    
    for i = 1:10:N
    x = v0 * d_exact(i);
    y = sensor_depth(i);
    
    [x1,y1] = ellipse(x,y,pi/1,0,0,'r');
    idx = find (y1>0);
    y1(idx) = 0;
    plot(x1,y1,'r-', 'LineWidth', 3)
    end
    plot([0 0],[0 -20], 'k-', 'LineWidth', 10);  
    axis equal
end



%---------function for plotting Ellipse----------
function [x,y]=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% usage exmple: the following produces a red ellipse centered at 1,1
% and tipped down at a 45 deg axis from the x axis
% ellipse(1,2,pi/4,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%
% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original 
% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch
% Check the number of input arguments 
if nargin<1,
  ra=[];
end;
if nargin<2,
  rb=[];
end;
if nargin<3,
  ang=[];
end;
if nargin<5,
  x0=[];
  y0=[];
end;
 
if nargin<6,
  C=[];
end
if nargin<7,
  Nb=[];
end
% set up the default values
if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;
% work on the variable sizes
x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);
if isstr(C),C=C(:);end;
if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;
% how many inscribed elllipses are plotted
if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;
% drawing loop
for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end;
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k)
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;
  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
  %p=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
  %set(p,'color',C(rem(k-1,size(C,1))+1,:));
  
  %if nargout > 0
  %  h(k)=p;
  end
  
end;
