%% EXAMPLE 1:  Trapezoidal rule

% Make an example Parabola
% y = ax^2 + bx + c

% Input coefficients
a = -1;
b = 4;
c = 0;

x = [-1:.1:5];

f = @(a,b,c,x) a*x.^2 + b*x + c;

y = f(a,b,c,x);

% Area using points
h = 1;
x1 = 1; x2 = x1 + h; x3 = x1 + 2*h;
% Y- points 
y1 = f(a,b,c,x1); y2 = f(a,b,c,x2); y3 = f(a,b,c,x3); 

% Area using TRAPEZOIDAL RULE
area1 = 1/2 * h * (y1 + 2*y2 + y3); % KEY

% plot
figure
plot(x,y,'k','LineWidth',2);
grid on, hold on

% plot dashed lines
plot([x1 x1],[0 y1],'--k');
plot([x2 x2],[0 y2],'--k');
plot([x3 x3],[0 y3],'--k');
plot([x1 x2],[y1 y2],'--k');
plot([x2 x3],[y2 y3],'--k');

% color area under the trapezoid
area([x1 x2],[y1 y2]);
area([x2 x3],[y2 y3]);

% plot circluar markers
plot([x1],[y1],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
plot([x2],[y2],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
plot([x3],[y3],'o','MarkerFaceColor','r','MarkerEdgeColor','k');

%% EXAMPLE 2:  Midpoint rule

% Make an example Parabola
% y = ax^2 + bx + c

% Input coefficients
a = -1;
b = 4;
c = 0;

x = [-1:.1:5];

f = @(a,b,c,x) a*x.^2 + b*x + c;

y = f(a,b,c,x);

% Area using points
h = 1;
x1 = 1; x2 = x1 + h; x3 = x1 + 2*h;
% Y- points 
y1 = f(a,b,c,x1); y2 = f(a,b,c,x2); y3 = f(a,b,c,x3); 

% Area using TRAPEZOIDAL RULE
area1 = 1/2 * h * (y1 + 2*y2 + y3); % KEY

% plot
figure
plot(x,y,'k','LineWidth',2);
grid on, hold on

% plot dashed lines
plot([x1 x1],[0 y1],'--k');
plot([x2 x2],[0 y2],'--k');
plot([x3 x3],[0 y3],'--k');
plot([x1 x2],[y1 y2],'--k');
plot([x2 x3],[y2 y3],'--k');

% color area under the trapezoid
area([x1 x2],[y1 y2]);
area([x2 x3],[y2 y3]);

% plot circluar markers
plot([x1],[y1],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
plot([x2],[y2],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
plot([x3],[y3],'o','MarkerFaceColor','r','MarkerEdgeColor','k');
