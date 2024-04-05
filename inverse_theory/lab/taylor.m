fx = @(x) log(1+x);

x= -1:.01:1;

y = fx(x);

plot(x,y,'.')
hold on
y_taylor = @(x) x - (x.^2/2) + (x.^3/3) - (x.^4/4);
plot(x,y_taylor(x),'.')