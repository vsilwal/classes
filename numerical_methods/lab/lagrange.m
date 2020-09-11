function [val] = lagrange(xi,yi,x)
%
% INPUT
%   x       interpolation points
%   f       function values
%   y       evaluate the interpolation polynomial at points y

%
n = length(xi);
%
% if x is a row vector, transform it to a column vector
%
if size(xi,1)==1
    xi = xi';
end
%
% compute weights using the interpolation points
%
for j=1:n
    dif = xi(j) - [xi(1:j-1);xi(j+1:n)];
    w(j) = 1/prod(dif);
end

m = length(x);
val = zeros(m,1);
%
%  evaluate the interpolation polynomial at y
%
for i = 1: m
    for j=1:n
        XXX
    end
end