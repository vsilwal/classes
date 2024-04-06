function [c,d,e,f] = add2(a,b)
%
% Perfoms arithmetic operations on two numbers

c = a + b;
d = a - b;
e = a/b;
f = a*b;

%% Example
if 0
    val1 = 2;
    val2 = 4;
    [out1,out2,out3] = add2(val1, val2)
end