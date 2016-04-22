function q = testing( p )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x = p(1);
y = p(2);
q = [x - 0.6*sin(x) - 0.3*cos(y);y - 0.6*cos(x)+ 0.3*sin(y)];

end

