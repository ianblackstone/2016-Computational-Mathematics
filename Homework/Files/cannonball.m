function f = cannonball(t,u)

g = 9.81;
m = 15;
c = 0.02;
p = 1.29;
s = 0.25;
w = 0;

f = [ u(3); u(4); -c*p*s/(2*m)*(u(3)-w)^2; -c*p*s/(2*m)*u(4)-g*u(2)/u(1)];

end