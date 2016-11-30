function rhs = dblPend(t,u)

g = 9.81;
c = cos(u(1)-u(2));
s = sin(u(1)-u(2));

M = [[1 0 0 0];[0 1 0 0];[0 0 6 2*c];[0 0 4*c 2]];
f = [u(3); u(4); -3*g*sin(u(1))-2*s*u(4)^2; -2*g*sin(u(2))+4*sin(u(3))^2];

rhs = M\f;

end