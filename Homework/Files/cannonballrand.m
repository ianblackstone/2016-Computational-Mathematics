function f = cannonballrand(t,u)

% Declare constants.
g = 9.81;
m = 15;
c = 0.02;
p = 1.29;
s = 0.25;

% Return the derivative of the vector u.
f = [u(3);u(4); -c*p*s/(2*m)*(u(3)-(10*randn))^2; -c*p*s/(2*m)*(u(4))^2-g];

end