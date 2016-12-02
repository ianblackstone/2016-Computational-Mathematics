function f = cannonballvar(t,u)

% Declare constants.
g = 9.81;
m = 15;
c = 0.02;
p = 1.29;
s = 0.25;

% Set the wind to 10 if the time (rounded down) is even.
if mod(floor(t),2)==0
    w=10;
else
    w=0;
end

% Return the derivative of the vector u.
f = [u(3);u(4); -c*p*s/(2*m)*(u(3)-(w))^2; -c*p*s/(2*m)*(u(4))^2-g];

end
