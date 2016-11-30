%% Problem 5 notes
%F = (Xf(t),Yf(t))
%R = (cos(t),sin(t))
%u = [Xf ; Yf] = [u1 ; u2]
%F' = 1/sqrt((u1-cos(t))^2 + (u2-sin(2))^2) * [ u1-cos(t) ; u2-sin(t) ]
%make F' a function to find the current velocity of the fox.

%% Paratrooper 7.17
fuction f = paratrooper(t,y)
K1 = 1/15;
K2 = 4/15;
m = 80;

if t< 5
    alpha = K1*y(2)^2;
else
    alpha = K2*y(2)^2;
end
f = [y(2); -9.81 + alpha/m];
end

y0 = [600; 0];
[t,y] = ode45(@paratrooper,[0 5],??]

% event detection
function [gstop,isstop,direction] = ground(t,y)

% gstop is 0 when we want to stop the simulaion.  y(1) stores the current
% position so when the paratrooper is on the ground the simuation will
% stop.
gstop = y(1);

% isstop is for hether we want to continue the smulation or not, boolean.
isstop = 1;

% Direction limits the directions we are allowed to approach 0 from.  Sqrt
% problems must come from the right.  For this the answer is the same
% regardless so it is direction invariant.
direction = [];

end
% save function as ground in its own file.

% detect event
opt = odeset('event',@ground);
[t,y] = ode45(@paratrooper,[0,inf],y0,opt);

% the last time entry is the time the event ocurred at.
t(end)

% y stores [position velocity] so the last velocity will be the impact
% velocity.
y(end,2)

%% Double Pendulum

x1 = l1*sin(theta1);
y1 = -l1*cos(theta1);
x2 = x1 + l2*sin(theta2);
y2 = y1 - l2*cos(theta2);

% eom1 = (m1+m2) * l1*theta1'' + m2*l2*theta2''*cos(theta1-theta2) =
% -g*(m1+m2)*sin(theta1) - m2*l2*(theta2')^2 * sin(theta1-theta2)

% eom2 = m2*l1*theta1''*cos(theta1-theta2) + m2*l2*theta2'' =
% -g*m2*sin(theta2) + m2*l1*(theta1')^2 * sin(theta1 - theta2)

% l1 = 2;
% l2 = 1;
% m1 = 1;
% m2 = 2;
% g = -9.81;
% c = cos(theta1-theta2);
% s = sin(theta1-theta2);

% % this is the same storage method used in RK4
% u = [u1 ; u2 ; u3; u4];
% u1 = theta1;
% u2 = theta2;
% u3 = theta1';
% u4 = theta2';

% in matrix vector form:
% [ 1 0 0 0 ; 0 1 0 0 ; 0 0 (m1+m2)*l1 m2*l2*c ; 0 0 m2*l1*c m2*l2 ] *
% [u1' ; u2' ; u3' ; u4'] =
% [u3 ; u4 ; -g*(m1+m2)*sin(u1)-m2*l2*s*u4^2 ; -g*m2*sin(u2)+m2*l1*s*u3^2]

%% the actual code for double pendulum
function rhs = dblPend(t,u)

g = 9.81;
c = cos(u1-u2);
s = sin(u1-u2);

M = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 (m1+m2)*l1 m2*l2*c ; 0 0 m2*l1*c m2*l2 ];
f = [u(3) ; u(4) ; -g*(m1+m2)*sin(u(1))-m2*l2*s*u(4)^2 ; -g*m2*sin(u(2))+m2*l1*s*u(3)^2];

rhs = M\f;

end

[t,u] = ode45(@dblPend,time,u);

%% Plot the results

function status = dblPendPlot(t,u,flag)

if stcmp(flag,'done')
    status = 0;
    return;
end

th1 = u(1,end);
th2 = u(2,end);

x1 = l1*sin(th1);
y1 = l2*cos(th1);

x2 = l1*sin(th1) + l2*sin(th2);
y2 = l1*cos(th1) + l2*cos(th2);

x = [0 x1 x2];
y = [0 y1 y2];

plot(x,y,'k-o');
xlim([-3 3]);
ylim([-3 3]);
daspect([1 1 1]);
pause(0.05);
drawnow;
status=0;

end

opt = odeset('OutputFcn',@dblPendPlot);
[t,u] = ode45(@db