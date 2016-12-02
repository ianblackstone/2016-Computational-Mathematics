%% Homework #7
%  Ian Blackstone, Helena-Nikolai Fujishin
%  Math 365, Fall 2016
%% Problems

function hmwk1()
hmwk_problem(@prob1,'prob1');
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
end
function hmwk_problem(prob,msg)
try
	prob()
	fprintf('% s : Success!\n',msg);
catch me
	fprintf('% s : Something went wrong.\n',msg);
	fprintf('% s\n',me.message);
end
fprintf('\n');
end

%% Problem #1 : Lotka–Volterra predator-prey model NCM 7.16
function prob1()

% For Nonmodified Predator versus prey model:
alpha = 0.01;

% Here we have our system where 
% F(1) = dr/dt, Change in population of rabbits
% F(2) = df/dt, Change in population of foxes
% y(1) = population of rabbits and y(2) = population of foxes
F = @(t,y) [(2.*y(1))-alpha.*(y(1).*y(2)); alpha.*y(1).*y(2)-y(2)];
[t y] = ode45(F,[0 10],[300 150]);
% This is a plot of the population changes over time
figure
plot(t,y(:,1),'b',t,y(:,2),'r') 
legend('Population of rabbits','Population of Foxes')
title('Unmodified Lotka–Volterra predator-prey model') 
xlabel('Time Units')
ylabel('Population')

figure
plot(y(:,2),y(:,1),'g') 
legend('Population of Foxes versus Rabbits')
title('Unmodified Lotka–Volterra predator-prey model') 
xlabel('Fox Population')
ylabel('Rabbit Population')

% For Modified Predator versus prey model:
alpha = 0.01;
R = 400;

% Here we have our system where 
% F(1) = dr/dt, Change in population of rabbits
% F(2) = df/dt, Change in population of foxes
% y(1) = population of rabbits and y(2) = population of foxes
F = @(t,y) [(2.*(1-(y(1)./R)).*y(1))-alpha.*(y(1).*y(2)); alpha.*y(1).*y(2)-y(2)];
[t y] = ode45(F,[0 10],[300 150]);

% This is a plot of the population changes over time
figure
plot(t,y(:,1),'b',t,y(:,2),'r') 
legend('Modified population of rabbits','Modified population of Foxes')
title('Modified Lotka–Volterra predator-prey model') 
xlabel('Time Units')
ylabel('Population')

figure
plot(y(:,2),y(:,1),'g') 
legend('Population of Foxes versus Rabbits')
title('Modified Lotka–Volterra predator-prey model') 
xlabel('Fox Population')
ylabel('Rabbit Population')
end

%% Problem #2 : Trajectory of a spherical cannonball NCM 7.18
function prob2()

% This problem uses several functions stored in external files.

% function f = cannonball0(t,u)

% % Declare constants.
% g = 9.81;
% m = 15;
% c = 0.02;
% p = 1.29;
% s = 0.25;

% % Return the derivative of the vector u.
% f = [u(3);u(4); -c*p*s/(2*m)*u(3)^2; -c*p*s/(2*m)*(u(4))^2-g];

% end

% function f = cannonballneg10(t,u)

% % Declare constants.
% g = 9.81;
% m = 15;
% c = 0.02;
% p = 1.29;
% s = 0.25;

% % Return the derivative of the vector u.
% f = [u(3);u(4); -c*p*s/(2*m)*(u(3)-(-10))^2; -c*p*s/(2*m)*(u(4))^2-g];

% end

% function f = cannonballvar(t,u)

% % Declare constants.
% g = 9.81;
% m = 15;
% c = 0.02;
% p = 1.29;
% s = 0.25;

% % Set the wind to 10 if the time (rounded down) is even.
% if mod(floor(t),2) =  = 0
%     w = 10;
% else
%     w = 0;
% end

% % Return the derivative of the vector u.
% f = [u(3);u(4); -c*p*s/(2*m)*(u(3)-(w))^2; -c*p*s/(2*m)*(u(4))^2-g];

% end

% function f = cannonballrand(t,u)

% % Declare constants.
% g = 9.81;
% m = 15;
% c = 0.02;
% p = 1.29;
% s = 0.25;

% % Return the derivative of the vector u.
% f = [u(3);u(4); -c*p*s/(2*m)*(u(3)-(10*randn))^2; -c*p*s/(2*m)*(u(4))^2-g];

% end

% Declare the initial speed, timespan, and options for the ODE solver.
v = 50;
t1 = linspace(0,10);
opt = odeset('RelTol',10e-2,'event',@ground);

% For our cannon ball without wind:
% I'd like to set different colors to use for legend:
colors = ['y-' 'm.' 'c-' 'r.' 'g-' 'b.' 'k-' 'y.' 'm-' 'c.' 'r-' 'g.' 'b-' 'k.' 'y*' 'mx' 'c+'];
labels = []; % empty matrix for angles
impactspeed = []; % Create empty matrix for my impact speeds to go into
flighttime = []; % Empty matrix for flight times
distance = []; % Empty matrix for downrange distances

% Using my for loop to loop over 17 values of theta.
for i = 1:17

	% Determine theta for this loop.
	theta = i*(pi/36);

	[t u] = ode45(@cannonball0,t1,[0 , 0 , v*cos(theta) , v*sin(theta)],opt); % Calling Ode with initials
	
	% Calculate the final speed and flight time for this loop.
	impactspeed(i) = sqrt(u(end,3)^2+u(end,4)^2);
	flighttime(i) = t(end);

	% Add the data for this loop to the plot.
	distance(i) = u(end,1);
	plot(u(:,1),u(:,2),colors(i))
	axis([0,250,0,140])
	title('Graph of CannonBall without wind (theta in degrees)')
	ylabel('y position')
	xlabel('x position')
	legendInfo{i} = ['theta = ' num2str(i*(pi/36)*(180/pi))];
	labels{i} = ['theta(in degrees) = ' num2str(i*(pi/36)*(180/pi))];
	hold on
end

% Last, set my legend up:
legend(legendInfo)

% Report table of values:
impactspeed = impactspeed.';
xdistance = distance.';
flighttime = flighttime.';
InitialDegrees = labels.';
NoWind = table(impactspeed,xdistance ,flighttime,'RowNames',InitialDegrees)

figure % So we can start a new figure

% For our cannon ball with wind -10m/s:
labels = []; % empty matrix for angles
impactspeed = []; % Create empty matrix for my impact speeds to go into
flighttime = []; % Empty matrix for flight times
distance = []; % Empty matrix for downrange distances

% Using my for loop to loop over 17 values of theta. 
for i = 1:17
	% Set theta for this loop
	theta = i*(pi/36);

	[t u] = ode45(@cannonballneg10,t1,[0 , 0 , v*cos(theta) , v*sin(theta)],opt); % Calling Ode with initials
	
	% Calculate the final speed and time of flight for this loop.
	impactspeed(i) = sqrt(u(end,3)^2+u(end,4)^2);
	flighttime(i) = t(end);

	% Add the data for this loop to the plot
	distance(i) = u(end,1);
	plot(u(:,1),u(:,2),colors(i))
	axis([0,250,0,140])
	title('Graph of CannonBall with wind -10m/s (theta in degrees)')
	ylabel('y position')
	xlabel('x position')
	legendInfo{i} = ['theta = ' num2str(i*(pi/36)*(180/pi))];
	labels{i} = ['theta(in degrees) = ' num2str(i*(pi/36)*(180/pi))];
	hold on
end

% Last, set my legend up:
legend(legendInfo)

% Report table of values:
impactspeed = impactspeed.';
xdistance = distance.';
flighttime = flighttime.';
InitialDegrees = labels.';
Wind10 = table(impactspeed,xdistance ,flighttime,'RowNames',InitialDegrees)

figure % So we can start a new figure

% For our cannon ball with wind 10 m/s if the integer part of t is even, and zero otherwise:
labels = []; % empty matrix for angles
impactspeed = []; % Create empty matrix for my impact speeds to go into
flighttime = []; % Empty matrix for flight times
distance = []; % Empty matrix for downrange distances

% Using my for loop to loop over 17 values of theta. 
for i = 1:17
	% Set theta for this loop
	theta = i*(pi/36);

	[t u] = ode45(@cannonballvar,t1,[0 , 0 , v*cos(theta) , v*sin(theta)],opt); % Calling Ode with initials
	
	% Calculate the final speed and time of flight for this loop.
	impactspeed(i) = sqrt(u(end,3)^2+u(end,4)^2);
	flighttime(i) = t(end);

	% Add the data for this loop to the plot
	distance(i) = u(end,1);
	plot(u(:,1),u(:,2),colors(i))
	axis([0,250,0,140])
	title('Graph of CannonBall with Variable Wind (0 or 10m/s)(theta in degrees)')
	ylabel('y position')
	xlabel('x position')
	legendInfo{i} = ['theta = ' num2str(i*(pi/36)*(180/pi))];
	labels{i} = ['theta(in degrees) = ' num2str(i*(pi/36)*(180/pi))];
	hold on
end

% Last, set my legend up:
legend(legendInfo)

% Report table of values:
impactspeed = impactspeed.';
xdistance = distance.';
flighttime = flighttime.';
InitialDegrees = labels.';
VariableWind = table(impactspeed,xdistance ,flighttime,'RowNames',InitialDegrees)

figure % So we can start a new figure

% For our cannon ball with wind 10 m/s*randn:
labels = []; % empty matrix for angles
impactspeed = []; % Create empty matrix for my impact speeds to go into
flighttime = []; % Empty matrix for flight times
distance = []; % Empty matrix for downrange distances

% Using my for loop to loop over 17 values of theta. 
for i = 1:17
	% Set theta for this loop
	theta = i*(pi/36);

	[t u] = ode45(@cannonballrand,t1,[0 , 0 , v*cos(theta) , v*sin(theta)],opt); % Calling Ode with initials
	
	% Calculate the final speed and time of flight for this loop.
	impactspeed(i) = sqrt(u(end,3)^2+u(end,4)^2);
	flighttime(i) = t(end);

	% Add the data for this loop to the plot
	distance(i) = u(end,1);
	plot(u(:,1),u(:,2),colors(i))
	axis([0,250,0,140])
	title('Graph of CannonBall with Random Wind(theta in degrees)')
	ylabel('y position')
	xlabel('x position')
	legendInfo{i} = ['theta = ' num2str(i*(pi/36)*(180/pi))];
	labels{i} = ['theta(in degrees) = ' num2str(i*(pi/36)*(180/pi))];
	hold on
end

% Last, set my legend up:
legend(legendInfo)

% Report table of values:
impactspeed = impactspeed.';
xdistance = distance.';
flighttime = flighttime.';
InitialDegrees = labels.';
RandomWind = table(impactspeed,xdistance ,flighttime,'RowNames',InitialDegrees)

end

%% Problem #3 : Satellite Orbit
function prob3()

% Part a

% The following function is called from an external file.
% function f = Orbit(t,u)

% mu = 0.012277471;
% mus = 1 - mu;

% D1 = ((u(1)+mu)^2 + u(2)^2)^(3/2);
% D2 = ((u(1)-mus)^2 + u(2)^2)^(3/2);

% f = [u(3); u(4); u(1)+2*u(4)-mus*(u(1)+mu)/D1-mu*(u(1)-mus)/D2; u(2)-2*u(3)-mus*u(2)/D1-mu*u(2)/D2];

% end

% Generate a list of time points for one full cycle.
b = 17.06521656015796;
t = linspace(0,b);

% Initial conditions
u0 = [0.994,0,0,-2.0015851063790825];

% Call the ODE solver using a given relative tolerance and output a phase diagram.
figure
opt = odeset('RelTol',10e-2,'OutputFcn',@odephas2);
[t1,u1] = ode45(@Orbit,t,u0,opt);
title('Graph of satelite position, RelTol = 10E-2')
ylabel('y position')
xlabel('x position')

figure
opt = odeset('RelTol',10e-4,'OutputFcn',@odephas2);
[t2,u2] = ode45(@Orbit,t,u0,opt);
title('Graph of satelite position, RelTol = 10E-4')
ylabel('y position')
xlabel('x position')

figure
opt = odeset('RelTol',10e-6,'OutputFcn',@odephas2);
[t3,u3] = ode45(@Orbit,t,u0,opt);
title('Graph of satelite position, RelTol = 10E-6')
ylabel('y position')
xlabel('x position')

% Part b

% Run the ODE solver for 3 periods using the same relative tolerance.
figure
[t4,u4] = ode45(@Orbit,3*t,u0,opt);
title('Graph of satelite position for 3 orbital periods')
ylabel('y position')
xlabel('x position')

% These results show the satelite to precess in its orbit, which is not the actual behavior.
% Changing this to a larger tolerance causes the orbit to be less table.
% Even the smallest possible tolerance (~2E-14) does not give accurate results over three periods.

% Part c

% Run the 113 ODE solver.
figure
[t5,u5] = ode113(@Orbit,3*t,u0,opt);
title('Graph of satelite position using ODE113 solver')
ylabel('y position')
xlabel('x position')

% These results are much more accurate, though there is still innacuracy in the last portion of the third orbit.

end
