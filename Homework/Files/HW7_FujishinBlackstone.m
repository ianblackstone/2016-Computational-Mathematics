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
    fprintf('%s : Success!\n',msg);
catch me
    fprintf('%s : Something went wrong.\n',msg);
    fprintf('%s\n',me.message);
end
fprintf('\n');
end

%% Problem #1 : Lotka–Volterra predator-prey model NCM 7.16
function prob1()
%For Nonmodified Predator versus prey model:
alpha = 0.01;
%Here we have our system where 
%F(1)=dr/dt, Change in population of rabbits
%F(2)=df/dt, Change in population of foxes
%y(1)=population of rabbits and y(2)=population of foxes
F=@(t,y) [(2.*y(1))-alpha.*(y(1).*y(2)); alpha.*y(1).*y(2)-y(2)];
[t y]=ode45(F,[0 10],[300 150]);
%This is a plot of the population changes over time
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

%For Modified Predator versus prey model:
alpha = 0.01;
R = 400;
%Here we have our system where 
%F(1)=dr/dt, Change in population of rabbits
%F(2)=df/dt, Change in population of foxes
%y(1)=population of rabbits and y(2)=population of foxes
F=@(t,y) [(2.*(1-(y(1)./R)).*y(1))-alpha.*(y(1).*y(2)); alpha.*y(1).*y(2)-y(2)];
[t y]=ode45(F,[0 10],[300 150]);
%This is a plot of the population changes over time
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
%Defining our constants:
g = 9.81;
m = 15;
v =50;
c = 0.02;
rho = 1.29;
s = 0.25;

for k = 1:17
    sprime= @(t,D) [v.*cos(k*pi/36);v.*sin(k*pi/36); (-g./v).*cos(k*pi/36); (-D./m-g).*sin(k*pi/36)]
    D = @(t,sprime) (c*rho*s/2)*((sprime(1)-0).^2+sprime(2).^2)
    [t sprime]=ode45(D,[0 17],[0]);
end

end
    
    
%% Problem #3 : 
function prob3()

b = 17.06521656015796;


t = linspace(0,b);

pos0 = [0.994,0,0,-2.0015851063790825];

function f = Orbit(t,pos)

mu = 0.012277471;
mus = 1 - mu;

D1 = ((pos(1)+mu)^2 + pos(2)^2)^(3/2);
D2 = ((pos(1)-mu)^2 + pos(2)^2)^(3/2);

f = [pos(3);pos(4);pos(1)+2*pos(4)-mus*(pos(1)+mu)/D1-mu*(pos(1)-mus)/D2;pos(2)-2*pos(3)-mus*pos(2)/D1-mu*pos(2)/D2];
end

[t,u] = ode45(@Orbit,t,pos0);

% x3' = x1 + 2x4 - mu* * (x+mu)/D1 - mu (x-mu*)/D2


end
