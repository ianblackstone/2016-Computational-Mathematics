%% Problem 1
% part a
x = [-1,0,1];

P2 = @(x) 6 .* (x)/(-1) .* (x - 1)/(-1 - 1) + 1 .* (x + 1)/(1).*(x - 1)/(-1) + 2 .* (x + 1)/(1+1).*(x)/(1)

% part b
% New! Improved! Now even gives correct answers!
A = vander(x)
b = [6;1;2];

a = A\b

P2V = @(x) a(1)*x.^2 + a(2)*x + a(3);


% part c
X = linspace(-1,1,20);

err = P2V(X) - P2(X);

plot(X,err)

%% Problem 2

% part a

n = 52;
j = [0:n];
xj = -cos(j.*pi()/n);

fx = abs(xj - 0.2) + (1-xj.^2).^3;

m = 1000;
i = [0:m];
u = -cos(pi()*i./m);

b = baryinterp(xj,fx,u);

figure
plot(xj,fx,'ro',u,b,'b-')
legend('function','interpolation')

% part b

p1 = polyfit(xj,fx,n);

p2 = polyval(p1,u);

figure
plot(xj,fx,'ro',u,p2,'b-')
legend('function','interpolation')

%% Problem 3

T = [1,2,2,4/3,2/3,4/15,4/45,88/315,2/315,4/2835,4/14175,4/467775,8/6081075,8/42567525];

x = linspace(-1,1,201);
p = polyval(T,x);
y = exp(2*x);

err = p - y;
n1 = norm(err,1);
n2 = norm(err,2);
ninf = norm(err,inf);
figure('name','Taylor expansion')
table(n1,n2,ninf)
plot(x,err)

% part b

i = 0:14;

xeq = -1 + 2*i./14;
xch = -cos(i*pi()./14);
xle = [-0.987992518020485 -0.394151347077563 0.570972172608539 -0.937273392400706 -0.201194093997435 0.724417731360170 -0.848206583410427 0 0.848206583410427 -0.724417731360170 0.201194093997435 0.937273392400706 -0.570972172608539 0.394151347077563 0.987992518020485];

Teq = polyfit(xeq,exp(2*xeq),14)
Tch = polyfit(xch,exp(2*xch),14)
Tle = polyfit(xle,exp(2*xle),14)

peq = polyval(Teq,xeq);
pch = polyval(Tch,xch);
ple = polyval(Tle,xle);
yeq = exp(2*xeq);
ych = exp(2*xch);
yle = exp(2*xle);

% I can't find what I'm doing wrong but these don't look right.

figure('name','equally spaced points')
erreq = peq - yeq;
n1eq = norm(erreq,1);
n2eq = norm(erreq,2);
ninfeq = norm(erreq,inf);
table(n1eq,n2eq,ninfeq)
plot(xeq,erreq)

figure('name','Chebyshev points')
errch = pch - ych;
n1ch = norm(errch,1);
n2ch = norm(errch,2);
ninfch = norm(errch,inf);
table(n1ch,n2ch,ninfch)
plot(xch,errch)

figure('name','Legendre points')
errle = ple - yle;
n1le = norm(errle,1);
n2le = norm(errle,2);
ninfle = norm(errle,inf);
table(n1le,n2le,ninfle)
plot(xle,errle)

% part c
% The chebyshev points produce the smallest error of all the points.

% part d
f = chebfun('exp(2*x)');
r = remez(f,14);


%% Problem 4

% part a

% Doesn't work
% x = [-1,0,2,4,5];
% y = [0,0,1,2,2+19/22];
% t = -2:.01:6;
% p = pchip(x,y,t);

% % Still not working.
% x1 = [-1 0];
% y1 = [0 0];

% p1 = polyfit(x1,y1,1);

% x2 = [0 2];
% y2 = [0 1];

% p2 = polyfit(x2,y2,1);

% x3 = [2 4];
% y3 = [1 2];

% p3 = polyfit(x3,y3,1);

% x4 = [4 5];
% y4 = [2 (2+19/22)];

% p4 = polyfit(x4,y4,1);

% p = [p1;p2;p3;p4];
% y = [0,0,1,2,2,2];
% x = [-1,0,2,4,5,6];

% d = [(y(2)-y(1))/(x(2)-x(1)), (y(3)-y(2))/(x(3)-x(2)), (y(4)-y(3))/(x(4)-x(3)), (y(5)-y(4))/(x(5)-x(4)),0,0];
% c = @(kn,k,kp) 3*k -2*kn - kp;
% b = @(kn,k,kp) kn - 2*k + kp;

% g1 = [y(1), d(2), c(d(1),d(2),d(3)), b(d(1),d(2),d(3))];
% g2 = [y(2), d(3), c(d(2),d(3),d(4))/2, b(d(2),d(3),d(4))/4];
% g3 = [y(3), d(4), c(d(3),d(4),d(5))/2, b(d(3),d(4),d(5))/4];
% g4 = [y(4), d(5), c(d(4),d(5),d(6)), b(d(4),d(5),d(6))];

% g = [g1;g2;g3;g4];


% Declare x and y points
y = [0,0,1,2,2];
x = [-1,0,2,4,5];

% declare the values of d, c, and b for each interval
d = [0, 0 , 0, 19/22 , 0];
c = [0 , 0 , (1.5-19/22)/2 , (1.5-19/11)/2 , 0];
b = [0 , 0.5 , (-1+19/22)/4 , (19/22-1)/4 , 0];

% create a matrix of each value
g = [y(1:5)' d' c' b'];

% create a piecewise function
P = mkpp(x,g);

% Evaluate the piecewise function from -1 to 5 and plot it
x2 = linspace(-1,5);
y2 = ppval(x2,P);
plot(x2,y2)

% These results are not correct, I do not understand how the mkpp function is working to create these results.

% d = zeros(4,1);
% c = zeros(4,1);
% b = zeros(4,1);

% % for i = 1:4
% % 	d(i) = (y(i+1) - y(i))/(x(i+1)-x(i));
% % end

% d = [0 19/22 19/22 0];

 for i = 1:3
 	c(i) = ( d(i) - d(i+1) )/( x(i+1) - x(i) );
 	b(i) = ( -d(i) + d(i+1) )/( x(i+1) - x(i) )^2;
 end


% c = [-0.5 0.5 0.25 0];
% b = [0.5 -0.5/4 -0.25 0.5]

g = [y(2:5)' d' c' b']


% part b

% The format for this is (breaks,coefficients), but I don't know the coefficients for the middle section.
% I think we get them from part a but I'm not sure I did that right so who knows.
P = mkpp([-1,0,2,4],g);
% make a matrix with [f0,d0,c0,b0;f1,d1,d2,c2,b2;...]

% part c

x2 = linspace(-1,5);
y2 = ppval(x2,P);
plot(x2,y2)

%% Problem 5

% part a

% see lab

load shark;
plot(x,y,'k.-','MarkerSize',10);

% part b

% this gives an error, the first and last x value are the same and it complains that the data is not a function because of this.
t = 0:0.01:1;
p3 = pchip(x,y,t);

% part c
% I think this will be solved by appending the second point to the end of the list so the curves overlap.
% Can't test that until I get part b working.


% part d



% part e



%% Problem 6

% Load data
load topodata;

% First figure, color coded scatter plot fo the data
figure
scatter3(x,y,z,40,z,'.');
colorbar

% call rbffit function
r = rbffit(x,y,z);

% generate a list of X and Y points
x1 = linspace(0,261,100);
y1 = linspace(0,399,153);

% Create a grid from the generated points
[X,Y] = meshgrid(x1,y1);

% Call rbval function to generate the matching Z points for the x and y points generated earlier.
Z = rbfval(r,x,y,X,Y);

% Second figure, a 3D surface map
figure
surf(X,Y,Z)
shading interp;
colormap(autumn);
lighting phong;
camlight right;

% define the contour lines to use.
v = 0:20:260;

% third figure, a contour plot.
figure
contour(X,Y,Z,v)