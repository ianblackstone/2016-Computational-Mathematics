%% Problem 1
% part a
x = [-1,0,1];

P2 = @(x) 6 .* (x)/(-1) .* (x - 1)/(-1 - 1) + 1 .* (x + 1)/(1).*(x - 1)/(-1) + 2 .* (x + 1)/(1+1).*(x)/(1)

% part b Not working, don't know why.
% I've tried manually inputting the matrix and generating it with the function and neither works.
A = vander(x)
b = [6;1;2];

a = A\b

P2V = @(x) a(1) + a(2)*x + a(3)*x.^2;

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

plot(xj,fx,'ro',u,b,'b-')
legend('function','interpolation')

% part b

% Not actually sure how to use these two functions, documentation is a bit confusing.
p1 = polyfit(xj,fx,n);

p2 = polyval(p1,fx);


%% Problem 3
% Don't know where to start, skipping for now.


%% Problem 4

% part a

% Seems to work, needs to have the plot properly bounded I think.
x = [-1,0,2,4,5];
y = [0,0,1,2,2+19/22];
t = -2:.01:6;
p = pchip(x,y,t);

plot(x,y,'o',t,p,'-')

% part b

% The format for this is (breaks,coefficients), but I don't know the coefficients for the middle section.
% I think we get them from part a but I'm not sure I did that right so who knows.
p2 = mkpp([0,2,4],[0,0;2,0;19/22,2]);

% part c

x2 = linspace(-1,5);
y2 = ppval(x2,p2);
plot(x2,y2)

%% Problem 5

% part a

load shark;
plot(x,y,'k.-','MarkerSize',10);

% part b

% this gives an error, the first and last x calue are the same and it complains that the data is not a function because of this.
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