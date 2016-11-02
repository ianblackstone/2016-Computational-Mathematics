% Adaptive Quadrature
% if we have some complicated function we will split it into some number of subintervals, but rather than just split into equal
% parts we should instead focus smaller section on more volatile sections of the plot.

% Adapt simpsons rule over the entire interval then compare and subdivide half into a further 2 peces.  Repeat until the error is within tolerance.
% Once you have an accurate enough answer move to the next unintegrated interval and repeat the process of subdividing until each section has error
% below the tolerance.

% There are data structures that make this code simple.  Apply Simpson's rule to the top level choice, then apply it to the next interval seeing each interval as a branching tree:


% 						---->(from a to [a+b/2)
% simp(from a to b) --->  
% 						---->(from (a+b)/2 to b)

% We need to write a function that represents the integrand

L = 10
f = @(x) = 50*x.*exp(-x/4)./(x+5/3);
% integral(function,start,stop,tolerance flag,tolerance value)
Q = integral(f,0,L,'RelTol',1e-10)
% This integral function does adaptive quadrature and Simpson's rule.  You don't have to pass an anonymous function.

%% Numerical integration in higher dimensions
% While adaptive methods can be done in 2D and higher the code becomes complex.
% Divide the domain into a grid, find the integral in 1d on each line of the grid.  Each point in the interior of the grid is weighted by 1,
% each edge is weighted by 1/2, and the corners are weighted by 1/4.

% we want to integrate this function from 0 to 1 on both axes.
f = @(x,y) exp(x-y)
n = 20;

% ndgrid() creates a grid similar to meshgrid.  This function has x on the vertical axis and y on the horizontal.
[xx,yy] = ndgrid(linspace(0,1,n+1));

h = 1/n;
wx = h*[0.5;ones(n-1,1);0.5];
wy = wx';
w = wx*wy;
g = f(xx,yy);
Q = sum(sum(w.*g))
exact = (1-exp(-1)*exp(1)-1)
abs(Q-exact)

% the error in trapazoidal rule is proportional to 1/n^2 , but the total number of points here is not n, but n^2.  So the error is proportional to 1/N, where N is the total number of points.
% As dimensionality increases the error scales to 1/N^(2/d), and this makes increasing the number of points by even a very large amount have rapidly
% diminishing returns.  So to get around this we can use Monte-Carlo integration where N pseudo random points are selected in our volume and the integral is evaluated
% based on these points.

N = 1000;
x = rand(N,2);
Q = sum(f( x(:,1) , x(:,2) ));
abs(Q-exact)
plot( x(:,1) , x(:,2) )

% Because this is random we get some variability in the answers, with errors ranging from 0.2 to 6E-6.  To help smooth out this we can use sets of numbers that cover our
% entire domain rather than allowing the numbers to clump up which throws off our results.

N = 1000;
x = net(haltonset(2),N);
Q = sum(f( x(:,1) , x(:,2) ));
abs(Q-exact)
plot( x(:,1) , x(:,2) )

%  This uses the Halton set of quasi random numbers to fill the entire domain with even distribution.  These random points always have the error scale with 1/N^(1/2), meaning that anything higher than 4 dimensions benefits from using them.