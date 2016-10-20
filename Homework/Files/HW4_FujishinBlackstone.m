%% Homework #4
%  Ian Blackstone
%  Math 365, Fall 2016

%% Problems

function hmwk1()

hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
hmwk_problem(@prob5,'prob5');
hmwk_problem(@prob6,'prob6');


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

%% Problem #1 : Polynomial Interpolation
function prob1()

% part a, given values of (-1,6),(0,1) and (1,2) determine unique 2nd degree
% polynomial p2(x) to the data.
p2 = @(x) 3.*x.*(x-1)-(x+1).*(x-1)+2.*(x+1)/2.*x;
x = [-1:1:1];

% part b) we can approach using Van Der Monde matricies where we have three
%equations and three unknowns to solve for our coefficients
A = vander(x);
B = [6;1;2];

coeff = A\B;
a0 = coeff(3)
a1 = coeff(2)
a2 = coeff(1)

p2new = @(x) a0+a1.*x+ a2.*x.^2

% part c: check that a and b give same result:
X = linspace(-1,1,20);

err = p2(X) - p2new(X);
fprintf('The maximum error is %.3f \n',max(err))

end

%% Problem #2 : Barycentric Interpolation
function prob2()

% part a
% create 52 points to evaluate the function at
n = 52;
xj = -cos((0:n)*pi/n);

% Delcare the function and create y values at the desired points
f = @(x) abs(x-0.2)+(1-x.^2).^3;
fj = f(xj);

% Create 1000 points to evaluate the interpolation at.
m = 1000; 
u = -cos(pi*(0:m)/m);

% Create the barycentered interpolation of the function.
qn = baryinterp(xj,fj,u);

% Plot the function and the interpolation of the function
figure
plot(u,qn,'r-',xj,fj,'s')
legend('function','interpolation')
title('Barycentric Interpolation')

%part b
% create an interpolation of the function using the polyfit function and evaluate the polynomial at the points created earlier using polyval
cj = polyfit(xj,fj,n);
pn = polyval(cj,u);

% plot the function and the interpolation
figure
plot(u,pn,'r-',xj,fj,'s')
end

%% Problem #3 : Approximation Theory
function prob3()

% List the coeficients of the 14th degree taylor expansion of e^(2x)
T = [1,2,2,4/3,2/3,4/15,4/45,88/315,2/315,4/2835,4/14175,4/467775,8/6081075,8/42567525];

% Create a polynomial using the coefficients given and define the exact function.
x = linspace(-1,1,201);
p = polyval(T,x);
y = exp(2*x);

% calculate the difference between the polynomial and the exact function and plot the error.
err = p - y;
n1 = norm(err,1);
n2 = norm(err,2);
ninf = norm(err,inf);
figure
table(n1,n2,ninf)
plot(x,err)
title('Taylor expansion of e^(2x)')

% part b

% Define equally spaced points, Chebyshev points, and legendre points within the range of -1 to 1.
i = 0:14;
xeq = -1 + 2*i./14;
xch = -cos(i*pi()./14);
xle = [-0.987992518020485 -0.394151347077563 0.570972172608539 -0.937273392400706 -0.201194093997435 0.724417731360170 -0.848206583410427 0 0.848206583410427 -0.724417731360170 0.201194093997435 0.937273392400706 -0.570972172608539 0.394151347077563 0.987992518020485];

% calculate the y values for each set of x and create polynomials to interpolate the points.
Teq = polyfit(xeq,exp(2*xeq),14)
Tch = polyfit(xch,exp(2*xch),14)
Tle = polyfit(xle,exp(2*xle),14)

% evaluate the created polynomials at each set of points and 
peq = polyval(Teq,x);
pch = polyval(Tch,x);
ple = polyval(Tle,x);

% plot the error for each method and report the norms in a table.
figure
erreq = peq - y;
n1eq = norm(erreq,1);
n2eq = norm(erreq,2);
ninfeq = norm(erreq,inf);
table(n1eq,n2eq,ninfeq)
plot(x,erreq)
title('Equally spaced points')

figure
errch = pch - y;
n1ch = norm(errch,1);
n2ch = norm(errch,2);
ninfch = norm(errch,inf);
table(n1ch,n2ch,ninfch)
plot(x,errch)
title('Chebyshev points')

figure
errle = ple - y;
n1le = norm(errle,1);
n2le = norm(errle,2);
ninfle = norm(errle,inf);
table(n1le,n2le,ninfle)
plot(x,errle)
title('Legendre points')

% part c
% The chebyshev and legendre points produce very small errors, legendre has a smaller norm1 and 2 but larger infinite norm.  I would select the chebyshev points to use however, as they were able to be built using a simple formula rather than typed manually.

% part d
% call the chebfun function to define the function to be interpolated and the remez function to find the best interpolation with 14 degrees.
f = chebfun('exp(2*x)');
[r,errr] = remez(f,14);

% Plot the error of both the Chebyshev point interpolation and the chebfun interpolation.
figure
plot(x,errch,'b',x,errr*x,'r')
title('error in Chebyshev points versus chebfun interpolation')
legend('Chebyshev','Chebfun')


end


%% Problem #4 : Piecewise Cubic Hermite Interpolation
function prob4()


end


%% Problem #5 : Closed Curve Interpolation
function prob5()
%part a
load shark;
plot(x,y,'k.-','MarkerSize',10);
N=100
t=1:length(x);
tt=linspace(t(1),t(end),101);
%part b
%Using Spline: create smooth parameteric line
xx=spline(t,x,tt);
yy=spline(t,y,tt);
px=pchip(t,x,tt);
py=pchip(t,y,tt);
plot(x,y,'bs',xx,yy,'k.-',px,py,'b-')
legend('Data Points','Spline Interpolation','PCHIP Interpolation')
%part c
x1=[x(end-1);x;x(2)];
y1=[y(end-1);y;y(2)];
t1=1:length(x1);
tt1=linspace(t1(1),t1(end),N+1);
xx1=spline(t1,x1,tt1);
yy1=spline(t1,y1,tt1);
px1=pchip(t1,x1,tt1);
py1=pchip(t1,y1,tt1);
plot(x1,y1,'bs',px1,py1,'r-',xx1,yy1,'g-')
legend('Data Points','Corrected PCHIP Interpolation','Corrected Spline Interpolation')
%part d
%The representation I like best is that of the corrected PCHIP
%Interpolation. The shark's curves are much smoother and there is less of
%an overshoot on the nose.
%part e

%I don't know how to do this. This is ridiculous and all sources I have
%tried have yeilded no results including NCM as well as MathWorks forums
%and even other students. I have spent about 2 hours on this question alone
%and speaking with the instructor doesn't help because he makes me feel
%stupid and offers only minimal support. 

end

%% Problem #6 : Scattered Data Interpolation
function prob6()
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
end