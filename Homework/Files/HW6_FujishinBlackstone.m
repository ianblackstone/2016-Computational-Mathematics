%% Homework #6
%  Ian Blackstone, Helena-Nikolai Fujishin
%  Math 365, Fall 2016
%% Problems
function hmwk1()
hmwk_problem(@prob1,'prob1');
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
hmwk_problem(@prob5,'prob5');
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

%% Problem #1 : Degree of Precision
function prob1()
% Part A
a=0;
b=1;
n = 100;
f = @(x,k) x^k;
exct = @(x,k) 1/(k+1)*x^(k-1);

for k = 0:n
T = (b-a)/2*(f(a,k) + f(b,k));
err = abs(exct(b,k)-T);
if err >= 1E-14
	fprintf('m for part a = % .i \n',k-1)
	break
end
end

% Part B: Simpson's rule
c =(a+b)/2;
h =(b-a);
for k = 0:n
    M = h*f((a+b)/2,k);
    T = h*(f(a,k)+f(b,k))/2;
S = ((2/3)*M)+((1/3)*T);
err = abs(exct(b,k)-S);
if err >= 1E-14
	fprintf('m for part b = % .i \n',k-1)
	break
end
end

end

%% Problem #2 :
function prob2()
% Our given/Known Values
xi = [-0.96816023950763 -0.83603110732664 -0.61337143270059  -0.32425342340381 0 0.32425342340381 0.61337143270059 0.83603110732664 0.96816023950763]';
ci = [0.081274388361574 0.18064816069486 0.26061069640294 0.31234707704 0.33023935500126 0.31234707704 0.26061069640294 0.18064816069486 0.081274388361574]';
n = 9; % Our iteration number

% Part A f1
f1= @(x) ((x-1).^2).*exp(-(x.^2))
I1 = integral(f1,-1,1); % Exact integral of f1
GQForm1=0; % Start our initial use of GQ formula
for i =1:n
GQForm1 = GQForm1+(ci(i)*f1(xi(i))); % Summation of the GQ formula
abserror1(i) = abs(I1-GQForm1);
end
iteration = [1:n]';
error = abserror1';
UsingGQf1 = table(iteration,xi,ci,error)

% Part A f2
f2 = @(x) 2*(1./(1+x.^2))
I2 = integral(f2,-1,1);
GQForm2=0;
for i =1:n
GQForm2 = GQForm2+(ci(i)*f2(xi(i))); % Summation of the GQ formula
abserror2(i) = abs(I2-GQForm2);
end
iteration = [1:n]';
error = abserror2';
UsingGQf2 = table(iteration,xi,ci,error)

% Part B
% Let us set-up the Simpson's Rule:
S1 = simps(-1,1,8,f1);
errorS1=abs(I1-S1)
S2 = simps(-1,1,8,f2);
errorS2=abs(I2-S2)
% To compare errors side by side:
errors1 = [errorS1;abserror1(end)];
errors2 = [errorS2;abserror2(end)];
Array = {'Simpsons Rule','GQ Formula'};
T=table(errors1,errors2,'RowNames',Array)
% From this we can compare and see the GQFormula is much more accurate than
% Simpsons

% Part C
% Determine degree of precision of 9-point GQ formula of f1
a = -1;
b = 1;
n = 9;
f1 = @(x,k) ((x-1).^2).*exp(-(x.^2)).^k
for k = 0:n
T = (b-a)/2*(f1(a,k) + f1(b,k));
exct = @(x) integral(@(x) f1(x,k),a,b); % Since the function alone won't integrate, I have to make bypass function
err = abs(exct(b)-T);
if err >= 1E-14
	fprintf('m for part f1 = % .i \n',k-1)
	break
end
end

% Determine degree of precision of 9-point GQ formula of f2
a = -1;
b = 1;
n = 9;
f2 = @(x,k) (2*(1./(1+x.^2))).^k;
for k = 0:n
T = (b-a)/2*(f2(a,k) + f2(b,k));
exct = @(x) integral(@(x) f2(x,k),a,b); % Since the function alone won't integrate, I have to make bypass function
err = abs(exct(b)-T);
if err >= 1E-14
	fprintf('m for part f2 = \n',k-1)
	break
end
end

end
    
    
%% Problem #3 : 
function prob3()

% Part a
% Define parameters and conditions.
N = 1000;
h = 1/N;
j = 0:N;
a = 293.15;
b = 293.15;
B = 1e-3;
S = 1e4;
Cp = 200;
rho = 10;
eps = 5e-2;

% Generate a list of zeros and set the first and last values to our end point conditions.
uj = zeros(N,1);
uj(1) = a;
uj(end) = b;

% Define a function for each integral in u.
f1 = @(x) x.*(-S/(rho*Cp*B)).*(exp(-(((x-(1/2))/(eps)).^2)));
f2 = @(x) (x-1).*(-S/(rho*Cp*B)).*(exp(-(((x-(1/2))/(eps)).^2)));

% Define a function to integrate over to find the temperature at that point.
u = @(x) (a-quad(f1,0,x))*(1-x) + (b+quad(f2,x,1))*x;

% Loop over each value of x to find the temperature at that point.
for x = j
	uj(x+1) = u(h.*x);
end

% Plot the data.
plot(h.*j,uj)

% Part b

% Define a function that will be zero at the point we will burn ourselves.
u2 = @(x) (a-quad(f1,0,x))*(1-x) + (b+quad(f2,x,1))*x - 343.15;

% Use fzero and an initial guess of 0.25 to find a zero.
burn = fzero(u2,0.25)

end


%% Problem #4 :
function prob4()

% Part a
% The exact answer is -512/pi^6 as given by the following code:
% syms a b c d e f
% exact = int(int(int(int(int(int(sin(pi/2*(a + b + c + d + e + f)),a,0,1),b,0,1),c,0,1),d,0,1),e,0,1),f,0,1)

% Part b
N = 10000;
rng(11032016);
MCpoints = rand(N,6);

g = @(a,b,c,d,e,f) sin(pi/2*(a + b + c + d + e + f))

a = MCpoints(:,1);
b = MCpoints(:,2);
c = MCpoints(:,3);
d = MCpoints(:,4);
e = MCpoints(:,5);
f = MCpoints(:,6);

MC = sum(g(a,b,c,d,e,f))/N

% Part c
HSpoints = net(haltonset(6),N);

a = HSpoints(:,1);
b = HSpoints(:,2);
c = HSpoints(:,3);
d = HSpoints(:,4);
e = HSpoints(:,5);
f = HSpoints(:,6);

HS = sum(g(a,b,c,d,e,f))/N

% Part d

n = 10;
h = 1/n;
d = 6;
% Trapezoidal rule weights in one dimension
w1 = h*[0.5;ones(n-1,1);0.5];
w = w1;
% Generate the trapezoidal rule weights in 6 dimensions
for j=1:d-1
w = bsxfun(@times,w,reshape(w1,[ones(1,j) n+1]));
end
% Grid of points
[x1,x2,x3,x4,x5,x6] = ndgrid(linspace(0,1,n+1));
g = w.*f(x1,x2,x3,x4,x5,x6);
Qtrap = sum(g(:));

% Part e



end

%% Problem #5 : 
function prob5()
% Determine duration and step size.
Tmax = 20;
h = 0.1;

% Loop to generate three plots.
for x = 0:2
	% Initial fox starting location.
	F = [4*rand() - 4,4*rand() - 4];

	% create two lists to hold each set of points.
	Rpoints = zeros(Tmax/h,2);
	Fpoints = zeros(Tmax/h,2);

	% calculate the path in steps of h.
	for t = 0:Tmax/h
		R = [cos(t*h),sin(t*h)];
		dR = [-sin(t*h),cos(t*h)];
		dF = norm(dR,2)*(R-F)/norm(R-F,2);
		F = F + h*dF;
		Rpoints(t+1,:) = R;
		Fpoints(t+1,:) = F;
	end

	% Plot the results.
	figure()
	plot(Rpoints(:,1),Rpoints(:,2),'r-+',Fpoints(:,1),Fpoints(:,2),'k-o')
end

end
