%% Homework #1
%  Ian Blackstone
%  Math 365, Fall 2016

%% Problems
% This file shows you the prefered way to format your homework for this 
% course using the Matlab publish command.

function hmwk1()
%%% 
% Every assignment will include a function like this which will serve as a
% "main" function for your homework. 
% 

% You will then call all of your homework problem "functions" like this. 
hmwk_problem(@prob1,'prob1');
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
hmwk_problem(@prob5,'prob5');
hmwk_problem(@prob6,'prob6');
hmwk_problem(@prob7,'prob7');
hmwk_problem(@prob8,'prob8');

end

function hmwk_problem(prob,msg)
%%%
% This function should be included in every assignment
try
    prob()
    fprintf('%s : Success!\n',msg);
catch me
    fprintf('%s : Something went wrong.\n',msg);
    fprintf('%s\n',me.message);
end
fprintf('\n');
end

%% Problem #1 : Surface area of a torus
% In this problem, we compute the surface area of a torus whose
% inner radius is 2.21 and whose outer radius is 4.36. 
% The result is included in the report.
function prob1()

% Declaring constants
R = 4.36;
r = 2.21;

% Calculate an display surface area
S = pi^2 * (R - r) * (R + r)

end

%% Problem #2 : Plotting a simple relationship
% In this problem, we plot the relationship between energy and magnitude of
% an earthquake according to the Richter scale.  The plot is included in
% the report.
function prob2()

% Create two lists, one from 1 through 9 an the second an empty list.
M = 1:9;
E = 10^4.4 * 10.^(3/2 * M);

% Plot using a semi logarthmic plot.
semilogy(M,E)

end

%% Problem #3 : Simple vectorization
% Compute the squared sum of a vector of entries.
function prob3()

%%%
% Create the arrays and then include the code for computing the sum of the
% squares of the entries here.

% Generate two anonymous functions of x.
X = @(x) x;
Y = @(x) 5-2*x;

% Create a blank list for values of Z.
Z = [];

% Fin the values of Z.
for I = 1:5
	Z = [Z,(X(I)+Y(I))^2];
end

% Print all values of Z.
Z

% This section generates Z in a different way without using a for loop by calling a list of values as our variable.

% Create a list of values for I.
I = [1:5];

% Generate values of Z.
Z = (X(I)+Y(I)).^2;

% Display all values stored in Z.
Z

end


%% Problem #4 : Simple anonymous function handle; plotting
function prob4()

%%%
% Create function handles for two functions and
% construct a third composite function. 

% Anonymous function handles go here

% Construct a vector of equally spaced points over [-3,3] 

% Plot the results

% Generate 3 anonymous functions.
f = @(x) cos(2*x);
g = @(x) exp(x);
h = @(x) g(f(g(x)));

% Create a list of points to plot over.
x = linspace(-3,3,500);

% Define a figure and plot it.
figure
plot(h(x))
title('Plot of h(x) from -3 to 3')
xlabel('X')
ylabel('h(x)')

% display the value of h at x = 4.3.
h(4.3)

end

%% Problem #5 : Loading data from a file, simple statistics, and using fprintf
function prob5()

%%%
% Here is how you load data from the file heights.dat.  Note you first need
% to download this file from the course website and save it in the working
% directory for this homewokr.
h = load('heights.dat');

% Compute the min, max, mean, and standard deviation.

% Print the results using fprintf

% open the .dat file and read the data into a variable.
f = fopen('heights.dat','r');
h = fscanf(f,'%f');

% Find the minimum, maximum an standard deviation of the data.
Max = max(h);
Min = min(h);
Sdev = std(h);

% Print the values of the maximum, minimum, and standard deviation.
fprintf('the maximum is %.2f \n',Max)
fprintf('the minimum is %.2f \n',Min)
fprintf('the standard deviation is %.4f \n',Sdev)

end

%% Problem #6 : Generating sequences of numbers with a for loop
% Use matrix multiplication to generate a famous sequence from mathematics.
function prob6()

%% Part a
% Declare two arrays.
A = [1,1;1,0];
X = [1,0;0,1];

% Generate values for X2 through X6 and display them.
for k = 1:5
	X = A * X;
	disp(['X(' num2str(k+1) ') = '])
	disp(X)
end

%% Part b


%% Part c

% Declare variables
A = [1,1;1,0];
X = [1,0;0,1];

for k = 1:29
	X = A * X;
end

disp(['X(30) = '])
disp(X)

end

%% Problem #7 : Quadratic formula
% NCM, problem 1.38. The quadratic formula is relatively simple, but to
% implement it properly on the computer in floating point arithmetic takes
% some care.

function prob7()

% Generate two anonymous functions to find each root.
X1 = @(a,b,c) (-b+sqrt(b^2 -4*a*c))/(2*a);
X2 = @(a,b,c) (-b-sqrt(b^2 -4*a*c))/(2*a);

% Declare variables
a = 1;
b = -100000000;
c = 1;

% Find values for each root.
X1(a,b,c)
X2(a,b,c)

% find roots using the roots function.
Xr = roots([a b c]);

Xr(1) * Xr(2)

X1(a,b,c) * X2(a,b,c)

c/a

end

%% Problem #6 : How to succeed in Math 365
function prob8()

%%%
% Publish allows you to create lists, use different font styles, and include
% preformatted code
%
% *How to succeed in Math 365*
% 
% * _Always_ start your homework early
% * _Don't_ spend too much time googling for answers
% * _Read_ the 
% <http://math.boisestate.edu/~wright/course/m365//homework_tips.pdf 
% homework tips>!
% 
% *Steps for getting help on homework problems.*
% 
% # Read the Matlab tutorials available on the course website
% # Read lecture notes and demo codes on the online website. 
% # Use Matlab online "help" system for help on Matlab commands. 
% # Read the <http://www.mathworks.com/moler/chapters.html Course textbook>
% # Email the professor for help, if you can't find answers in the above. 
% # Do not spend too much time with Prof. Google or Dr. YouTube.  This 
% is likely going to be a waste of time!  Spend more time thinking about
% what you have learned in class, and debugging your own code.

%%%
% Include sample code that you don't want run by "formatting" the code
% like this.  Use exactly three spaces between the percent sign and 
% your code.
% 
%   curly = 4*pi;
%   larry = sin(curly);
%   moe = tan(curly + larry);
%   

%%%
% There are lots of helpful hints for publishing by issuing the command
% 
%   >> doc publishing markup

end
