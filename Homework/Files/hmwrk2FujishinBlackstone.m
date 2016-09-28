%% Homework #2
%  Ian Blackstone and Helena-Nikolai Fujishin
%  Math 365, Fall 2016

%% Problems
% Every assignment will include a function like this which will serve as a
% "main" function for your homework. 
% 
function hmwrk2FujishinBlackstone()
%%%

% You will then call all of your homework problem "functions" like this. 
hmwk_problem(@prob1,'prob1');
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
hmwk_problem(@prob5,'prob5');
hmwk_problem(@prob6,'prob6');

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

%% Problem #1 : Computing pi
function prob1()
%Write a MatLab function for approximating pi using Machin's series. 
%(a) Compute an approximation for pi using N=1,2...16 and report results in
%a table
format long
N= 16; 
Sum=0;
Our_Pie = zeros(1,16);
percenterror = zeros(1,16);

for k = 0:N
     Sum = Sum + 4*((-1)^k/(2*k+1))*((4*(1/5)^(2*k+1))-((1/239)^(2*k+1)));
     Our_Pie(k+1)= Sum;
     percenterror(k+1) = abs(((pi-Our_Pie(k+1))/pi))*100;
     fprintf('Our iternation number k: %d\n',k)
     fprintf('Our approximation for pi: %22f\n', Our_Pie(k+1))
     fprintf('Our error for pi: %2f\n', percenterror(k+1))
end
fprintf('Our ending pi: %d\n', Our_Pie(16))
%Plotting the absolute error as a function of the number of terms:
semilogy(percenterror)
title('Graph of percent error as function of N')
xlabel('Number of terms used in sum')
ylabel('Percent Error')

end

%% Problem #2 : Continued Fractions
% 
function prob2()
%(a)Write a function in Matlab that computes the 
%continued fraction expression for given array of numbers a(j), j=1:n.
    function x=contfraction(a)
        n=length(a);
        x=a(n);
        for j = n-1:-1:1
            x=a(j)+1/x;
        end
    end
%(b) Use your function to approximate sqrt(2)
l=20;
a=zeros(1,l);
a(1)=1;
a(2:l)=2;
y = contfraction(a);
error=((sqrt(2)-y)/y)*100;
fprintf('Our approximation for sqrt(2) is %17.16f\n',y)
fprintf('Our percent error is %20.19f\n',error)
%(c) approximate e
a=[2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 14, 1, 1, 16, 1, 1, 18, 1, 1, 20, 1, 1, 22, 1, 1, 24, 1, 1, 26, 1, 1, 28, 1, 1, 30, 1, 1, 32, 1, 1, 34, 1, 1, 36, 1, 1, 38, 1, 1, 40, 1, 1, 42, 1, 1, 44, 1, 1, 46, 1, 1, 48, 1, 1, 50, 1, 1, 52, 1, 1, 54, 1, 1, 56, 1, 1, 58, 1, 1, 60, 1, 1, 62, 1, 1, 64, 1, 1, 66];
y=contfraction(a);
fprintf('Our approximation for e is %17.16f\n',y)
fprintf('We would need n to be at least %1d to obtain full approximation \n',n)
end

%% Problem #3 : 
function prob3()
%Use Newton's Method to calculate sqrt(5)
a=5;
y(1)= 0.5; %this is our initial guess
% iterate using Newton's method to make an estimate of the square root
for n = 1:10
    y = [y,0.5*y(n)*(3-a*y(n)^2)];
end
% Calculate the square root
sq = a*y(n);
% Declare the results of the estimate
fprintf('the square root is approximately %.8f\n',sq)
end

%% Problem #4 : 
function prob4()
%(a)Use fzero to find the friction factor f corresponding to parameter values
%E/D = 0.0001 and ReD = 3 × 10^5. Use a tolerance 10^?8 with fzero. 

% Part a
% Set options for tolerance
options = optimset('TolX',10^-8);

% Create a function to find the zeros of
g = @(f) -2*log10(0.0001/(3.7) + 2.51/(3*10^5 * sqrt(f))) -1/sqrt(f);

% Find the zeroes of the function within a range
fzero(g,[0.001 1],options)

% Part b
% Set options for tolerance
options = optimset('TolX',10^-12);

% Create a function to find the zeros of
f = @(c) -1 + 8*c^2 - 8*c^4 - cos(4*c*sqrt(1-c^2));

% Find the zeroes of the function within a range
fzero(f,[0 0.5],options)
end

%% Problem #5 : 
function prob5()
%Find the minimum depth for a pipe to be buried so it doesn't freeze for 60
%days.
Ti = 20;
Ts = -15;
a = 0.138E-6; %Thermal Conductivity
t = 5.184E6; %t=60days(but this is in seconds)
T = @(x) erf(x/(2*sqrt(a*t)))*(Ti-Ts) + Ts;
min = fzero(T,1);
fprintf('The minimum depth to bury the pipe is %f meters\n',min);
end

%% Problem #6 : 
function prob6()
% Creating a list to hold all operations to be done, all stored as anonymous functions with A,B,C, and D as the inputs.  All four variables must be used as inputs for each item so the for loop can pass the proper inputs to each list item
problemlist = {@(A,B,C,D) 2*A + C, @(A,B,C,D) C - 3*B, @(A,B,C,D) 3*B - 2*D,@(A,B,C,D) A*D, @(A,B,C,D) C*A, @(A,B,C,D) A*C, @(A,B,C,D) B*D, @(A,B,C,D) D*B, @(A,B,C,D) B*C, @(A,B,C,D) C*B, @(A,B,C,D) D*A*B, @(A,B,C,D) 2*D.' + B, @(A,B,C,D) D*D, @(A,B,C,D) A*A, @(A,B,C,D) C.' * D, @(A,B,C,D) B*A.', @(A,B,C,D) -2*A.' + 5*C, @(A,B,C,D) B.' + A*D, @(A,B,C,D) 0.5 * (B + B.'), @(A,B,C,D) 0.5 * (B-B.'), @(A,B,C,D) C*C.',@(A,B,C,D) C.' * C, @(A,B,C,D) (C*C.')^(-1), @(A,B,C,D) (C.'*C)^(-1), @(A,B,C,D) B* (A*D).', @(A,B,C,D) A*D*B.'};

% Generate 4 matrices using random numbers
A = rand(2,3);
B = rand(3);
C = rand(3,2);
D = rand(3);

% Determine the length of the problem set to be done to iterate over the correct number
[~,length] = size(problemlist);

% Iterate over each item in the list.
for n = 1:length
    % Try the nth item in the problem set, catching any errors and outputting the reason for the error.  Problem w throws a warning that is not caught by this code, but it does not prevent the code from running.
    % We know that the results for problem w are innacurate because (C*C.')^-1 * (C*C.') should give us a 3x3 identity matrix and it does not.
    % I could not find a method of catching warnings and errors at the same time and it appears the only way to silence the warning is a user setting.
    % The second elseif term is to catch problem n if it is entered as A^2 rather than A*A.  A^2 gives an error because A is not a square matrix, while A*A warns about inner dimensions not matching .
    try
        problemlist{n}(A,B,C,D);
    catch MExc
        if strcmp(MExc.identifier,'MATLAB:dimagree')
            fprintf('%c cannot be performed, matrix dimensions do not agree \n',n+96)
        elseif strcmp(MExc.identifier,'MATLAB:innerdim')
            fprintf('%c cannot be performed, inner matrix dimensions do not agree \n',n+96)
        elseif strcmp(MExc.identifier,'MATLAB:mpower:notScalarAndSquareMatrix')
            fprintf('%c cannot be performed, inputs are not a scalar and a square matrix \n',n+96)
        else
            fprintf('An unexpected error ocurred when calculating %c',n+96)
        end
    end
end
end