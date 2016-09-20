%% Problem 1
% initialize a variable and a list of zeros
X = 0;
pie = zeros(1,16);

% Sum each term up to 16 terms, storing each succesive approximation in the pie list
for k = 0:16
	X = X + (-1)^k / (2 * k + 1) * ((4 * (1 / 5)^(2 * k + 1) - (1/239)^(2 * k + 1)));
    pie(k+1) = X;
end

% Multiply every approximation by 4.
pie = 4*pie;

% Create a new array with the error between pi and our approximation
pieerror = abs(pie/pi - 1)*100;

% Plot the error on a logarithmic scale
semilogy(pieerror)


%% Problem 2
% Declare the number of items to sum
n = 50;

% Create a list of zeros and set the last item in that list equal to 1/n
X = zeros(1,n);
X(end) = 1/n;

% Starting at the first nonzero item in the list (n-k+1) generate the next item in the list
for k = 1:n
    X(n-k) = n-k+1 + 1/X(n-k+1);
end

%% Problem 3
% Declare a, the number we will be taking the square root of
a = 5;

% Declare an initial guess 
y = [0.5];

% Iterate using Newton's method to make an estimate of the square root
for n = 1:10
    y = [y,0.5*y(n)*(3-a*y(n)^2)];
end

% Calculate the square root
sq = a*y(n);

% Declare the results of the estimate
fprintf('the square root is approximately %.8f\n',sq)


%% Problem 4

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


%% Problem 5

% Declare constants
Ti = 20;
Ts = -15;
a = 0.138E-6;
t = 5.184E6;

% Create a function for temperature
T = @(x) erf(x/(2*sqrt(a*t)))*(Ti-Ts) + Ts;

% Calculate the minimum depth
mindepth = fzero(T,10);

% Output the minimum depth
fprintf('The minimum depth to bury the pipe is %f meters\n',mindepth);


%% Problem 6

% Creating a list to hold all operations to be done, all stored as anonymous functions with A,B,C, and D as the inputs.  All four variables must be used as inputs for each item so the for loop can pass the proper inputs to each list item
problemlist = {@(A,B,C,D) 2*A + C, @(A,B,C,D) C - 3*B, @(A,B,C,D) 3*B - 2*D,@(A,B,C,D) A*D, @(A,B,C,D) C*A, @(A,B,C,D) A*C, @(A,B,C,D) B*D, @(A,B,C,D) D*B, @(A,B,C,D) B*C, @(A,B,C,D) C*B, @(A,B,C,D) D*A*B, @(A,B,C,D) 2*D.' + B, @(A,B,C,D) D*D, @(A,B,C,D) A*A, @(A,B,C,D) C.' * D, @(A,B,C,D) B*A.', @(A,B,C,D) -2*A.' + 5*C, @(A,B,C,D) B.' + A*D, @(A,B,C,D) 0.5 * (B + B.'), @(A,B,C,D) 0.5 * (B-B.'), @(A,B,C,D) C*C.',@(A,B,C,D) C.' * C, @(A,B,C,D) (C*C.')^(-1), @(A,B,C,D) (C.'*C)^(-1), @(A,B,C,D) B* (A*D).', @(A,B,C,D) A*D*B.'};

% Generate 4 matrices using random numbers
A = rand(2,3);
B = rand(3);
C = rand(3,2);
D = rand(3);

% Determine the length of the problem set to be done to iterate over the correct number
[rows,length] = size(problemlist);

% Iterate over each item in the list.
for n = 1:length
    % Try the nth item in the problem set, catching any errors and outputting the reason for the error.  Problem w throws a warning that is not caught by this code, but it does not prevent the code from running.
    % We know that the results for problem w are innacurate because (C*C.')^-1 * (C*C.') should give us a 3x3 identity matrix and it does not.
    % I could not find a method of catching warnings and errors at the same time and it appears the only way to silence the warning is a user setting.
    % The second elseif term is to catch problem n if it is entered as A^2 rather than A*A.  A^2 gives an error because A is not a square matrix, while A*A warns about inner dimensions not matching .
    try
        problemlist{n}(A,B,C,D);
    catch MExc
        if MExc.identifier == 'MATLAB:dimagree'
            fprintf('%c cannot be performed, matrix dimensions do not agree \n',n+96)
        elseif MExc.identifier == 'MATLAB:innerdim'
            fprintf('%c cannot be performed, inner matrix dimensions do not agree \n',n+96)
        elseif MExc.identifier == 'MATLAB:mpower:notScalarAndSquareMatrix'
            fprintf('%c cannot be performed, inputs are not a scalar and a square matrix \n',n+96)
        else
            fprintf('An unexpected error ocurred when calculating %c',n+96)
        end
    end
end