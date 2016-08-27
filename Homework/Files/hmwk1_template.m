%% Homework #1
%  Names of team members
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
% inner radius is 3.21 and whose outer radius is 3.56. 
% The result is included in the report.
function prob1()

%  Your work goes here

end

%% Problem #2 : Plotting a simple relationship
% In this problem, we plot the relationship between energy and magnitude of
% an earthquake according to the Richter scale.  The plot is included in
% the report.
function prob2()

%  Your work goes here

end

%% Problem #3 : Simple vectorization
% Compute the squared sum of a vector of entries.
function prob3()

%%%
% Create the arrays and then include the code for computing the sum of the
% squares of the entries here.

end


%% Problem #4 : Simple anonymous function handle; plotting
function prob4()

%%%
% Create function handles for two functions and
% construct a third composite function. 

% Anonymous function handles go here

% Construct a vector of equally spaced points over [-3,3] 

% Plot the results

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

end

%% Problem #6 : Generating sequences of numbers with a for loop
% Use matrix multiplication to generate a famous sequence from mathematics.
function prob6()

%  Your work goes here

end

%% Problem #7 : Quadratic formula
% NCM, problem 1.38. The quadratic formula is relatively simple, but to
% implement it properly on the computer in floating point arithmetic takes
% some care.
function prob7()

%  Your work goes here

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
