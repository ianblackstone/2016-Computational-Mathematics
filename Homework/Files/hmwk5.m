%%  Homework #5
%   Ian Blackstone, Helena-Nikolai Fujishin
%   Math 365, Fall 2016
%%  Problems
function hmwk1()
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
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

%%  Problem #1 :Oscillatory Data Fitting NCM 5.8
function prob1()
% In the problem we were given the data:
tdata = [1:25];
ydata = [5.0291 6.5099 5.3666 4.1272 4.2948 6.1261 12.5140 10.0502 9.1614 7.5677 7.2920 10.0357 11.0708 13.4045 12.8415 11.9666 11.0765 11.7774 14.5701 17.0440 17.0398 15.9069 15.4850 15.5112 17.6572];
y1 = ydata';
% We know that the system is set up like y=Bt where B is a column vector of 
% our coeffiencts Beta(1) and Beta(2) to fit the linear equation y= Beta(1)+Beta(2)*t

% PART A
t1 = ones(length(tdata),1); %  We create a new vector to be our t for the y=Bt 
%  by creating a column of ones...
tnew =[t1 , tdata']; % ...and then setting it so that we have a column of 
%  ones next to a column made of tdata
B = tnew\y1; %  To solve for our coefficients we use the backslash operator
y2 =B(1)+(B(2).*tdata); %  And plug it into our equation setting the new equation
plot(tdata,ydata,'bo',tdata,y2,'r-'); %  We plot our data points and our given line
title('Fitting Oscillatory Data')
ylabel('y values')
xlabel('t data')
legend('Original Data','Fit Line')

% Need to plot residuals 
plot(tdata,ydata-y2,'r*')
title('Residuals')

% PART B
% We discard the outlier by hand, note changes to some variables:
tdata2=tdata;
ydata2=ydata;
ydata2(7)=[];
tdata2(7) = [];
y1 = ydata2';
t1 = ones(length(tdata2),1);
tnew =[t1 , tdata2'];
B = tnew\y1;
y2new =B(1)+(B(2).*tdata2);
% We plot the residuals again and note the differences:
plot(tdata2,ydata2,'bo',tdata2,y2new,'r-');% We plot our data points and our given line without outlier
title('Fitting Oscillatory Data Excluding Outliers')
ylabel('y values')
xlabel('t data')
legend('Original Data- Without Outlier','Fit Line')
% Need to plot residuals 
plot(tdata2,ydata2-y2new,'r*')
title('Residuals Excluding Outlier')

% PART C
% We fit the data, excluding the outlier, to
% y(t)=Beta(1)+Beta(2)t+Beta(3)sin(t)
tnow = [ones(length(tdata2),1), tdata2', sin(tdata2')]; % We make our new matrix without the outlier for t
B = tnow\y1;
y3 =B(1)+(B(2).*tnow)+(B(3).*sin(tnow));
plot(tdata2,ydata2,'bo',tnow,y3,'r-');% We plot our data points and our given line without outlier
title('Fitting Oscillatory Data Excluding Outliers')
ylabel('y values')
xlabel('t data')
legend('Original Data- Without Outlier','Fit Line')

% PART D
% We evaluate our equation above over [0,26], including the outlier marked
% with *
tt=linspace(0,26,100);
yy =B(1)+(B(2).*tt)+(B(3).*sin(tt));
% The question says to plot the third fit on a finer grid, but it doesn't
% say whether or not to include the outliner, only to mark it so that is
% what I have done:
plot(tdata(7),ydata(7),'b*',tdata2,ydata2,'o',tt,yy,'-')
legend('outlier','data excluding outlier','best fit line without outlier')
end

%%  Problem #2 : Longley Data Set NCM 5.11
function prob2()
load longley.dat
y = longley(:,1);
X = longley(:,2:7);
X1 =[ones(length(X),1),X];
% PART A 
% Use backslash operator to compute our Bs
B=X1\y;
% PART B
% Compare our Bs with the known Bs
% I will do this by calculating the errors:
KnownB= [-3482258.63459582; 15.0618722713733; -0.358191792925910e-01; -2.02022980381683;  -1.03322686717359;-0.511041056535807e-01; 1829.15146461355];
error = ((KnownB-B)./KnownB)*100;
plot(0:6,error','ro')
legend('Percent Errors of B versus Known B')
xlabel('Corresponds to B(j) where j=0:6')
ylabel('((KnownB-B)./KnownB)*100')
% PART C
% Use errorbar to plot y with error bars whose magnitude is the difference
% between y and the least squares fit
lsf=X1*B;
errorbar(y,y-lsf);
title('Error between y and least squares fit')

% PART D
% Use corrcoef to compute the correlation coefficients for X without the 
% column of 1’s. Which variables are highly correlated?
coefs = corrcoef(X);
% I dont exactly understand what this is asking for or how to do it.


% PART E
% Normalize the vector y so that its mean is zero and its standard deviation
% is one.
y = y - mean(y(:));
y = y/std(y(:));
% Do the same thing to the columns of X. Now plot all seven normalized
% variables on the same axis. Include a legend.
X = X - mean(X(:));
X = X/std(X(:));
n=linspace(-10,10,16);
plot(n,y,'r*',n,X(:,1),'bo',n,X(:,2),'go',n,X(:,3),'b+',n,X(:,4),'co',n,X(:,5),'kd',n,X(:,6),'k*')
title('Normalized Variables..Problem 5.11 NCM')
ylabel('Variables from columns of y and X')
xlabel('n=linspace(-5,5,16)')
legend('y','X(:,1)','X(:,2)','X(:,3)','X(:,4)','X(:,5)','X(:,6)')
end

%%  Problem #3 : Fitting Planetary Orbits NCM 5.12
function prob3()
x = [1.02 .95 .87 .77 .67 .56 .44 .30 .16 .01]';
y = [0.39 .32 .27 .22 .18 .15 .13 .12 .13 .15]';

% PART A
% Determine the coefficients in the quadratic form that fits these data in 
% the least squares sense by setting one of the coefficients equal to one and
% solving a 10-by-5 overdetermined system of linear equations for the other
% five coefficients. Plot the orbit with x on the x-axis and y on the y-axis.
% Superimpose the ten data points on the plot.
M=[(x.^2), (x.*y) (y.*y) x y]
t=ones(10,1)*-1;
c=M\t;
Z= c(1)*x.^2 + c(2)*x.*y + c(3)*y.^2 + c(4)*x + c(5)*y+1% Maybe you can help me figure out why meshplot and contour isn't working?
[X,Y] = meshgrid(min(x):10:max(x),min(y):10:max(y));
figure
plot(x,y,'bo')

% PART B
% (b) This least squares problem is nearly rank deficient. To see what effect this
% has on the solution, perturb the data slightly by adding to each coordinate
% of each data point a random number uniformly distributed in the interval
% [?.0005, .0005]. Compute the new coefficients resulting from the perturbed
% data. Plot the new orbit on the same plot with the old orbit. Comment on
% your comparison of the sets of coefficients and the orbits.
xnew= x+((rand)-(0.5));
ynew= y+((rand)-(0.5));
Mnew=[(xnew.^2), (xnew.*ynew) (ynew.*ynew) xnew ynew]
tnew=ones(10,1)*-1;
cnew=Mnew\tnew;
Znew= cnew(1)*xnew.^2 + cnew(2)*xnew.*ynew + cnew(3)*ynew.^2 + cnew(4)*xnew + cnew(5)*ynew+1
% I think I need to do Meshgrid again but it isn't working for the first one
% so I am just going to consult the all-powerful Ian..
[Xnew,Ynew] = meshgrid(min(xnew):10:max(xnew),min(ynew):10:max(ynew));
figure
plot(x,y,'bo',xnew,ynew,'ro')
legend('original data','pertubated data')
end


%%  Problem #4 : Simpson's Rule
function prob4()
% PART A
% Use Matlab’s symbolic tool box to obtain exact answers to I(f1) (and
% I(f2))
f1 = @(x) (-1.+x).^2.*exp(-x.^2)
f2 = @(x) 2*(1./(1.+x).^2)
I1a = integral(f1,-1,1)
I2a = integral(f2,-1,1)

% PART B
% Write Matlab Function for approximating a general integral. 
function I = simps(f,a,b,n)
%  The function implements Simpson's rule
h = (b-a)/n;
x = zeros(1,n+1);
x(1) = a;
x(n+1) = b;
p = 0;
q = 0;
%  Define the x-vector
for i = 2:n
    x(i) = a + (i-1)*h;
end
%  Define the terms to be multiplied by 4
for i = 2:((n+1)/2)
    p = p + (f(x(2*i -2)));
end
%  Define the terms to be multiplied by 2
for i = 2:((n-1)/2)
    q = q + (f(x(2*i -1)));
end
%  Calculate final output
x = (h/3)*(f(a) + 2*q + 4*p + f(b));
end

% I1b = simps(f1,-1,1,4)
I2b = simps(f1,-1,1,30)

% PART C
% Using your Simpson’s rule method from part (a), compute an approximation to I(f1) and
% I(f2) for n = 4, 8, 16, 32, 64, 128.

% Report the magnitude of the error in these approximations for each n in a nice table.

% Produce a plot of the magnitude of the error vs. 1/n (on a log-log scale).

% For the I(f1) verify that the error is decreasing like O(1/n4).

% Does this rate of decrease in the error appear to be true for I(f2)? 

% What rate does the error appear to decrease for this integral?
end