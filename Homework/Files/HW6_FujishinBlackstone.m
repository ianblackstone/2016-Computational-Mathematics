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
    fprintf('%s : Success!\n',msg);
catch me
    fprintf('%s : Something went wrong.\n',msg);
    fprintf('%s\n',me.message);
end
fprintf('\n');
end

%% Problem #1 : Degree of Precision
function prob1()
%Part A
a=0;
b=1;
n = 100;
f = @(x,k) x^k;
exct = @(x,k) 1/(k+1)*x^(k-1);

for k = 0:n
T = (b-a)/2*(f(a,k) + f(b,k));
err = abs(exct(b,k)-T);
if err >= 1E-14
	fprintf('m for part a = %.i \n',k-1)
	break
end
end

%Part B: Simpson's rule
c =(a+b)/2;
h =(b-a);
for k = 0:n
    M = h*f((a+b)/2,k);
    T = h*(f(a,k)+f(b,k))/2;
S = ((2/3)*M)+((1/3)*T);
err = abs(exct(b,k)-S);
if err >= 1E-14
	fprintf('m for part b = %.i \n',k-1)
	break
end
end

end

%% Problem #2 :
function prob2()
%Our given/Known Values
xi = [-0.96816023950763 -0.83603110732664 -0.61337143270059  -0.32425342340381 0 0.32425342340381 0.61337143270059 0.83603110732664 0.96816023950763]';
ci = [0.081274388361574 0.18064816069486 0.26061069640294 0.31234707704 0.33023935500126 0.31234707704 0.26061069640294 0.18064816069486 0.081274388361574]';
n = 9; %Our iteration number

%Part A f1
f1= @(x) ((x-1).^2).*exp(-(x.^2))
I1 = integral(f1,-1,1); %Exact integral of f1
GQForm1=0; %Start our initial use of GQ formula
for i =1:n
GQForm1 = GQForm1+(ci(i)*f1(xi(i))); %Summation of the GQ formula
abserror1(i) = abs(I1-GQForm1);
end
iteration = [1:n]';
error = abserror1';
UsingGQf1 = table(iteration,xi,ci,error)

%Part A f2
f2 = @(x) 2*(1./(1+x.^2))
I2 = integral(f2,-1,1);
GQForm2=0;
for i =1:n
GQForm2 = GQForm2+(ci(i)*f2(xi(i))); %Summation of the GQ formula
abserror2(i) = abs(I2-GQForm2);
end
iteration = [1:n]';
error = abserror2';
UsingGQf2 = table(iteration,xi,ci,error)

%Part B
%Let us set-up the Simpson's Rule:
S1 = simps(-1,1,8,f1);
errorS1=abs(I1-S1)
S2 = simps(-1,1,8,f2);
errorS2=abs(I2-S2)
%To compare errors side by side:
errors1 = [errorS1;abserror1(end)];
errors2 = [errorS2;abserror2(end)];
Array = {'Simpsons Rule','GQ Formula'};
T=table(errors1,errors2,'RowNames',Array)
%From this we can compare and see the GQFormula is much more accurate than
%Simpsons

%Part C
%Determine degree of precision of 9-point GQ formula of f1
a = -1;
b = 1;
n = 9;
f1 = @(x,k) ((x-1).^2).*exp(-(x.^2)).^k
for k = 0:n
T = (b-a)/2*(f1(a,k) + f1(b,k));
exct = @(x) integral(@(x) f1(x,k),a,b);%Since the function alone won't integrate, I have to make bypass function
err = abs(exct(b)-T);
if err >= 1E-14
	fprintf('m for part f1 = %.i \n',k-1)
	break
end
end
%Determine degree of precision of 9-point GQ formula of f2
a = -1;
b = 1;
n = 9;
f2 = @(x,k) (2*(1./(1+x.^2))).^k;
for k = 0:n
T = (b-a)/2*(f2(a,k) + f2(b,k));
exct = @(x) integral(@(x) f2(x,k),a,b);%Since the function alone won't integrate, I have to make bypass function
err = abs(exct(b)-T);
if err >= 1E-14
	fprintf('m for part f2 = \n',k-1)
	break
end
end

end
    
    
%% Problem #3 : 
function prob3()

end


%% Problem #4 :
function prob4()

end
%% Problem #5 : 
function prob5()

end
