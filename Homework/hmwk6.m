%% Problem 1
a = 0;
b = 1;
n = 100;

f = @(x,k) x^k;
exct = @(x,k) 1/(k+1) * x^(k-1);

for k = 0:n
Q = 1/2 * (f(a,k) + f(b,k));

err = abs(exct(a,k)-Q);

if err < 1E-14
	fprintf('m = %.i \n',k-1)
	break
end
end

%% problem 2

