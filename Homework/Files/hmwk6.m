%% Problem 1

% part a
a = 0;
b = 1;
n = 10;

f = @(x,k) x^k;
exct = @(x,k) 1/(k+1) * x^(k+1);

for k = 0:n
Q = (b-a)/2 * (f(a,k) + f(b,k));

err = abs(exct(b,k)-Q)

if err >= 1e-14
	fprintf('m = %.i \n',k-1)
	break
end

end

% part b
a = 0;
b = 1;
n = 10;

f = @(x,k) x^k;
exct = @(x,k) 1/(k+1) * x^(k+1);

for k = 0:n
Q = (b-a)/6 * (f(a,k) + 4*f((b+a)/2,k) + f(b,k));

err = abs(exct(b,k)-Q)

if err >= 1e-14
	fprintf('m = %.i \n',k-1)
	break
end

end

%% Problem 2
