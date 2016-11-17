%% Problem 1

% part a
a = 0;
b = 1;
n = 10;

f = @(x,k) x^k;
exct = @(x,k) 1/(k+1) * x^(k+1);

for k = 0:n
Q = (b-a)/2 * (f(a,k) + f(b,k));

err = abs(exct(b,k)-Q);

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

err = abs(exct(b,k)-Q);

if err >= 1e-14
	fprintf('m = %.i \n',k-1)
	break
end

end

%% Problem 2



%% Problem 3
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

uj = zeros(N,1);
uj(1) = a;
uj(end) = b;

f = @(x) (S/(rho*Cp*B))*(exp(-(((x-(1/2))/(eps))^2)));

u = @(x) (a-quad(@f,0,x))*(1-x) - (b+quad(@f,x,1))*x;

for x = j.*h
	uj(x) = u(x)


%% Problem 4

