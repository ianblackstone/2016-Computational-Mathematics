%% Problem 1
n = input('input the number of terms to sum: ');

X = 0;
pie = zeros(1,n);

for k = 0:n
	X = X + (-1)^k / (2 * k + 1) * ((4 * (1 / 5)^(2 * k + 1) - (1/239)^(2 * k + 1)));
    pie(k+1) = X;
end

pie = 4*pie;

pieerror = abs(pie/pi - 1)*100;

semilogy(pieerror)

%% Problem 2
n = 50;

X = zeros(1,n);
X(end) = 1/n;

for k = 1:n
    X(n-k) = n-k+1 + 1/X(n-k+1);
end

%% Problem 3
