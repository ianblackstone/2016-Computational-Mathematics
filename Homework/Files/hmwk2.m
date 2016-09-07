n = input('input the number of terms to sum')

X = zeros(1,n);

for k = 0:n
	X(k + 1) = (-1)^k / (2 * k + 1) * ((4 * (1 / 5)^(2 * k + 1) - (1/239)^(2 * k + 1)));
end