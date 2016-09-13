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



%% Problem 6

A = rand(2,3);
B = rand(3);
C = rand(3,2);
D = rand(3);

try
    a = 2*A + C;
    c = 3*B - 2*D;
    d = A*D;
    e = C*A;
    f = A*C;
    g = B*D;
    h = D*B;
    i = B*C;
    j = C*B;
    k = D*A*B;
    l = 2*D^2 + B;
    m = D*D;
    n = A^2;
    o = C.' * D;
    p = B*A.';
    q = -2*A.' + 5*C;
    r = B.' + A*D;
    s = 0.5 * (B + B.');
    t = 0.5 * (B-B.');
    u = C*C.';
    v = C.' * C;
    w = (C*C.')^(-1);
    x = (C.' *C)^(-1);
    y = B* (A*D).';
    z = A*D*B.';
catch error
    a = NaN;
end

try
    b = C - 3*B;
catch error
    b = NaN;
end

