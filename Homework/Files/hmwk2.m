%% Problem 1
n = inut('input the number of terms to sum: '),

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

problems = {@(A,B,C,D) 2*A + C, @(A,B,C,D) C - 3*B, @(A,B,C,D) 3*B - 2*D,@(A,B,C,D) A*D, @(A,B,C,D) C*A, @(A,B,C,D) A*C, @(A,B,C,D) B*D, @(A,B,C,D) D*B, @(A,B,C,D) B*C, @(A,B,C,D) C*B, @(A,B,C,D) D*A*B, @(A,B,C,D) 2*D.' + B, @(A,B,C,D) D*D, @(A,B,C,D) A*A, @(A,B,C,D) C.' * D, @(A,B,C,D) B*A.', @(A,B,C,D) -2*A.' + 5*C, @(A,B,C,D) B.' + A*D, @(A,B,C,D) 0.5 * (B + B.'), @(A,B,C,D) 0.5 * (B-B.'), @(A,B,C,D) C*C.',@(A,B,C,D) C.' * C, @(A,B,C,D) (C*C.')^(-1), @(A,B,C,D) (C.' *C)^(-1), @(A,B,C,D) B* (A*D).', @(A,B,C,D) A*D*B.'};

A = rand(2,3);
B = rand(3);
C = rand(3,2);
D = rand(3);

[size,length] = size(problems);

for n = 1:length
    try
        problems{n}(A,B,C,D);
        fprintf('%c \n',n+96)
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