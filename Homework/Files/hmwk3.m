%% Problem 1

A = diag(1:5);

D = diag(1./diag(A))

A = diag(2.^(1:5));

D = diag(1./diag(A))

%% Problem 2

% part a
A = reshape(1:25,[5,5]);

A1 = A - triu(A)

A = repmat(1:5,[5,1]);

A2 = A - tril(A,1)

A3 = A1 + A2

% part b

B1 = repmat(2:2:10,[5,1])

B2 = reshape(25:-1:1,[5,5])

B3 = B2 - diag(diag(B2))

B4 = B3 - eye(5)

% part c

C1 = 

% part d
d = [1:9].^2

D = reshape(d,[3,3])

D1 = repmat(D,[2,2])

%% Problem 3

% part a
A = toeplitz([-2,1,0,0,0,0,0,0,0,0])

B = toeplitz([-2,1,0,0,0,0,0,0,0,1])

% part b

A = toeplitz([1,0,0,0,0,0,0,0,0],[1:8])

r = [1:6:19;-4:-6:-25];
c = [1:4:13;-3:-4:-17];
% I don't think this is how he wants it done, but it could still be written in the way he wants I think.
B = toeplitz(c(1:end-1),r(1:end-1))

%% Problem 4

A = rand(9);

% part a

A1 = max(A,[],2)

% part b

A2 = min(A,[],1)

% part c

A3 = max(max(A))

% part d

A4 = 