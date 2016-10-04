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

A(find(A<0.5)) = 0

%% Problem 5

% part a

user = memory;
maxmem = user.MaxPossibleArrayBytes;

%Nmax = sqrt(maxmem/8)
Nmax = 3000;

% part b

A = rand(int16(Nmax));
B = rand(int16(Nmax));

tic
A*B;
t = toc;

Nops = Nmax^2 * (Nmax-1) + Nmax^3;
Flops = Nops/t;
Cores = feature('NumCores');
Hz = Flops/(4*Cores);
GHz = Hz / 10e9;
fprintf('This computer operates at %.2f GHz \n',GHz)

% part c

lutimes = zeros(50,1);

N = 2;
Nmaxlog = log10(Nmax);
Nvec = logspace(N,Nmaxlog);

for n = 1:numel(Nvec)
	A = rand(int16(Nvec(n)));
	b = rand(int16(Nvec(n)),1);
	tic
	A\b;
	lutimes(int16(n)) = toc;
end

% Need to find this formula
Ttimes = 3/2 * Nvec.^3 / Flops;

figure
loglog(Nvec,lutimes,'o',Nvec,Ttimes,'-')

%% Problem 6

T0 = 293.15
N = 100
T = zeros(N,1)
T(1) = T0
T(end) = T0

f = @(x,S,p,Cp,B,L,ep) S*(2*p*Cp*B)^(-1) * exp(-((x-L/2)/ep)^2)

%% Problem 7



%% Problem 8

load('C:\Users\Ian\Documents\GitHub\2016-Computational-Mathematics\Homework\Files\bcsstk38.mat')
A = Problem.A;

[row1,col1] = size(A);

spar1 = 1 - nnz(A)/(row1*col1);

figure
title('Sparsity of the matrix A')
spy(A)
legend(sprintf('Sparsity %f',spar1))

B = chol(A,'lower');

[row2,col2] = size(B);

spar2 = 1 - nnz(B)/(row2 *col2);

compspar = spar2/spar1;

figure
title('Fill in of the lower Cholesky decomposition of A')
spy(B)
legend(sprintf('Fill in: %f',1 - compspar))

c = symamd(A);
C = A(c,c);

[row3,col3] = size(C);

spar3 = 1 - nnz(C)/(row3*col3);
compspar2 = spar3/spar1;

figure
title('Fill in of the symamd of A')
spy(C)
legend(sprintf('Fill in: %f',1 -compspar2))