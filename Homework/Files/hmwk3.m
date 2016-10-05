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

n=7;

r = [1:3:3*n].*[(-1).^(0:n-1)];
c = [1:2:2*n].*[(-1).^(0:n-1)];

B = toeplitz(c,r)

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

Nmax = sqrt(maxmem/8)

% part b

A = rand(int16(Nmax));
B = rand(int16(Nmax));

tic
A*B;
t = toc;

Nops = Nmax^2 * (Nmax-1) + Nmax^3;
Flops = Nops/t;
Hz = Flops/4;
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
Ttimes = 2/3 * Nvec.^3 / Flops;

figure
loglog(Nvec,lutimes,'o',Nvec,Ttimes,'-')

%% Problem 6

T0 = 293.15;
N = 200;
T = zeros(N,1);
T(1) = T0;
T(end) = T0;

b = zeros(N,1);

h = 1/N;

S = 10^4;
p = 10;
Cp = 200;
ep = 5*10^2;
B = 10^-3;
L = 1;

f = @(x,S,p,Cp,B,L,ep) S*(2*p*Cp*B)^(-1) * exp(-((x-L/2)/ep)^2)

a = 0
for num = 0:h:1
	a = a+1
	b(a) = h^2 * f(num,S,p,Cp,B,L,ep);
end

s = ones(N,1);

A = spdiags([-0.5*s 0*s -0.5*s], -1:1,N,N);

B = A\T

%% Problem 7

load BSUSurferResults;

[row,col] = size(G);

spar = 1 - nnz(G)/(row * col);

spy(G)
title(sprintf('%f',spar))

%  The following code was added to the end of pagerank.m to outpout the bottom 10 pages as well as the top 10
%[~,id] = sort(x,1,'ascend');
%fprintf('\n#     PageRank     Page\n'); 
%for j=1:10
%   fprintf('%02d    %1.2e     %s\n',n - 10 + j,x(id(j)),U{id(j)});
%end

pagerank(U,G,0.85);

pagerank(U,G,0.95);

% The top rated sites are mostly unchanged between the two weightings with only two sites swapping positions
% but the bottom ten are almost completely different, with only index.boisestate.edu/feed remaining in the
% lowest ranked pages in both weightings.

%% Problem 8

load('C:\Users\Ian\Documents\GitHub\2016-Computational-Mathematics\Homework\Files\bcsstk38.mat')
A = Problem.A;

[row1,col1] = size(A);

spar1 = 1 - nnz(A)/(row1*col1);

spy(A)
title('Sparsity of the matrix A')
legend(sprintf('Sparsity %f',spar1))

B = chol(A,'lower');

[row2,col2] = size(B);

spar2 = 1 - nnz(B)/(row2 *col2);

compspar = spar2/spar1;

spy(B)
title('Fill in of the lower Cholesky decomposition of A')
legend(sprintf('Fill in: %f',1 - compspar))

c = symamd(A);
C = A(c,c);

[row3,col3] = size(C);

spar3 = 1 - nnz(C)/(row3*col3);
compspar2 = spar3/spar1;

spy(C)
title('Fill in of the symamd of A')
legend(sprintf('Fill in: %f',1 - compspar2))