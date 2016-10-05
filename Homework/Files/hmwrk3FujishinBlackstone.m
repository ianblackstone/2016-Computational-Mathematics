%% Homework #3
%  Ian Blackstone, Helena Fujishin
%  Math 365, Fall 2016

function hmwk1()

% You will then call all of your homework problem "functions" like this. 
hmwk_problem(@prob1,'prob1');
hmwk_problem(@prob2,'prob2');
hmwk_problem(@prob3,'prob3');
hmwk_problem(@prob4,'prob4');
hmwk_problem(@prob5,'prob5');
hmwk_problem(@prob6,'prob6');
hmwk_problem(@prob7,'prob7');
hmwk_problem(@prob8,'prob8');

end

function hmwk_problem(prob,msg)
%%%
% This function should be included in every assignment
try
    prob()
    fprintf('%s : Success!\n',msg);
catch me
    fprintf('%s : Something went wrong.\n',msg);
    fprintf('%s\n',me.message);
end
fprintf('\n');
end

%% Problem #1 : Diagonal matrices and their inverse
function prob1()

% Here we generate a matrix using the diag function then generate a new matrix based off
% the diagonal line of our first matrix.
A = diag(1:5);

D = diag(1./diag(A))

A = diag(2.^(1:5));

D = diag(1./diag(A))

end

%% Problem #2 : Manipulating matrices
function prob2()

% part a
% We generate a matrix using the reshape function then we subtract the upper half of
% that matrix from itself to form a new matrix.
A = reshape(1:25,[5,5]);

A1 = A - triu(A)

% Generate a matric using a repeated pattern then subtract the superdiagonal and
% below of the matrix from the original matrix.
A = repmat(1:5,[5,1]);

A2 = A - tril(A,1)

% Combine two matrices to form a new one.
A3 = A1 + A2

% part b

% Generate a matrix using a pattern.
B1 = repmat(2:2:10,[5,1])

% Generate a matrix from a list of reserved numbers.
B2 = reshape(25:-1:1,[5,5])

% take the diagonal line of our matrix B2 and subtract it from itself to make a new matrix.
B3 = B2 - diag(diag(B2))

% subtract a 5x5 identity matrix from B3.
B4 = B3 - eye(5)

% part c



% part d

% Create a list of the squares of 1 through 9.
d = [1:9].^2

% turn our list into a matrix.
D = reshape(d,[3,3])

% Repeat our matrix to make a new larger matrix.
D1 = repmat(D,[2,2])

end

%% Problem #3 : Toeplitz matrices
function prob3()

% part a
% Create two toeplitz matrices from vectors.
A = toeplitz([-2,1,0,0,0,0,0,0,0,0])

B = toeplitz([-2,1,0,0,0,0,0,0,0,1])

% part b

% create a Toeplitz matrix from two vectors, verifying that the column command takes precedence.
A = toeplitz([1,0,0,0,0,0,0,0,0],[1:8])

% Determine the size of our matrix.
n=7;

% Create two vectors, one for the row and anoither for the column.
r = [1:3:3*n].*[(-1).^(0:n-1)];
c = [1:2:2*n].*[(-1).^(0:n-1)];

% Generate a Toeplitz matrix using our vectors.
B = toeplitz(c,r)

end


%% Problem #4 : Max, min, and logical indexing
function prob4()

% Generate a random 9x9 matrix.
A = rand(9);

% part a

% Find the highest value in each row of our random matrix.
A1 = max(A,[],2)

% part b

% Find the minimum value in each row of our random matrix.
A2 = min(A,[],1)

% part c

% Find the maximum value in all of the random matrix.
A3 = max(max(A))

% part d

% Set all values in the random matrix that are less than 0.5 to 0.
A(find(A<0.5)) = 0

end

%% Problem #5 : Loading data from a file, simple statistics, and using fprintf
function prob5()

% Determine the maximum possible array that can fit in available memory.
% This command will not function on Mac OS computers.
user = memory;
maxmem = user.MaxPossibleArrayBytes;

% A matrix has N^2 components, and each component uses 8 bytes of memory.
% To find the maximum possible array size we divide the maximum array size by 8 and take the square root.
Nmax = sqrt(maxmem/8)

% If using Mac OS comment out the above lines and remove the comment from the line below.
% Nmax = 3000;

% Create two random matrices, Nmax must be converted to a 16 bit integer for large memory systems.
% Larger memory system with more than 68 GB of memory available may need to convert to int32.
A = rand(int16(Nmax));
B = rand(int16(Nmax));

% Calculate the time taken to multiply the two random matrices.
tic
A*B;
t = toc;

% Calculate the number of operations performed in A*B and use the time it took to calculate that operation
% to find the floating point operations per second and then estimate the operating speed of the processor.
Nops = Nmax^2 * (Nmax-1) + Nmax^3;
Flops = Nops/t;
Hz = Flops/4;
GHz = Hz / 10e9;
fprintf('This computer operates at %.2f GHz \n',GHz)

% part c

% Initialize a list.
lutimes = zeros(50,1);

% create a list of logarithmically spaced points.
N = 2;
Nmaxlog = log10(Nmax);
Nvec = logspace(N,Nmaxlog);

% For each item in our list of points generate random matrices and a random vecotr then time how long it takes to
% solve each matrix/vector combination.
for n = 1:numel(Nvec)
	A = rand(int16(Nvec(n)));
	b = rand(int16(Nvec(n)),1);
	tic
	A\b;
	lutimes(int16(n)) = toc;
end

% Generate a set of theoretical data points to be plotted based on the calculated speed of the processor.
Ttimes = 2/3 * Nvec.^3 / Flops;

% Plot both the theoretical line and the actual times over the matrix size on a log-log plot.
figure
loglog(Nvec,lutimes,'o',Nvec,Ttimes,'-')
title('Time taken to solve an NxN matrix by matrix size.')
xlabel('Matrix size')
ylabel('Time taken')
legend('A\\b (experimental)','A\\b (theoretical)')

% part e
% The fastest computer at this time is the Sunway TaihuLight, with more than 10 million cores and an operating
% capacity of 93 petaflops.  My computer operating at 91 gigaflops would have been the fastest computer in the
% world until late 1993 when it would have been overtaken by the Numerical Wind Tunnel.  It is likely that we
% will see exascale machines within 10 years given the pace that China has been on to reach this milestone.

end

%% Problem #6 : Heat Transfer
function prob6()

%  Your work goes here

end

%% Problem #7 : PageRank
function prob7()

% Load the results of BSUSurfer.m.
load BSUSurferResults;

% find the size of the matrix G.
[row,col] = size(G);

% Determine the sparsity of matrix G.
spar = 1 - nnz(G)/(row * col);

% Generate the connectivity matrix with title.
figure
spy(G)
title(sprintf('%f',spar))

%  The following code was added to the end of pagerank.m to outpout the bottom 10 pages as well as the top 10
%[~,id] = sort(x,1,'ascend');
%fprintf('\n#     PageRank     Page\n'); 
%for j=1:10
%   fprintf('%02d    %1.2e     %s\n',n - 10 + j,x(id(j)),U{id(j)});
%end

% run pagerank twice, with two different weightings.
pagerank(U,G,0.85);

pagerank(U,G,0.95);

% The top rated sites are mostly unchanged between the two weightings with only two sites swapping positions
% but the bottom ten are almost completely different, with only index.boisestate.edu/feed remaining in the
% lowest ranked pages in both weightings.

end

%% Problem #8 : Reducing fill-in
function prob8()

% Load the data from file, storing it as matrix A.
load('bcsstk38.mat')
A = Problem.A;

% Determine the size of matrix A.
[row,col] = size(A);


% Calculate the sparsity of A.
spar1 = 1 - nnz(A)/(row*col);

% Create the connectivity matrix.
figure
spy(A)
title('Sparsity of the matrix A')
legend(sprintf('Sparsity %f',spar1))

% Create the Cholesky decomposition of matrix A.
B = chol(A,'lower');

% Calculate the sparsity of B.
spar2 = 1 - nnz(B)/(row*col);

% Compare the sparsity of A and B.
compspar = spar2/spar1;

% Plot the connetivity matrix of B and display the fill in.
figure
spy(B)
title('Fill in of the lower Cholesky decomposition of A')
legend(sprintf('Fill in: %f',1 - compspar))

% Use the symand function and generate the Cholesky decomposition of the rearranged A to reduce fill in.
c = symamd(A);
C = chol(A(c,c),'lower');

% Calculate the sparsity of C and compare it to the sparsity of A.
spar3 = 1 - nnz(C)/(row*col);
compspar2 = spar3/spar1;

% Plot the connectivity matrix of C and display the fill in.
figure
spy(C)
title('Fill in of the symamd of A')
legend(sprintf('Fill in: %f',1 - compspar2))

end
