%% ncm 7.18

% [x' , y' , theta' , v'] = [ v*cos(theta) , v*sin(theta) , -g/v * cos(theta) , -D/m -g*sin(theta)]

% D = cps/2 * ((x' - w(t))^2 + y'^2)

% w = 0;
% w = 10;

% w = -10 for even t, otherwise 0.
% if mod(floor(2),2) == 0
% 	w = -10;
% else
% 	w=0;

% w = 10*randn

% use event detection like in the paratrooper problem.

%% Eigenvalues and vectors
% A = [4/5 , 3/10 ; 1/5 , 7/10];
% If we multiply A by some vector we get a vector of different length and direction.
% A*[1 ; 2] = [7/5 , 8/5]

x = [1;2];
A = [4/5 , 3/10 ; 1/5 , 7/10];
A*x

quiver(x(1),x(2),'b-')
quiver(0,0,x(1),x(2),'b-')
hold on
u = A*x;
quiver(0,0,u(1),u(2),'r-')
shg

u - x

x = [3/5 ; 2/5];

u = A*x

u - x

x = [1 ; -1]

u = A*x

% These last two values of x just scale x, so these vectors are eigenvectors of A.
% A*x = c*x , where c is a constant known as an eigenvalue.  Traditionally c is lambda.

% A*x - c*x = 0
% (A-c*I)*x = 0
% det(A-c*I) = 0
% This is the characteristic equation

% A = [1 , 4 ; 2 , 3];
% A - c*I = [1-c , 4 ; 2 , 3-c]
% (1-c)(3-c) - (4)(2)
% (c^2 - 4*c - 5) = 0
% (c-5)(c+1) = 0
% c = 5 , c = -1

% We have a few options to do this for large arrays:
% Power method
% Rayleigh Quotient iteration
% Inverse iteration
% QR algorithm

% we have to use one of these methods because beyond a 6x6 matrix
% there are no more equations available to locate the roots of an equation.

% The power method:
% set x as some vector
x = rand(2,1);
% repeat this operation : x = A * x;

for k = 1:200;
	x = A*x;
end

x

% Rayleigh quotient

% If x is an eigenvector of A then
% c = x^T * A * x / (x^T * x)

x' * A * x / (x' * x)

%%  eig function gives the eigenvalues

eig(A)

% To get the vectors we can call the function like this:

[V , c] = eig(A)

% V will have the vectors stored as columns and the
% diagonal elements of c are the values.

%% Suppose A is an n by n matrix with n linearly independent
% eigenvectors X1, X2, X3,..., Xn.
% letting S = [x1|x2|x3|...|xn] we have
% A = S^-1 * c * S , where c is the eigenvalues arranged as a diagonal matrix.
% A^2 = S^-1 * c * S * S^-1 * c * S;
% since S * S^-1 is the identity matrix this can be collapsed into
% A^2 = S^-1 * c^2 * S^2
% This can be generalized for all powers:
% A^n =  S^-1 * c^n * S

%% Singular value decompisition
% let A be an m*n matrix then there exists matrices U(m*m) , V(n*n) , and Sigma(m*n) (Sigma is a diagonal matrix)
% such that
% A = U*Sigma*V^T
% where the diagonal elements of Sigma are in descending order and all greater than or equal to 0.

% Alternatively we can use the outer products:
% A = sum( U_j * Sigma_j * V_j^T )
% Since Sigma is ordered by descending elements we can approximate this by just taking the largest k values since other values will be small.
