function res = findvar(w,t,y,h)
% Make an array of zeros.
m = length(t);
n = length(w);
X = zeros(m,n);

% fill the array with the values of sin(w1t) and sin(w2t).
for j = 1:n
X(:,j) = sin(w(j)*t);
end

% Solve the linear equations to find the amplitude of each wave function.
a = X\y;
z = X*a;

% find the residual.
res = norm(z-y);

% plot the function as it evolves.
set(h(2),'ydata',z);
set(h(3),'string',sprintf('a1: %.3f  a2: %.3f  w1: %.3f   w2: %.3f',a,w))
end