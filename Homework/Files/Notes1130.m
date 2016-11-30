% Truncated SVD to k terms (k<n)

% A_k = sum from j = 1 to k of sigma_j * u_j * v_j^T

% Property: Of all rank k matrices, the one that minimizes no ||A-B_k||_2 is the truncated svd to k terms.

% SVD as a copression method
% A is n by n.
% Suppose we approximate A by A_k
% A_k = sum _j=1 ^k of sigma_j * u_j * V_j^T
% How much information is required to store A_k?
% Each u_j has m entries.
% Each sigma_j * V_j has n entries.
% Total: km*kn = k(m+n)
% Total for non approximation: mn
% This is a useful compression is k < mn/(m+n)

% Compression ratio: (m+n)*k / mn

% Reading an image

A = imread('bsubronco.bmp');
size(A)

A = double(A);
imagesc(A), axis image, colormap(gray);

[U,S,V] = svd(A);

% S (sigma) is the singular values of the matri.  We can see the size of these values by plotting
semilogy(diag(S),'x-')

k = 1;
Ak = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
imagesc(Ak), axis image, colormap(gray)

% This is a bad compression because we lost so much of the data.  We can try higher terms by increasing k.

% to lose the gray background we can convert back to white by scaling between 0 and 255.

Ak = round(max(min(Ak,255),0));

[m,n] = size(Ak);

for k = 1:m
	 Ak = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
	 Ak = round(max(min(Ak,255),0));
	 imagesc(Ak), axis image, colormap(gray), drawnow

end