%% Problem 1

% Part a

A = imread('bsubronco.png');

[m,n,p] = size(A);

B = reshape(A,[m,3*n]);
B = double(B)/255;

imagesc(B), axis image, colormap(gray);


% Part b
[U,S,V] = svd(B);

CR = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75];

k = floor(CR.*(3*m*n)/(m+3*n));

for p = 1:8
	 Ak = U(:,1:k(p))*S(1:k(p),1:k(p))*V(:,1:k(p))';
	 Ak = min(1,max(0,Ak));
	 Ak = reshape(Ak,[m,n,3]);
	 figure('Name',sprintf('Compression: %.f',100*CR(p)))
	 imagesc(Ak), axis image, colormap(jet), drawnow
end

% Part c
semilogy(diag(S),'-')

%% Problem 2

%% Load in the training images
load mnist_training_data;

%% Load in the test images
load mnist_test_data;

%% Determine a "good" basis for the images using the SVD
[U0,S0] = mnist_svd(images0);
[U1,S1] = mnist_svd(images1);
[U2,S2] = mnist_svd(images2);
[U3,S3] = mnist_svd(images3);
[U4,S4] = mnist_svd(images4);
[U5,S5] = mnist_svd(images5);
[U6,S6] = mnist_svd(images6);
[U7,S7] = mnist_svd(images7);
[U8,S8] = mnist_svd(images8);
[U9,S9] = mnist_svd(images9);

%% Determine what number the test images are

correct = zeros(1,40);

for k = 1:40
	for test = 1:length(testImages)
		B = double(testImages{test});
		[m,n] = size(B);
		% Flatten B into a vector b
		b = B(:);

		% Compute the residual between the "fit" of test image and the numbers in
		% the database.
		residual = zeros(10,1);
		residual(1) = norm(B(:)-U0(:,1:k)*U0(:,1:k)'*b);
		residual(2) = norm(B(:)-U1(:,1:k)*U1(:,1:k)'*b);
		residual(3) = norm(B(:)-U2(:,1:k)*U2(:,1:k)'*b);
		residual(4) = norm(B(:)-U3(:,1:k)*U3(:,1:k)'*b);
		residual(5) = norm(B(:)-U4(:,1:k)*U4(:,1:k)'*b);
		residual(6) = norm(B(:)-U5(:,1:k)*U5(:,1:k)'*b);
		residual(7) = norm(B(:)-U6(:,1:k)*U6(:,1:k)'*b);
		residual(8) = norm(B(:)-U7(:,1:k)*U7(:,1:k)'*b);
		residual(9) = norm(B(:)-U8(:,1:k)*U8(:,1:k)'*b);
		residual(10) = norm(B(:)-U9(:,1:k)*U9(:,1:k)'*b);

		[p , prediction] = min(residual);
		if prediction -1 == labels(test)
			correct(k) = correct(k) + 1;
		end
	end
end

plot(1:40,correct)
