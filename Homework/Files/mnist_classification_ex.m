%% Load in the training images
load mnist_training_data;

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

%% Load in the test images
load mnist_test_data;

%% Determine what number the test images are
k = 10;     % Number of basis vectors to use
test = randperm(numel(testImages),1);   % Image to test from the test images array.
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

% Output the results in a nice table.
numbers = (0:9)';
minRes = min(residual);
likelyNumber = residual == minRes;
correctNumber = numbers == labels(test);
table(numbers,residual,likelyNumber,correctNumber)
% Get the predicted and correct numbers.
predicted = numbers(likelyNumber);
correct = labels(test);
if predicted == correct
    fprintf('The prediction is correct\n');
else
    fprintf('The prediction is incorrect\n');
end

% Display the image with the predicted and correct image
imagesc(B), colormap(gray), axis off;
title(sprintf('Predicted=%d, Actual=%d',predicted,correct));
