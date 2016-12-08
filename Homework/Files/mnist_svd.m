% Computes the SVD for the given MNIST training images.
function [U,S,V] = mnist_svd(mnist_image)

% Figure out how many images there are
numImages = numel(mnist_image);
% Figure out the size of each image
[m,n] = size(mnist_image{1});

% Loop through the images storing them in a m*n-by-numImages matrix
A = zeros(m*n,numImages);
for k = 1:numImages
    A(:,k) = double(reshape(mnist_image{k},m*n,1));
end

% Compute the SVD
[U,S,V] = svd(A);

end