function p = baryinterp(x,f,u)
% BARYINTERP Polynomial interpolation using barycentric interpolation
% formula.
% 
% P = BARYINTERP(X,F,U) Evaluates the polynomial interpolant of the data
% {X,F} at the points in U.
%
% Example:
% x = -cos((0:20)*pi/20);
% f = 1./(1 + 25*x.^2);
% u = linspace(-1,1,101);
% p = baryinterp(x,f,u);
% plot(x,f,'ko',u,p,'b-')

[mu,nu] = size(u);

% Make everything a column vector
x = x(:);
u = u(:);
f = f(:);

% Get the number of entries in t
m = numel(u);

% Compute the Barycentric weights
w = barywghts(x);

% Initialize values to store the polynomial
p = zeros(size(u));

for i=1:m
    %
    % Check if the point t(i) is any of the interpolation points.
    % If it is then set p(i) equal to the function value corresponding
    % to the point the are equal.
    %
    interpPt = u(i) == x;
    if any(interpPt)
        temp = f(interpPt);
        p(i) = temp(1);
    else
        % Compute the denominator
        denominator = w./(u(i) - x);

        % Compute the numerator
        numerator = bsxfun(@times,f,denominator);

        % Sum the numerator divided by denominator
        p(i) = sum(numerator)/sum(denominator);
    end
end

p = reshape(p,mu,nu);
    
% Function to compute the barycentric weights
function w = barywghts( x )
n = length(x);
z= (eye(n) + repmat(x,1,n) - repmat(x,1,n).');
w = 1./prod(z,2);
end

end