function sinfit

% Load the data.
load('PeriodicData.mat')

% Draw a figure.
clf
shg
set(gcf,'doublebuffer','on')
h = plot(t,y,'o',t,0*t,'-');
h(3) = title('');
ylabel('y')
xlabel('t')
legend('y','fit')

% Declare initial guess for omega.
w0 = [3 7]';

% Find new values for omega and amplitude by minimizing the residual residual.
w = fminsearch(@findvar,w0,[],t,y,h);

% Draw the final function.
set(h(2),'color','black')
end