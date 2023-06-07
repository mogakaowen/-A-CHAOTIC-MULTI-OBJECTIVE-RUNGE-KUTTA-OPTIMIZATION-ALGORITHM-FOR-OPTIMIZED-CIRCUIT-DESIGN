% Parameters
r = 3.7;          % Control parameter
N = 1000;         % Number of iterations
x0 = 0.5;         % Initial condition

% Initialize arrays
x = zeros(N+1,1);
x(1) = x0;

% Iterate the logistic map
for n = 1:N
    x(n+1) = r*x(n)*(1-x(n));
end

% Plot the results
plot(1:N+1, x);
xlabel('Iteration number');
ylabel('x');
title(['Logistic map (r = ' num2str(r) ', x_0 = ' num2str(x0) ')']);