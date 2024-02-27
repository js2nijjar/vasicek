close all;
clear all;

% nsims - number of simultations
% r_init - initial interest rate
% c - [kappa theta sigma]
% N - number of grid points
% dt - time step

nsims = 1000;
r_init = 0.0525;
tbounds = [0 1];
c = [0.1 0.02 0.01];
N = 1000;
dt = (tbounds(2)-tbounds(1))/(N-1);

% vector of times
tvec = linspace(tbounds(1), tbounds(2), N);

% pre-allocate output
rvec = zeros(nsims,N);

% set initial interest rate
rvec(:,1) = r_init;

% Weiner increments
weiner_incr = sqrt(dt)*randn(1,N-1,nsims);

% Euler-Maruyama method
for i = 1:nsims
    for j=2:numel(tvec)
        t = tbounds(1)+(j-1).*dt;
        r = rvec(i,j-1);
        mu = c(1).*(c(2)-r);
        sig = c(3);
        dW = weiner_incr(:,:,i);
        rvec(i,j) = r + mu.*dt + sig.*dW(j-1);
    end
end

figure;
hold on;

% plot all simulations
for i = 1:nsims
    plot(tvec, rvec(i,:));
end

xlabel('Time');
ylabel('Interest Rate');
title('Interest Rate Paths');
grid on;