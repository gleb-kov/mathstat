pkg load statistics

n = 100;
mu = 1;
sigma = 1;
X = sort(normrnd(mu, sigma, n, 1));

% F_x by normcdf

F_n = 1 / n : 1 / n : 1;
[a, b] = stairs(X, F_n);

u = 1.36;
delta = T * u / sqrt(n)
% plot
