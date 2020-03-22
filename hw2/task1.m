pkg load statistics

mu = -1;
sigma = 0.5;
t0 = -0.25;
F = normcdf(t0, mu, sigma);
y = 0.95;
T = norminv((y + 1) / 2);
n = 10^4;
m = 100;
X = normrnd(mu, sigma, n, m);
p = mean(X < t0);
delta = T * sqrt(f .* (1 - f) / n);
plot(1 : 1 : m, p - delta, ".-", 1 : 1 : m, p + delta, ".-", 1 : 0.1 : m, F)
