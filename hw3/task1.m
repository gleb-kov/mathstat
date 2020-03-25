pkg load statistics

n = 100;
mu = 1;
sigma = 1;

t = (mu - 3 * sigma) : 0.5 : mu + 3 * sigma;
F_x = normcdf(t, mu, sigma);

X = sort(normrnd(mu, sigma, n, 1));
F_n = 1 / n : 1 / n : 1;
[a, b] = stairs(X, F_n);

u = 1.36;
delta = u / sqrt(n); % ? * T
plot(a, b, t, F_x, a, max(b - delta, 0), a, min(b + delta, 1))
