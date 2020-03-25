pkg load statistics

function res = test(mu, sigma, n)
    X = sort(normrnd(mu, sigma, n, 1));
    dist = -1;
    for i=1:n
        X_i = X(i, 1);
        F_X_i = normcdf(X_i, mu, sigma);
        current_val = max(abs(F_X_i  - i / n), abs(F_X_i - (i - 1) / n));
        dist = max(dist, current_val);
    endfor
    res = dist * sqrt(n);
endfunction

n = 100;
mu = 1;
sigma = 1;

t = mu - 3 * sigma : 0.5 : mu + 3 * sigma;
F_x = normcdf(t, mu, sigma);

X = sort(normrnd(mu, sigma, n, 1));
F_n = 1 / n : 1 / n : 1;
[a, b] = stairs(X, F_n);

u = 1.36;
delta = u / sqrt(n); % ? * T
plot(a, b, t, F_x, a, max(b - delta, 0), a, min(b + delta, 1))

total = 0;
% kolmogorov
for i=1:10^2
    a = test(mu, sigma, 100);
    if (a >= 0.95)
        total += 1;
    endif
endfor
