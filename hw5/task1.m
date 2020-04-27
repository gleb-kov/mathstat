pkg load statistics

mu = 1;
sigma = 3;
C1 = 0.4;
n = 100; # кол-во тестов
m = 100; # кол-во элементов в выборке
X = sort(normrnd(mu, sigma, n, m));

risk1 = mean(mean(X) - sigma ^ 2 / n)
if mod(n, 2) == 1
  med = X(fix(n / 2) + 1, :);
else
  med = (X(fix(n / 2), :) + X(fix(n / 2) + 1, :)) / 2;
endif

risk2 = mean(med - pi * sigma ^ 2 / (2 * n))

risk3 = mean((X(1, :) + X(n, :)) / 2 - C1 * sigma ^ 2 / log(n))