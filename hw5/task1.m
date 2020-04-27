pkg load statistics

function count_risks(mu, sigma, n, m)
  C1 = 0.4;
  X = sort(normrnd(mu, sigma, m, n));
  #mean(X)
  risk1 = mean(mean(X) - sigma ^ 2 / m)
  if mod(m, 2) == 1
    med = X(fix(m / 2) + 1, :);
  else
    med = (X(fix(m / 2), :) + X(fix(m / 2) + 1, :)) / 2;
  endif

  risk2 = mean(med - pi * sigma ^ 2 / (2 * m))

  risk3 = mean((X(1, :) + X(m, :)) / 2 - C1 * sigma ^ 2 / log(m))
endfunction

mu = 1;
sigma = 3;
m = 100;
count_risks(mu, sigma, 1000, m);
count_risks(mu, sigma, 10000, m);