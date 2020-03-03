pkg load statistics;

function res = g(x)
  res = sqrt(1 + x * x);
endfunction

function res = g1(x)
  res = g(x) * sqrt(2 * pi * 2);
endfunction

function res = f(x)
  res = g(x) * exp(-(x + 2) ^ 2 / 4);
endfunction

function calc_value(n)
  mu = -2;
  sigma = sqrt(2);
  y = 0.95;
  T = norminv((y + 1) / 2);
  X = normrnd(mu, sigma, 1, n);
  F_x = arrayfun(@g1, X);
  v = mean(F_x);
  delta = (std(F_x) * T) / sqrt(n);
  printf("%g %g %g\n", v, v - delta, v + delta);
  printf("Delta is %g\n\n", delta); 
endfunction

printf("Sample answer = %g\n\n", quad(@f, -inf, inf));
calc_value(10^4);
calc_value(10^6);
