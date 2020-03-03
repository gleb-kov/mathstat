pkg load statistics;

function res = f(x)
  res = sin(x) / (x^2 + 1);
endfunction

function monte_carlo(n)
  l = 0;
  r = 5;
  y = 0.95;
  Q = norminv((y + 1) / 2);
  X = unifrnd(l, r, 1, n);
  F_x = arrayfun(@f, X) * (r - l);
  v = mean(F_x);
  delta = (std(F_x) * Q) / sqrt(n);
  printf("%g %g %g\n", v, v - delta, v + delta);
  printf("Delta is %g\n\n", delta); 
endfunction

printf("Sample answer = %g\n\n", quad(@f, 0, 5));
monte_carlo(10^4);
monte_carlo(10^6);
