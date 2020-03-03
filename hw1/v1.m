pkg load statistics;

function monte_carlo(n)
  k = 10;
  y = 0.95;
  c = 2.21;
  a = 3; 
  Q = norminv((y + 1) / 2);
  X = rand(k, n);
  F_x = sum(X.^a);
  v = mean(F_x <= c);
  delta = Q * sqrt(v * (1 - v) / n);
  printf("%g %g %g\n\n", v, v - delta, v + delta);
  printf("Delta is %g\n", delta);
endfunction

monte_carlo(10^4);
monte_carlo(10^6);
