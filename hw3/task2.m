pkg load statistics

function p = test(n, a, b)
  m = 100;
  X = unifrnd(a, b, m, n);
  X = sort(X);
  X_1 = X(1, :);
  F_X_1 = unifcdf(X_1, a, b);
  res = max(abs(F_X_1  - 1 / n), F_X_1);
  for i=2:m
    #printf("%d %d\n", edf_x(1, i), unifcdf(edf_x(1, i), a, b));
    X_i = X(i, :);
    F_X_i = unifcdf(X_i, a, b);
    current_val = max(abs(F_X_i  - i / m), abs(F_X_i - (i - 1) / m));
    res = max(res, current_val);
  endfor
  gamma = 0.95;
  u_gamma = 1.36;
  p = mean((sqrt(m) * res) > u_gamma)
endfunction 

n = 10;
a = 20;
b = 80;

df_x = 0:0.01:n;
df_y = unifcdf(df_x, a, b);

edf_x = sort(unifrnd(a, b, 1, n));
edf_y = 1/n:1/n:1;

test(10000, a, b);
test(1000000, a, b);

[st_a, st_b] = stairs(edf_x, edf_y);

gamma = 0.95;
u = 1.36; 
delta = u / sqrt(n);
edf_y_minus = max(0, st_b - delta);
edf_y_plus = min(1, st_b + delta);
plot(df_x, df_y, st_a, st_b, st_a, edf_y_minus, st_a, edf_y_plus);