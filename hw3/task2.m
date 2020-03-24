pkg load statistics

n = 100;
a = 20;
b = 80;

df_x = 0:0.01:n;
df_y = unifcdf(df_x, a, b);

edf_x = sort(unifrnd(a, b, 1, n));
edf_y = 1/n:1/n:1;

[st_a, st_b] = stairs(edf_x, edf_y);

gamma = 0.95;
u = 1.36; 
delta = u / sqrt(n);
edf_y_minus = max(0, st_b - delta);
edf_y_plus = min(1, st_b + delta);
plot(df_x, df_y, st_a, st_b, st_a, edf_y_minus, st_a, edf_y_plus);