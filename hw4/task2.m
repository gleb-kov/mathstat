pkg load statistics

clc;
clear all;

function res = test_Chi2_1(tests, n, m)
  res = 0;
  a = 20;
  b = 80;
  alpha = 0.05;
  for t = 1 : tests    
    X = unifrnd(a, b, n, 1);

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    
    P = [];
    for i = 1 : m
      P(i) = unifcdf(l + i * delta, l, r) - unifcdf(l + (i - 1) * delta, l, r);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 >= chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Равномерное распределение проходит проверку гипотезы о равномерном распределении\n");
  printf("Для alpha = %d, вероятность ошибки первого рода получается %d\n", alpha, res / tests);
endfunction

function res = test_Chi2_2(tests, n, m, d)
  res = 0;
  a = 20;
  b = 80;
  alpha = 0.05;
  for t = 1 : tests    
    X = unifrnd(a, b, n, 1);

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    
    P = [];
    for i = 1 : m
      P(i) = unifcdf(l + i * delta, l - d, r + d) - unifcdf(l + (i - 1) * delta, l - d, r + d);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 >= chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Равномерное распределение проходит проверку гипотезы о равномерном распределении\n");
  printf("Левая и правая граница изменены на %d\n", d)
  printf("Для alpha = %d, вероятность ошибки второго рода получается %d\n", alpha, res / tests)
endfunction

n = 10 ^ 6;
m = 10 ^ 2;
a = 20;
b = 80;

# PART 1

X = unifrnd(a, b, n, 1);

l = min(X);
r = max(X);
delta = (r - l) / m;
[y_coords x_coords] = hist(X, m);

real_y = 1 / (b - a);
bar(x_coords, y_coords / (n * delta));
hold on;
plot([a b], [real_y real_y], "linewidth", 1);

printf("Размер выборки = %d\n", n);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", m);

printf("\n");

# PART 2

test_Chi2_1(10 ^ 3, 10 ^ 4, m);

printf("\n");

test_Chi2_2(10 ^ 3, 10 ^ 4, m, 0.25);

printf("\n");

test_Chi2_2(10 ^ 3, 10 ^ 4, m, 0.75);

printf("\n");

test_Chi2_2(10 ^ 3, 10 ^ 4, m, 1.5);

printf("\n");
