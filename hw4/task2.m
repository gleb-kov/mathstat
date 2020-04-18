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
    
    #Выборочное среднее
    E = mean(X);
    #Выборочная отклонение
    SQRT_D = std(X);
    
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

function res = test_Chi2_2(tests, n, m)
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
    
    #Выборочное среднее
    E = mean(X);
    #Выборочная отклонение
    SQRT_D = std(X);

    P = [];
    for i = 1 : m
      P(i) = normcdf(l + delta * i, E, SQRT_D) - normcdf(l + delta * (i - 1), E, SQRT_D);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 < chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Возьмём данные из равномерного распределения и вероятности из нормального\n");
  printf("Тогда для alpha = %d, вероятность ошибки второго рода получается - %d\n", alpha, res / tests);
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

# Проверка 1
test_Chi2_1(10 ^ 3, 10 ^ 4, m);

printf("\n");

# Проверка 2
test_Chi2_2(10 ^ 2, 10 ^ 4, m);
