pkg load statistics

clc;
clear all;

function res = test_Chi2_1(tests, n, m)
  res = 0;
  a = 20;
  b = 80;
  alpha = 0.05;
  for t = 1 : tests    
    X = sort(unifrnd(a, b, n, 1));

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    walls = [];
    x_coords = [];
    for i = 1 : m
      walls(i, 1) = l + (i - 1) * delta;
      walls(i, 2) = l + i * delta;
      x_coords(i) = (walls(i, 2) + walls(i, 1)) / 2;
    endfor
    
    #Выборочное среднее
    E = sum(x_coords .* cnt_in_bucket) / n;
    #Выборочная дисперсия
    D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
    SQRT_D = sqrt(D);

    P = [];
    for i = 1 : m
      P(i) = unifcdf(walls(i, 2), l, r) - unifcdf(walls(i, 1), l, r);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 < chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Равномерное распределение проходит проверку гипотезы о равномерном распределении\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

function res = test_Chi2_2(tests, n, m)
  res = 0;
  a = 20;
  b = 80;
  alpha = 0.05;
  for t = 1 : tests    
    X = sort(unifrnd(a, b, n, 1));

    l = min(X);
    r = max(X);

    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    walls = [];
    x_coords = [];
    for i = 1 : m
      walls(i, 1) = l + delta * (i - 1);
      walls(i, 2) = l + delta * i;
      x_coords(i) = (walls(i, 2) + walls(i, 1)) / 2;
    endfor
    
    #Выборочное среднее
    E = sum(x_coords .* cnt_in_bucket) / n;
    #Выборочная дисперсия
    D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
    SQRT_D = sqrt(D);

    P = [];
    for i = 1 : m
      P(i) = normcdf(walls(i, 2), E, SQRT_D) - normcdf(walls(i, 1), E, SQRT_D);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 < chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Равномерное распределение проходит проверку гипотезы о нормальном распределении\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

n = 10 ^ 5;
m = 10 ^ 2;
a = 20;
b = 80;


# PART 1

X = sort(unifrnd(a, b, 1, n));

l = min(X);
r = max(X);
delta = (r - l) / m;
y_coords = hist(X, m) / (n * delta);

x_coords = [];
for i = 1 : m
  cur_l = (i - 1) * delta + l;
  cur_r = cur_l + delta;
  x_coords(i) = (cur_r + cur_l) / 2;
endfor

real_y = 1 / (b - a);
stairs(x_coords, y_coords);
hold on;
plot([a b], [real_y real_y], "linewidth", 1);

printf("Размер выборки = %d\n", n);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", m);

printf("\n");

# PART 2

test_Chi2_1(10 ^ 3, 10 ^ 4, m);

printf("\n");

test_Chi2_2(10 ^ 3, 10 ^ 4, m);