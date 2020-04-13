pkg load statistics

clc;
clear all;

function res = test_Chi2_1(tests, n, m)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests    
    X = sort(normrnd(mu, sigma, n, 1));

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
  printf("Нормальное распределение проходит проверку гипотезы о нормальном распределении\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

function res = test_Chi2_2(tests, n, m)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests
    X = sort(normrnd(mu, sigma, n, 1));

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    x_coords = [];
    walls = [];
    for i = 1 : m
      cur_l = (i - 1) * delta + l;
      cur_r = cur_l + delta;
      walls(i, 1) = cur_l;
      walls(i, 2) = cur_r;
      x_coords(i) = (cur_r + cur_l) / 2;
    endfor

    fixed_cnt_in_bucket = [];
    fixed_walls = [];
    fixed_x_coords = [];
    fixed_m = 0;
    for i = 1 : m
      if (i == 1 || fixed_cnt_in_bucket(fixed_m) >= 6)
        fixed_m = fixed_m + 1;
        fixed_walls(fixed_m, 1) = walls(i, 1);
        fixed_cnt_in_bucket(fixed_m) = 0;
        fixed_x_coords(fixed_m) = 0;
      endif
      fixed_walls(fixed_m, 2) = walls(i, 2);
      fixed_cnt_in_bucket(fixed_m) = fixed_cnt_in_bucket(fixed_m) + cnt_in_bucket(i);
      fixed_x_coords(fixed_m) = fixed_x_coords(fixed_m) + cnt_in_bucket(i) * x_coords(i);
    endfor
    
    fixed_x_coords = fixed_x_coords ./ fixed_cnt_in_bucket;
    
    #Выборочное среднее
    E = sum(fixed_x_coords .* fixed_cnt_in_bucket) / n;
    #Выборочная дисперсия
    D = sum((fixed_x_coords - E) .^ 2 .* fixed_cnt_in_bucket) / n;
    SQRT_D = sqrt(D);

    P = [];
    for i = 1 : fixed_m
      P(i) = normcdf(fixed_walls(i, 2), E, SQRT_D) - normcdf(fixed_walls(i, 1), E, SQRT_D);
    endfor

    hi2 = sum(((fixed_cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 < chi2inv(1 - alpha, fixed_m - 1 - 2));
  endfor
  printf("Нормальное распределение проходит проверку гипотезы о нормальном распределении\n");
  printf("Данные сгруппированы, чтобы выполнялось nj >= 6\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

function res = test_Chi2_3(tests, n, m)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests
    X = sort(normrnd(mu, sigma, n, 1));

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    x_coords = [];
    walls = [];
    for i = 1 : m
      walls(i, 1) = l + delta * (i - 1);
      walls(i, 2) = l + delta * i;
      x_coords(i) = (walls(i, 1) + walls(i, 2)) / 2;
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
  printf("Нормальное распределение проходит проверку гипотезы о равномерном распределении\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

n = 10 ^ 6;
mu = 1;
sigma = 1;
m = 10 ^ 2;
tests = 10 ^ 3;

X = sort(normrnd(mu, sigma, n, 1));

l = min(X);
r = max(X);
delta = (r - l) / m;
x_coords = [];
y_coords = hist(X, m)' / (n * delta);

for i = 1 : m
  cur_l = (i - 1) * delta + l;
  cur_r = cur_l + delta;
  x_coords(i) = (cur_r + cur_l) / 2;
endfor

x_coords_for_normpdf = l:0.1:r;
stairs(x_coords, y_coords);
hold on;
plot(x_coords_for_normpdf, normpdf(x_coords_for_normpdf, mu, sigma));

printf("Размер выборки = %d\n", n);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", m);

printf("\n");

# PART 2

test_Chi2_1(tests, 10 ^ 4, m);

printf("\n");

test_Chi2_2(tests, 10 ^ 4, m);

printf("\n");

test_Chi2_3(tests, 10 ^ 4, m);

