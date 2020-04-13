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
    cnt_in_bucket = zeros(m, 1);
    walls = zeros(m, 2);
    for i = 1 : m
      cur_l = (i - 1) * delta + l;
      cur_r = cur_l + delta;
      walls(i, 1) = cur_l;
      walls(i, 2) = cur_r;
      cnt_in_bucket(i) = sum(X <= cur_r) - sum(X < cur_l);
    endfor

    x_coords = zeros(m, 1);
    for i = 1 : m
      x_coords(i) = (walls(i, 2) + walls(i, 1)) / 2;
    endfor
    
    #Выборочное среднее
    E = sum(x_coords .* cnt_in_bucket) / n;
    #Выборочная дисперсия
    D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
    SQRT_D = sqrt(D);

    P = zeros(m, 1);

    for i = 1 : m
      cur_l = walls(i, 1);
      cur_r = walls(i, 2);
      P(i) = normcdf(cur_r, E, SQRT_D) - normcdf(cur_l, E, SQRT_D);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    # TODO: исправить подсчет значения chi2inv
    res = res + (hi2 < chi2inv(1 - alpha, m - 1));
  endfor
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
    cnt_in_bucket = zeros(m, 1);
    walls = zeros(m, 2);
    for i = 1 : m
      cur_l = (i - 1) * delta + l;
      cur_r = cur_l + delta;
      walls(i, 1) = cur_l;
      walls(i, 2) = cur_r;
      cnt_in_bucket(i) = sum(X <= cur_r) - sum(X < cur_l);
    endfor

    fixed_cnt_in_bucket = zeros(m, 1);
    fixed_walls = zeros(m, 2);
    cur = 0;
    for i = 1 : m
      if (i == 1 || fixed_cnt_in_bucket(cur) >= 6)
        cur = cur + 1;
        fixed_walls(cur, 1) = walls(i, 1);
      endif
      fixed_walls(cur, 2) = walls(i, 2);
      fixed_cnt_in_bucket(cur) = fixed_cnt_in_bucket(cur) + cnt_in_bucket(i);
    endfor
    m = cur;
    walls = zeros(m, 1);
    cnt_in_bucket = zeros(m, 1);
    for i = 1 : m
      walls(i, 1) = fixed_walls(i, 1);
      walls(i, 2) = fixed_walls(i, 2);
      cnt_in_bucket(i) = fixed_cnt_in_bucket(i);
    endfor
    
    x_coords = zeros(m, 1);
    for i = 1 : m
      x_coords(i) = (walls(i, 2) + walls(i, 1)) / 2;
    endfor
    
    #Выборочное среднее
    E = sum(x_coords .* cnt_in_bucket) / n;
    #Выборочная дисперсия
    D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
    SQRT_D = sqrt(D);

    P = zeros(m, 1);

    for i = 1 : m
      cur_l = walls(i, 1);
      cur_r = walls(i, 2);
      P(i) = normcdf(cur_r, E, SQRT_D) - normcdf(cur_l, E, SQRT_D);
    endfor

    hi2 = sum(((cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    # TODO: исправить подсчет значения chi2inv
    res = res + (hi2 < chi2inv(1 - alpha, m - 1));
  endfor
  printf("Данные сгруппированы, чтобы выполнялось nj >= 6\n");
  printf("Для alpha = %d, успешно %d из %d\n", alpha, res, tests);
endfunction

n = 10 ^ 6;
mu = 1;
sigma = 1;
m = 100;

X = sort(normrnd(mu, sigma, n, 1));

l = min(X);
r = max(X);
delta = (r - l) / m;
x_coords = zeros(m, 1);
cnt_in_bucket = zeros(m, 1);
y_coords = zeros(m, 1);

for i = 1 : m
  cur_l = (i - 1) * delta + l;
  cur_r = cur_l + delta;
  x_coords(i) = (cur_r + cur_l) / 2;
  cnt_in_bucket(i) = sum(X <= cur_r) - sum(X < cur_l);
  y_coords(i) = cnt_in_bucket(i) / (n * delta);
endfor

# TODO: научиться строить график
x_coords_for_normpdf = l:0.1:r;
stairs(x_coords, y_coords);
hold on;
plot(x_coords_for_normpdf, normpdf(x_coords_for_normpdf, mu, sigma));

printf("Размер выборки = %d\n", n);
printf("Границы = [%d; %d]\n", l, r);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", m);

printf("\n");

# PART 2

test_Chi2_1(10 ^ 3, 10 ^ 4, m);

printf("\n");

# nj >= 6
test_Chi2_2(10 ^ 3, 10 ^ 4, m);

