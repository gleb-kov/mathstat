pkg load statistics

clc;
clear all;

function res = test_Chi2_1(tests, n, m)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests    
    X = normrnd(mu, sigma, n, 1);

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
    res = res + (hi2 >= chi2inv(1 - alpha, m - 1 - 2));
  endfor
  printf("Normal distribution satisfies the hypothesis about normal distribution\n")
  printf("For alpha = %d, probability of type I error is %d\n", alpha, res / tests)
endfunction

function res = test_Chi2_2(tests, n, m)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests
    X = normrnd(mu, sigma, n, 1);

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    walls = [];
    for i = 1 : m
      cur_l = (i - 1) * delta + l;
      cur_r = cur_l + delta;
      walls(i, 1) = cur_l;
      walls(i, 2) = cur_r;
    endfor

    #f means fixed, nj >= 6
    f_cnt_in_bucket = [];
    f_walls = [];
    f_m = 0;
    for i = 1 : m
      if (i == 1 || f_cnt_in_bucket(f_m) >= 6)
        f_m = f_m + 1;
        f_walls(f_m, 1) = walls(i, 1);
        f_cnt_in_bucket(f_m) = 0;
      endif
      f_walls(f_m, 2) = walls(i, 2);
      f_cnt_in_bucket(f_m) = f_cnt_in_bucket(f_m) + cnt_in_bucket(i);
    endfor
    
    #Выборочное среднее
    E = mean(X);
    #Выборочная отклонение
    SQRT_D = std(X);

    P = [];
    for i = 1 : f_m
      P(i) = normcdf(f_walls(i, 2), E, SQRT_D) - normcdf(f_walls(i, 1), E, SQRT_D);
    endfor

    hi2 = sum(((f_cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 >= chi2inv(1 - alpha, f_m - 1 - 2));
  endfor
  printf("Normal distribution satisfies the hypothesis about normal distribution\n")
  printf("Sample grouped to met n_j >= 6\n")
  printf("For alpha = %d, probability of type I error is %d\n", alpha, res / tests)
endfunction

function res = test_Chi2_3(tests, n, m, d)
  res = 0;
  mu = 1;
  sigma = 1;
  alpha = 0.05;
  for t = 1 : tests
    X = normrnd(mu, sigma, n, 1);

    l = min(X);
    r = max(X);
    
    delta = (r - l) / m;
    cnt_in_bucket = hist(X, m);
    walls = [];
    for i = 1 : m
      cur_l = (i - 1) * delta + l;
      cur_r = cur_l + delta;
      walls(i, 1) = cur_l;
      walls(i, 2) = cur_r;
    endfor

    #f means fixed, nj >= 6
    f_cnt_in_bucket = [];
    f_walls = [];
    f_m = 0;
    for i = 1 : m
      if (i == 1 || f_cnt_in_bucket(f_m) >= 6)
        f_m = f_m + 1;
        f_walls(f_m, 1) = walls(i, 1);
        f_cnt_in_bucket(f_m) = 0;
      endif
      f_walls(f_m, 2) = walls(i, 2);
      f_cnt_in_bucket(f_m) = f_cnt_in_bucket(f_m) + cnt_in_bucket(i);
    endfor
    
    #Выборочное среднее
    E = mean(X) + d;
    #Выборочная отклонение
    SQRT_D = std(X) + d;

    P = [];
    for i = 1 : f_m
      P(i) = normcdf(f_walls(i, 2), E, SQRT_D) - normcdf(f_walls(i, 1), E, SQRT_D);
    endfor

    hi2 = sum(((f_cnt_in_bucket - n .* P) .^ 2) ./ (n .* P));
    res = res + (hi2 < chi2inv(1 - alpha, f_m - 1 - 2));
  endfor
  printf("Normal distribution satisfies the hypothesis about normal distribution\n")
  printf("Выборочное среднее и выборочная дисперсия изменены на %d\n", d)
  printf("For alpha = %d, probability of type II error is %d\n", alpha, res / tests)
endfunction

n = 10 ^ 6;
mu = 1;
sigma = 1;
m = 10 ^ 2;
tests = 10 ^ 3;

X = normrnd(mu, sigma, n, 1);

l = min(X);
r = max(X);
delta = (r - l) / m;
[y_coords x_coords] = hist(X, m);

x_coords_for_normpdf = l:0.1:r;
bar(x_coords, y_coords / (n * delta));
hold on;
plot(x_coords_for_normpdf, normpdf(x_coords_for_normpdf, mu, sigma));

printf("Sample size = %d\n", n)
printf("Length of intervals = %d\n", delta)
printf("Number of intervals = %d\n", m)

printf("\n")

# PART 2

test_Chi2_1(10 ^ 3, 10 ^ 4, m);

printf("\n")

test_Chi2_2(10 ^ 3, 10 ^ 4, m);

printf("\n")

test_Chi2_3(10 ^ 3, 10 ^ 4, m, 0.005); #0.005

printf("\n")

test_Chi2_3(10 ^ 3, 10 ^ 4, m, 0.05); #0.1

printf("\n")
