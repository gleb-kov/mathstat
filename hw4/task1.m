pkg load statistics

clc;
clear all;

function res = f(x)
  res = exp(- x ^ 2 / 2);
endfunction

function res = laplace(x)
  res = quad(@f, 0, x) * 1 / sqrt(2 * pi);
endfunction

n = 10 ^ 5;
mu = 1;
sigma = 1;
alpha = 0.95;
real_hi2 = 100734.7;

X = sort(normrnd(mu, sigma, n, 1));

l = min(X);
r = max(X);
buckets = ceil((r - l) * n ^ (1 / 3));
delta = (r - l) / buckets;
x_coords = zeros(buckets, 1);
cnt_in_bucket = zeros(buckets, 1);
y_coords = zeros(buckets, 1);
walls = zeros(buckets, 1);

for i = 1 : buckets
  cur_l = (i - 1) * delta + l;
  cur_r = cur_l + delta;
  walls(i) = cur_r;
  x_coords(i) = (cur_r + cur_l) / 2;
  cnt_in_bucket(i) = sum(X <= cur_r) - sum(X < cur_l);
  y_coords(i) = cnt_in_bucket(i) / (n * delta);
endfor

x_coords_for_normpdf = -4:0.1:6;
bar(x_coords, y_coords);
hold on;
plot(x_coords_for_normpdf, normpdf(x_coords_for_normpdf, mu, sigma), "linewidth", 1);

printf("Размер выборки = %d\n", n);
printf("Границы = [%d; %d]\n", l, r);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", buckets);

printf("\n");

# PART 2

E = sum(x_coords .* cnt_in_bucket) / n;
D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
FIXED_D = n / (n - 1) * D;
SQ_D = sqrt(FIXED_D);

printf("Предполагаемое матожидание = %d\n", E);
printf("Предполагаемая дисперсия = %d\n", FIXED_D);

# suffix P means "Предполагаемое"

z = (walls - E) / SQ_D;

P = zeros(buckets, 1);

P(1) = laplace(z(1)) - (-0.5);
for i = 2 : buckets
  P(i) = laplace(z(i)) - laplace(z(i - 1));
endfor

nP = P * n;

hi2 = sum(((cnt_in_bucket - nP) .^ 2) ./ nP);

printf("Полученное значение Хи квадрат = %d\n", hi2);
printf("Теоретическое значение Хи квадрат = %d\n", real_hi2);
printf("Заданное распределение есть равномерное = %d\n", hi2 < real_hi2);

