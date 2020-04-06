pkg load statistics

clc;
clear all;

n = 10 ^ 5;
a = 20;
b = 80;
alpha = 0.95;
real_hi2 = 100734.7;

X = sort(unifrnd(a, b, 1, n));

# PART 1

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

real_y = 1 / (b - a);
bar(x_coords, y_coords);
hold on;
plot([a b], [real_y real_y], "linewidth", 1);

printf("Размер выборки = %d\n", n);
printf("Границы = [%d; %d]\n", l, r);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", buckets);

printf("\n");

# PART 2

E = sum(x_coords .* cnt_in_bucket) / n;
D = sum((x_coords - E) .^ 2 .* cnt_in_bucket) / n;
SQ_D = sqrt(D);

printf("Предполагаемое матожидание = %d\n", E);
printf("Предполагаемая дисперсия = %d\n", D);

# suffix P means "Предполагаемое"
aP = E - sqrt(3) * SQ_D;
bP = E + sqrt(3) * SQ_D;
fP = 1 / (bP - aP);

printf("Предполагаемые границы = [%d; %d]\n", aP, bP);
printf("Предполагаемое значение функции = %d\n", fP);

nP = zeros(buckets, 1);
nP(1) = n * fP * (walls(1) - aP);
for i = 1 : buckets - 2
  nP(i + 1) = n * fP * (walls(i + 1) - walls(i));
endfor
nP(buckets) = n * fP * (bP - walls(buckets - 1));

hi2 = sum(((cnt_in_bucket - nP) .^ 2) ./ nP);

printf("Полученное значение Хи квадрат = %d\n", hi2);
printf("Теоретическое значение Хи квадрат = %d\n", real_hi2);
printf("Заданное распределение есть равномерное = %d\n", hi2 < real_hi2);
