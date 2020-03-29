pkg load statistics

n = 5 * 10 ^ 5;
mu = 1;
sigma = 1;

X = normrnd(mu, sigma, n, 1);

l = min(X);
r = max(X);
buckets = ceil((r - l) * n ^ (1 / 3));
delta = (r - l) / buckets;
x_coords = zeros(buckets, 1);
y_coords = zeros(buckets, 1);

for i=1:buckets
  cur_l = i * delta + l;
  cur_r = cur_l + delta;
  x_coords(i) = (cur_r + cur_l) / 2;
  y_coords(i) = (sum(X <= cur_r) - sum(X < cur_l)) / (n * delta);
endfor

bar(x_coords, y_coords);

printf("Размер выборки = %d\n", n);
printf("Границы = [%d; %d]\n", l, r);
printf("Выбранная длина интервалов = %d\n", delta);
printf("Количество интервалов = %d\n", buckets);
