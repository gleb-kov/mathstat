pkg load statistics

mu = -1;
sigma = 0.5;
n = 10^4;
m = 100;
t0 = -0.25;
y_true = normcdf(t0, mu, sigma);
gamma = 0.95;
T = norminv((gamma + 1) / 2);
X = normrnd(mu, sigma, n, m);
p = mean(X < t0);

delta = T * sqrt(p .* (1 - p) / n);
left = p - delta;
right = p + delta;
correct = (left <= y_true) & (y_true <= right);
correct = sum(correct) / m;

printf("N = %d\n", n);
printf("Real value of Normal distribution function = %g\n", y_true);
printf("Correct simulation = %g\n", correct);

y_vector = repmat(y_true, m);
x = 1 : m
plot(x, left, x, right, x, y_vector)