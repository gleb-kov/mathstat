// TODO
pkg load statistics

n = 100
mu =  1
sigma =  1
X = sort(normrnd(mu, sigma, n, 1))

F_n = 1 / n : 1 / n : 1
[a, b] = stairs(X, F_n)

y = 0.95
u = 1.36
