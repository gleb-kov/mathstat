pkg load statistics

function count_risks(a, sigma, m, n)
  C1 = 0.4;
  X = sort(normrnd(a, sigma, n, m));
  med = median(X);
  
  std1_practical = std(mean(X))
  std1_theoretical = sigma / sqrt(n)
 
  std2 = std(med)
  std2_theoretical = sqrt(pi) * sigma / sqrt(2 * n)
  
  std3 = std((X(1, :) + X(n, :)) / 2)
  std3_theoretical = C1 * sigma / sqrt(log(n))
  printf('\n');
endfunction

count_risks(1, 3, 100, 100);
count_risks(1, 3, 100, 10000);
