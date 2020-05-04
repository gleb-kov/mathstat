pkg load statistics

function count_risks(a, delta, m, n)
  X = sort(unifrnd(a - delta / 2, a + delta / 2, n, m));
  med = median(X);
  
  std1_practical = std(mean(X))
  std1_theoretical = delta / sqrt(12 * n)
 
  std2_practical = std(med)
  std2_theoretical = delta / sqrt(4 * n)
  
  std3_practical = std((X(1, :) + X(n, :)) / 2)
  std3_theoretical = delta / (sqrt(2) * n)
  printf('\n');
endfunction

count_risks(1, 10, 100, 100);
count_risks(1, 10, 100, 10000);