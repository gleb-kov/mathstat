pkg load statistics

function count_risks(a, u, m, n)
  C2 = 0.9;
  X = sort(a + exprnd(u, n, m) - exprnd(u, n, m));
  med = median(X);
  
  std1_practical = std(mean(X))
  std1_theoretical = sqrt(2) * u / sqrt(n)
 
  std2_practical = std(med)
  std2_theoretical = u / sqrt(n)
  
  std3_practical = std((X(1, :) + X(n, :)) / 2)
  std3_theoretical = sqrt(C2) * u
  printf('\n');
endfunction

count_risks(1, 10, 100, 100);
count_risks(1, 10, 100, 10000);