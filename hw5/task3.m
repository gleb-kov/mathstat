pkg load statistics

function count_risks(a, u, n, m)
  C2 = 0.9;
  std1 = std2 = std3 = [];
  for i = 1 : n
    X = sort(laplace_rnd(m, n));
    % X += u ???
    med = median(X);
    std1(i) = std(mean(X));
    std2(i) = std(med);
    std3(i) = std((X(1, :) + X(m, :)) / 2);
  endfor
  risk1 = mean((std1 -  u * sqrt(2) / sqrt(m)) .^ 2);
  risk2 = mean((std2 - u / sqrt(m)) .^ 2);
  risk3 = mean((std3 - sqrt(C2) * u) .^ 2);
  printf("Risk1 = %d, risk2 = %d, risk3 = %d\n", risk1, risk2, risk3)
endfunction

count_risks(1, 3, 100, 100);  
count_risks(1, 3, 100, 10000);  