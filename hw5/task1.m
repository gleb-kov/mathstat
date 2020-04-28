pkg load statistics
#{
function count_risks(mu, sigma, n, m)
  C1 = 0.4;
  X = sort(normrnd(mu, sigma, m, n));
  #mean(X)
  risk1 = mean(mean(X) - sigma ^ 2 / m)
  if mod(m, 2) == 1
    med = X(fix(m / 2) + 1, :);
  else
    med = (X(fix(m / 2), :) + X(fix(m / 2) + 1, :)) / 2;
  endif

  risk2 = mean(med - pi * sigma ^ 2 / (2 * m))

  risk3 = mean((X(1, :) + X(m, :)) / 2 - C1 * sigma ^ 2 / log(m))
endfunction
#}

function count_risks(mu, sigma, n, m)
  C1 = 0.4;
  std1 = std2 = std3 = [];
  for i = 1 : n
    X = sort(normrnd(mu, sigma, m, n));
    med = median(X);
    std1(i) = std(mean(X));
    std2(i) = std(med);
    std3(i) = std((X(1, :) + X(m, :)) / 2);
  endfor
  risk1 = mean((std1 - sigma / sqrt(m)) .^ 2)
  risk2 = mean((std2 - sqrt(pi) * sigma / sqrt(2 * m)) .^ 2)
  risk3 = mean((std3 - sqrt(C1) * sigma / sqrt(log(m))) .^ 2)
endfunction

count_risks(1, 3, 100, 100);  
count_risks(1, 3, 100, 10000);   
  