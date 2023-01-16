f(x) = alpha * x * log(x)
alpha = 1
fit f(x) 'issue8.dat' using 1:2 via alpha

plot 'issue8.dat' using 1:2, f(x)
