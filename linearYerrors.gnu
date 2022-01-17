f(x) = alpha*x + beta
fit f(x) 'linear.dat' using 1:2:3 yerror via alpha, beta

plot 'linear.dat' using 1:2:3 with yerror, f(x)
