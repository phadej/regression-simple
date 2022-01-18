f(x) = alpha*x*x + beta *x + gamma
fit f(x) 'quad.dat' using 1:2:3 yerror via alpha, beta, gamma

plot 'quad.dat' using 1:2:3 with yerror, f(x)
