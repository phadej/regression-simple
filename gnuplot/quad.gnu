f(x) = alpha*x*x + beta *x + gamma
fit f(x) 'quad.dat' using 1:2 via alpha, beta, gamma

plot 'quad.dat' using 1:2, f(x)
