f(x) = alpha*x + beta
fit f(x) 'linear.dat' using 1:2 via alpha, beta
stats 'linear.dat'

plot 'linear.dat' using 1:2, f(x)
