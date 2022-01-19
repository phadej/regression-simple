stats 'linear.dat'
f(x) = alpha*x + beta
fit f(x) 'linear.dat' using 1:2 via alpha, beta

plot 'linear.dat' using 1:2, f(x)
