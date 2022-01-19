f(x) = alpha*x*x + beta *x + gamma
fit f(x) 'quad.dat' using 1:2:3:4 xyerror via alpha, beta, gamma

plot 'quad.dat' using 1:2:3:4 with xyerror, f(x), 0.1 * x * x - 3 * x + 5 title "True"

set key left top
set grid
replot
