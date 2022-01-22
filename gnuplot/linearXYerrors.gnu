g(x) = u*x + v
fit g(x) 'linear.dat' using 2:1:4:3 xyerror via u, v

f(x) = alpha*x + beta
fit f(x) 'linear.dat' using 1:2:3:4 xyerror via alpha, beta

plot 'linear.dat' using 1:2:3:4 with xyerror lt 1 title "Data", f(x) lt 1 title "Fit", 3 * x + 5 lt 2 title "True"

set key left top
set grid
replot
