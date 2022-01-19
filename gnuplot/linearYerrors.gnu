f(x) = alpha*x + beta
fit f(x) 'linear.dat' using 1:2:4 yerror via alpha, beta

plot 'linear.dat' using 1:2:4 with yerror lt 1 title "Data", f(x) lt 1 title "Fit", 3 * x + 5 lt 2 title "True"

set key left top
set grid
replot
