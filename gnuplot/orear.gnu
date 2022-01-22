alpha = 1e-3
beta = 6e5
f(x) = alpha*x - beta/x
fit f(x) 'orear.dat' using 1:2:4 yerror via alpha, beta
fit f(x) 'orear.dat' using 1:2:3:4 xyerror via alpha, beta

plot 'orear.dat' using 1:2:3:4 with xyerror lt 1 title "Data", f(x) lt 1 title "Fit"

set key left top
set grid
replot
