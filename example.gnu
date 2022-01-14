# > load 'example.gnu'
# Results should be the same as the library produces.

f(x) = alpha*x + beta
fit f(x) 'example.dat' via alpha, beta

g(x) = alpha*x*x + beta * x + gamma
fit g(x) 'example.dat' via alpha, beta, gamma

ff(x) = alpha*x + beta
fit ff(x) 'example.dat' using 1:2:3 yerror via alpha, beta

gg(x) = alpha*x*x + beta * x + gamma
fit gg(x) 'example.dat' using 1:2:3 yerror via alpha, beta, gamma
