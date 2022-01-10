# > load 'example.gnu'
# Results should be the same as the library produces.

f(x) = alpha*x + beta
fit f(x) 'example.dat' via alpha, beta

g(x) = alpha*x*x + beta * x + gamma
fit g(x) 'example.dat' via alpha, beta, gamma
