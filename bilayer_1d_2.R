library(tidyverse)
library(ggplot2)
library(geigen)

N = 50
X = 1.0
xs = seq(0.0, X, length = N)
m = as.integer(0.2307692 * N)
densities = c(1.0, 2.0)
cs = c(1.0, 5.0)

K = matrix(0, nrow = N, ncol = N)
M = matrix(0, nrow = N, ncol = N)

for(n in 1:(N - 1)) {
  j = if(n <= m) 1 else 2

  dx = xs[n + 1] - xs[n]
  K[n, n] = K[n, n] + cs[j] / dx
  K[n, n + 1] = K[n, n + 1] - cs[j] / dx
  K[n + 1, n] = K[n, n + 1]
  K[n + 1, n + 1] = K[n + 1, n + 1] + cs[j] / dx
  
  M[n, n] = M[n, n] + densities[j] * (xs[n + 1] - xs[n]) / 3.0
  M[n, n + 1] = M[n, n + 1] + densities[j] * (xs[n + 1] - xs[n]) / 6.0
  M[n + 1, n] = M[n, n + 1]
  M[n + 1, n + 1] = M[n + 1, n + 1] + densities[j] * (xs[n + 1] - xs[n]) / 3.0
}

r = geigen(K, M, TRUE)

print(r$values[1:5])

plot(xs, r$vectors[,2])

#0.00000    18.38083    64.44047   145.21138   270.67956
#0.00000    37.08421   109.52220   295.85138   519.00046
#0.0000    37.1330   109.5271   296.9175   526.5467   779.4810  1214.9543  1680.6284  2343.0971  4953.5159  8777.8619 42772.0322

