library(tidyverse)
library(ggplot2)
library(geigen)
library(Ryacas)

N = 3
X = 1.0
xs = seq(0.0, X, length = N)
m = as.integer(0.2307692 * N)
densities = c(1.0, 2.0)
cs = c(1.0, 5.0)

x = Sym("x")
fs = list("x", "1 - x")#list("1 - 3 * x + 2 * x * x", "4 * x - 4 * x * x", "-x + 2 * x * x")
bp = array(0, c(length(fs), 2, length(fs), 2))
for(i in 1:length(fs)) {
  for(di in 1:2) {
    for(j in 1:length(fs)) {
      for(dj in 1:2) {
        f = Integrate(paste0('(', ifelse(di == 2, "D(x) ", ""), fs[[i]], ') * (', ifelse(dj == 2, "D(x) ", ""), fs[[j]] ,')'), x) %>% as.expression
        bp[i, di, j, dj] = eval(f, list(x = 1.0)) - eval(f, list(x = 0.0))
        #cat(paste(i, di, j, dj, bp[i, di, j, dj], "\n"))
      }
    }
  }
}

Nf = (length(fs) - 1) * (N - 1) + 1
K = matrix(0, nrow = Nf, ncol = Nf)
M = matrix(0, nrow = Nf, ncol = Nf)

for(n in 1:(N - 1)) {
  j = if(n <= m) 1 else 2
  
  dx = xs[n + 1] - xs[n]
  r = ((length(fs) - 1) * (n - 1) + 1) : ((length(fs) - 1) * n + 1)
  K[r, r] = K[r, r] + cs[j] * bp[, 2,, 2] / dx
  
  M[r, r] = M[r, r] + densities[j] * bp[, 1,, 1] * dx
}

r = geigen(K, M, TRUE)

print(r$values[1:5])

plot(seq(0.0, X, length = dim(K)[1]), r$vectors[,2])

#0.00000    18.38083    64.44047   145.21138   270.67956
#0.00000    37.08421   109.52220   295.85138   519.00046
#0.0000    37.1330   109.5271   296.9175   526.5467   779.4810  1214.9543  1680.6284  2343.0971  4953.5159  8777.8619 42772.0322

