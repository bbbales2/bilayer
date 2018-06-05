library(tidyverse)
library(ggplot2)
library(geigen)

N = 20
X = 1.0
densities = c(1.0, 1.0)
cs = c(1.0, 2.0)
B = 0.2307692

inner = function(i0, i1, a, b) {
  # int_0_X x^i0 * x^i1 dx
  b^(i0 + i1 + 1) / (i0 + i1 + 1) - a^(i0 + i1 + 1) / (i0 + i1 + 1)
}

inner_dx = function(i0, i1, a, b) {
  ret = 0.0
  
  if(i0 == 0 || i1 == 0) {
    ret = 0.0
  } else {
    ret = i0 * i1 * inner(i0 - 1, i1 - 1, a, b)
  }
  
  ret
}

K = matrix(0, nrow = N, ncol = N)
M = matrix(0, nrow = N, ncol = N)

for(i in 1:N) {
  for(j in 1:N) {
    K[i, j] = cs[1] * inner_dx(i - 1, j - 1, 0, B) + cs[2] * inner_dx(i - 1, j - 1, B, X)
    M[i, j] = densities[1] * inner(i - 1, j - 1, 0, B) + densities[2] * inner(i - 1, j - 1, B, X)
  }
}

r = geigen(K, M, TRUE)

print(r$values)

xs = seq(0.0, X, length = 200)
map(1:N, ~ r$vectors[.,3] * xs^(. - 1)) %>%
  Reduce('+', .) %>% plot
#4.292280e-11 9.870433e+00 3.949167e+01 8.889355e+01 1.581258e+02
