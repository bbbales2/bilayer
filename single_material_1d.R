library(tidyverse)
library(ggplot2)
library(geigen)

N = 10
X = 1.0
density = 1.0
c = 1.0

inner = function(i0, i1) {
  # int_0_X x^i0 * x^i1 dx
  X^(i0 + i1) / (i0 + i1 + 1)
}

inner_dx = function(i0, i1) {
  ret = 0.0
  
  if(i0 == 0 || i1 == 0) {
    ret = 0.0
  } else {
    ret = i0 * i1 * inner(i0 - 1, i1 - 1)
  }
  
  ret
}

K = matrix(0, nrow = N, ncol = N)
M = matrix(0, nrow = N, ncol = N)

for(i in 1:N) {
  for(j in 1:N) {
    K[i, j] = c * inner_dx(i - 1, j - 1)
    M[i, j] = density * inner(i - 1, j - 1)
  }
}

r = geigen(K, M, TRUE)

xs = seq(0.0, X, length = 200)
map(1:N, ~ r$vectors[.,5] * xs^(. - 1)) %>%
  Reduce('+', .) %>% plot

