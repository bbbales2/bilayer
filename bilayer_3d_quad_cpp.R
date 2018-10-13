library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)
library(Ryacas)
library(RcppEigen)
library(Rcpp11)
library(Rcpp)
library(RSpectra)
sourceCpp("rus_dev.cpp")

IN = 5
JN = 5
KN = 7
X = 0.07
Y = 0.05
Z = 0.018
B = 0.014
zs = seq(0.0, Z, length = KN)
i = 1
while(zs[i] < B) {
  i = i + 1
}
zint = i
zs = c(zs[1:(i - 1)], B, zs[i:KN])
KN = length(zs)

M = calc_M(IN, JN, KN)
dKhatdcijii = array(build_dKdci_M(X, Y, zs, IN, JN, zint, densities), c(3 * M, 3 * M, 6, 6, 2))

densities = c(8700, 8700)

buildcm = function(c11, c12, c44) {
  cm = matrix(0, nrow = 6, ncol = 6)
  cm[1, 1] = c11
  cm[2, 2] = cm[1, 1]
  cm[3, 3] = cm[1, 1]
  cm[1, 2] = c12
  cm[2, 1] = cm[1, 2]
  cm[1, 3] = cm[1, 2]
  cm[2, 3] = cm[1, 2]
  cm[3, 1] = cm[1, 2]
  cm[3, 2] = cm[1, 2]
  cm[4, 4] = c44
  cm[5, 5] = cm[4, 4]
  cm[6, 6] = cm[4, 4]
  
  cm
}

cm2 = buildcm(250, 150, 140)
cm1 = buildcm(150, 100, 40)

cms = list(cm1, cm2)

Kc = matrix(0, 3 * M, 3 * M)
for(ii in 1:length(densities)) {
  for(i in 1:6) {
    for(j in 1:6) {
      Kc = Kc + cms[[ii]][i, j] * dKhatdcijii[,,i, j, ii]
    }
  }
}

nnz = sum(Kc != 0.0)
is = rep(0, nnz)
js = rep(0, nnz)
vs = rep(0, nnz)
n = 1
for(i in 1:nrow(Kc)) {
  for(j in 1:nrow(Kc)) {
    if(Kc[i, j] != 0.0) {
      is[n] = i - 1
      js[n] = j - 1
      vs[n] = Kc[i, j]
      n = n + 1
    }
  }
}

system.time(rev(eigen(Kc, TRUE)$values)[1:25])

system.time(calc_evals(Kc)[1:25])
system.time(calc_evals_spectra(Kc, 25, 0))
system.time(calc_evals_sparse(is, js, vs, nrow(Kc), 25, 0))
calc_evals(Kc)[1:25]
sort(calc_evals_spectra(Kc, 25, 0))
sort(calc_evals_sparse(is, js, vs, nrow(Kc), 25, 0))

print(c(IN, JN, KN, sqrt(calc_evals(Kc)[7:14] * 1e9) / (pi * 2)))

# 6.00 6.00 8.00 10335.06 11002.16 17565.35 19192.09 21049.15 21236.04 22134.66 26321.04
# 6.00 6.00 5.00 10323.04 10993.06 17569.18 19201.36 21050.77 21235.02 22138.73 26330.65
# 6.00 6.00 4.00 10326.18 11002.58 17569.73 19204.38 21063.43 21235.39 22149.41 26342.29
# 6.00 6.00 3.00 10332.87 11023.64 17571.69 19224.54 21125.32 21237.33 22199.95 26405.96
# 5.00 5.00 5.00 10355.02 11022.51 17677.25 19240.18 21237.74 21250.58 22288.96 26759.38
# 4.00 4.00 5.00 10415.65 11116.00 17702.84 19414.46 21311.04 22756.69 24091.98 28508.34