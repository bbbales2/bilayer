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
KN = 5
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

densities = c(8700, 8700)

M = calc_M(IN, JN, KN)
dKhatdcijii = array(build_dKdci_M(X, Y, zs, IN, JN, zint, densities), c(3 * M, 3 * M, 6, 6, 2))

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
c11 = 250
a = 2.8
c44 = 140
c12 = -(c44 * 2.0 / a - c11)
cm2 = buildcm(250, 150, 140)
cm1 = buildcm(250, 150, 140)

#cm2 = 100 * matrix(sapply(strsplit("3.69762	-0.429256	-0.86379	-0.337691	-1.55673	0.0315084	-0.429256	3.24644	-0.265911	0.186542	0.00307126	0.105072	-0.86379	-0.265911	4.14444	-0.289584	0.0512344	-0.835652	-0.337691	0.186542	-0.289584	4.35637	0.172839	-1.02667	-1.55673	0.00307126	0.0512344	0.172839	3.96744	-0.147433	0.0315084	0.105072	-0.835652	-1.02667	-0.147433	4.00015", "[ \t]"), as.numeric), nrow = 6)
#cm1 = 100 * matrix(sapply(strsplit("5.77877	-1.57676	1.51418	-1.38104	1.60109	0.162966	-1.57676	4.6482	-1.46722	2.83411	0.223867	-2.81944	1.51418	-1.46722	4.29902	-2.31215	-2.39193	1.49025	-1.38104	2.83411	-2.31215	4.55989	0.553293	-1.33558	1.60109	0.223867	-2.39193	0.553293	4.28768	0.638419	0.162966	-2.81944	1.49025	-1.33558	0.638419	9.8315", "[ \t]"), as.numeric), nrow = 6)

buildKc = function(cm1, cm2) {
  cms = list(cm1, cm2)
  Kc = matrix(0, 3 * M, 3 * M)
  for(ii in 1:length(densities)) {
    for(i in 1:6) {
      for(j in 1:6) {
        Kc = Kc + cms[[ii]][i, j] * dKhatdcijii[,,i, j, ii]
      }
    }
  }

  Kc  
}


dx = 1e-4
cm1m = cm1
cm1m[1, 1] = cm1m[1, 1] - dx
cm1p = cm1
cm1p[1, 1] = cm1p[1, 1] + dx

print((c(sqrt(calc_evals(buildKc(cm1p, cm2))[7:16] * 1e9) / (pi * 2)) - c(sqrt(calc_evals(buildKc(cm1m, cm2))[7:16] * 1e9) / (pi * 2))) / (2 * dx))

dx = 1e-1
cm1m = cm1
cm1m[1, 2] = cm1m[1, 2] - dx
cm1m[2, 1] = cm1m[1, 2]
cm1p = cm1
cm1p[1, 2] = cm1p[1, 2] + dx
cm1p[2, 1] = cm1p[1, 2]

print((c(sqrt(calc_evals(buildKc(cm1p, cm2))[7:16] * 1e9) / (pi * 2)) - c(sqrt(calc_evals(buildKc(cm1m, cm2))[7:16] * 1e9) / (pi * 2))) / (2 * dx))

#print((c(sqrt(calc_evals(buildKc(cm1, cm2p))[7:16] * 1e9) / (pi * 2)) - c(sqrt(calc_evals(buildKc(cm1, cm2m))[7:16] * 1e9) / (pi * 2))) / (2 * dx))

Kc = buildKc(cm1, cm2)

system.time(rev(eigen(Kc, TRUE)$values)[1:25])

system.time(calc_evals(Kc)[1:25])
system.time(calc_evals_spectra(Kc, 25, 0))
system.time(calc_evals_sparse(is, js, vs, nrow(Kc), 25, 0))
calc_evals(Kc)[1:25]
sort(calc_evals_spectra(Kc, 25, 0))
sort(calc_evals_sparse(is, js, vs, nrow(Kc), 25, 0))

print(c(sqrt(calc_evals(Kc)[7:16] * 1e9) / (pi * 2)))
print(c(sqrt(sort(calc_evals_sparse(is, js, vs, nrow(Kc), 25, 0))[7:14] * 1e9) / (pi * 2)))

# 6.00 6.00 8.00 10335.06 11002.16 17565.35 19192.09 21049.15 21236.04 22134.66 26321.04
# 6.00 6.00 5.00 10323.04 10993.06 17569.18 19201.36 21050.77 21235.02 22138.73 26330.65
# 6.00 6.00 4.00 10326.18 11002.58 17569.73 19204.38 21063.43 21235.39 22149.41 26342.29
# 6.00 6.00 3.00 10332.87 11023.64 17571.69 19224.54 21125.32 21237.33 22199.95 26405.96
# 5.00 5.00 5.00 10355.02 11022.51 17677.25 19240.18 21237.74 21250.58 22288.96 26759.38
# 4.00 4.00 5.00 10415.65 11116.00 17702.84 19414.46 21311.04 22756.69 24091.98 28508.34