library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 5
JN = 4
KN = 3
LN = 3
X = 0.2
Y = 0.3
Z = 0.4
B = Z /2

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
Cvoigt = function(cm) {
  C = array(0, c(3, 3, 3, 3))
  
  voigt = list(list(c(1, 1)),
               list(c(2, 2)),
               list(c(3, 3)),
               list(c(2, 3), c(3, 2)),
               list(c(1, 3), c(3, 1)),
               list(c(1, 2), c(2, 1)))
  
  for(i in 1:6) {
    for(j in 1:6) {
      for(vi in voigt[[i]]) {
        k = vi[1]
        l = vi[2]
        for(vj in voigt[[j]]) {
          n = vj[1]
          m = vj[2]
          C[k, l, n, m] = cm[i, j]
        }
      }
    }
  }
  
  C
}

cm1 = buildcm(250, 150, 140)
cm2 = buildcm(150, 100, 40)

cs = list(Cvoigt(cm1), Cvoigt(cm1))

inner = function(i0, i1, j0, j1, k0, k1, l0, l1, dx0, dx1, dy0, dy1, dz0, dz1) {
  # int_0_X x^i0 * x^i1 dx
  ret = 0.0

  pk = ifelse(dz0, k0, 1) * ifelse(dz1, k1, 1)
  pl0 = ifelse(dz0, l0, 1)
  pl1 = ifelse(dz1, l1, 1)
  pl = pl0 * pl1
  B0 = ifelse(dz0, 0, B)
  B1 = ifelse(dz1, 0, B)
  k0 = ifelse(dz0, k0 - 1, k0)
  k1 = ifelse(dz1, k1 - 1, k1)
  l0 = ifelse(dz0, l0 - 1, l0)
  l1 = ifelse(dz1, l1 - 1, l1)

  px = ifelse(dx0, i0, 1) * ifelse(dx1, i1, 1)
  i0 = ifelse(dx0, i0 - 1, i0)
  i1 = ifelse(dx1, i1 - 1, i1)

  py = ifelse(dy0, j0, 1) * ifelse(dy1, j1, 1)
  j0 = ifelse(dy0, j0 - 1, j0)
  j1 = ifelse(dy1, j1 - 1, j1)

  pre = px * (X^(i0 + i1 + 1) / (i0 + i1 + 1)) *
    py * (Y^(j0 + j1 + 1) / (j0 + j1 + 1))

  if(min(i0, i1, j0, j1, k0, k1) >= 0) {
    ret = ret + pre *
      pk * (B^(k0 + k1 + 1) / (k0 + k1 + 1))
  }

  if(min(i0, i1, j0, j1, k0, k1) >= 0) {
    ret = ret + pre *
      (Z - B) * B0^k0 * B1^k1
  }

  if(min(i0, i1, j0, j1, k1, l0) >= 0) {
    ret = ret + pre *
      pl0 * B1^k1 * (Z^(l0 + 1) - B^(l0 + 1)) / (l0 + 1)
  }

  if(min(i0, i1, j0, j1, k0, l1) >= 0) {
    ret = ret + pre *
      pl1 * B0^k0 * (Z^(l1 + 1) - B^(l1 + 1)) / (l1 + 1)
  }

  if(min(i0, i1, j0, j1, l0, l1) >= 0) {
    ret = ret + pre *
      pl * (Z^(l0 + l1 + 1) - B^(l0 + l1 + 1)) / (l0 + l1 + 1)
  }

  ret
}

idxs = list()
n = 1
for(i in 0:IN) {
  for(j in 0:JN) {
    for(k in 0:KN) {
      for(l in 0:LN) {
        if((i + j + max(k, l)) <= max(IN, JN, KN, LN)) {
          idxs[[n]] = c(i, j, k, l)
          n = n + 1
        }
      }
    }
  }
}
N = length(idxs)

dinp = array(0, c(N, N, 3, 3))
inp = array(0, c(N, N))

for(n0 in 1:N) {
  for(n1 in 1:N) {
    i0 = idxs[[n0]][1]
    j0 = idxs[[n0]][2]
    k0 = idxs[[n0]][3]
    l0 = idxs[[n0]][4]
    i1 = idxs[[n1]][1]
    j1 = idxs[[n1]][2]
    k1 = idxs[[n1]][3]
    l1 = idxs[[n1]][4]
      
    dinp[n0, n1, 1, 1] = inner(i0, i1, j0, j1, k0, k1, l0, l1, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    dinp[n0, n1, 1, 2] = inner(i0, i1, j0, j1, k0, k1, l0, l1, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)
    dinp[n0, n1, 1, 3] = inner(i0, i1, j0, j1, k0, k1, l0, l1, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
    dinp[n0, n1, 2, 1] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
    dinp[n0, n1, 2, 2] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
    dinp[n0, n1, 2, 3] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)
    dinp[n0, n1, 3, 1] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
    dinp[n0, n1, 3, 2] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
    dinp[n0, n1, 3, 3] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
    inp[n0, n1] = inner(i0, i1, j0, j1, k0, k1, l0, l1, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  }
}

K = matrix(0, nrow = 3 * N, ncol = 3 * N)
M = matrix(0, nrow = 3 * N, ncol = 3 * N)
for(ii in 1:length(densities)) {
  for(n0 in 1:N) {
    #cat(n0, "\n");
    for(n1 in 1:N) {
      k0 = idxs[[n0]][3]
      k1 = idxs[[n1]][3]
      
      for(i in 1:3) {
        for(k in 1:3) {
          total = 0.0
          
          for(j in 1:3) {
            for(l in 1:3) {
              total = total + cs[[ii]][i, j, k, l] * dinp[n0, n1, j, l]
            }
          }
          
          K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] = K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] + total
        }
        M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] + densities[ii] * inp[n0, n1]
      }
    }
  }
}

r = geigen(K, M, TRUE)

print(r$values[1:25])
print(sqrt(r$values[7:14] * 1e9) / (pi * 2))

# 0.000000e+00 0.000000e+00 0.000000e+00 3.226860e-15 2.765512e-14 2.086823e-13 5.138824e-01 6.912472e-01 7.008247e-01 8.574773e-01 1.308624e+00 1.376295e+00 1.404012e+00 1.738664e+00
# 1.757523e+00 1.963727e+00 2.119848e+00 2.939276e+00 3.012592e+00 3.089201e+00 3.411555e+00 3.470279e+00 3.533797e+00 3.577868e+00 3.744839e+00
#> print(sqrt(r$values[7:14] * 1e9) / (pi * 2))
# [1] 3607.879 4184.435 4213.324 4660.489 5757.414 5904.401 5963.560 6636.331

