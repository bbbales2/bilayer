library(lattice)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 10
JN = 10
KN = 10
X = 0.05
Y = 0.07
Z = 0.018
density = 8700.0
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
c44 = 140
a = 2.8

c12 = -(c44 * 2.0 / a - c11)
cm1 = buildcm(c11, c12, c44)

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

cs = Cvoigt(cm1)

inner = function(i0, i1, j0, j1, k0, k1) {
  # int_0_X x^i0 * x^i1 dx
  ret = 0.0
  if(min(i0, i1, j0, j1, k0, k1) >= 0) {
    ret = (X^(i0 + i1 + 1) / (i0 + i1 + 1)) *
      (Y^(j0 + j1 + 1) / (j0 + j1 + 1)) *
      (Z^(k0 + k1 + 1) / (k0 + k1 + 1))
  }
  ret
}

idxs = list()
n = 1
for(i in 0:IN) {
  for(j in 0:JN) {
    for(k in 0:KN) {
      if((i + j + k) <= max(IN, JN, KN)) {
        idxs[[n]] = c(i, j, k)
        n = n + 1
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
    i1 = idxs[[n1]][1]
    j1 = idxs[[n1]][2]
    k1 = idxs[[n1]][3]
      
    dinp[n0, n1, 1, 1] = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1, k0, k1)
    dinp[n0, n1, 1, 2] = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1, k0, k1)
    dinp[n0, n1, 1, 3] = i0 * k1 * inner(i0 - 1, i1, j0, j1, k0, k1 - 1)
    dinp[n0, n1, 2, 1] = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1, k0, k1)
    dinp[n0, n1, 2, 2] = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1, k0, k1)
    dinp[n0, n1, 2, 3] = j0 * k1 * inner(i0, i1, j0 - 1, j1, k0, k1 - 1)
    dinp[n0, n1, 3, 1] = k0 * i1 * inner(i0, i1 - 1, j0, j1, k0 - 1, k1)
    dinp[n0, n1, 3, 2] = k0 * j1 * inner(i0, i1, j0, j1 - 1, k0 - 1, k1)
    dinp[n0, n1, 3, 3] = k0 * k1 * inner(i0, i1, j0, j1, k0 - 1, k1 - 1)
    inp[n0, n1] = inner(i0, i1, j0, j1, k0, k1)
  }
}

K = matrix(0, nrow = 3 * N, ncol = 3 * N)
M = matrix(0, nrow = 3 * N, ncol = 3 * N)
for(n0 in 1:N) {
  for(n1 in 1:N) {
    for(i in 1:3) {
      for(k in 1:3) {
        total = 0.0
        
        for(j in 1:3) {
          for(l in 1:3) {
            total = total + cs[i, j, k, l] * dinp[n0, n1, j, l]
          }
        }
        
        K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] = K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] + total
      }
      M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] + density * inp[n0, n1]
    }
  }
}

r = geigen(K, M, TRUE)
print(1e-3 * sqrt(r$values[7:16] * 1e9) / (pi * 2))

dK = array(0, c(3 * N, 3 * N, 6, 6))
for(ii in 1:6) {
  for(jj in 1:6) {
    cm2 = matrix(0, nrow = 6, ncol = 6)
    cm2[ii, jj] = 1.0
    cs2 = Cvoigt(cm2)
    for(n0 in 1:N) {
      for(n1 in 1:N) {
        for(i in 1:3) {
          for(k in 1:3) {
            total = 0.0
            
            for(j in 1:3) {
              for(l in 1:3) {
                total = total + cs2[i, j, k, l] * dinp[n0, n1, j, l]
              }
            }
            
            dK[3 * (n0 - 1) + i, 3 * (n1 - 1) + k, ii, jj] = dK[3 * (n0 - 1) + i, 3 * (n1 - 1) + k, ii, jj] + total
          }
        }
      }
    }
  }
}

U = chol(M)
dK2 = array(0, dim(dK))
for(i in 1:6) {
  for(j in 1:6) {
    P = solve(t(U), t(solve(t(U), t(dK[,, i, j]))))
    dK2[,,i, j] = (P + t(P)) / 2.0
  }
}

K2 = matrix(0, 3 * N, 3 * N)
for(i in 1:6) {
  for(j in 1:6) {
    K2 = K2 + cm1[i, j] * dK2[,,i, j]
  }
}

r2 = eigen(K2, TRUE)
print(1e-3 * sqrt((r2$values %>% sort)[7:14] * 1e9) / (pi * 2))
