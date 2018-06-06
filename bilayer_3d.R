library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 10
JN = 10
KN = 10
X = 0.7
Y = 1.2
Z = 0.3
densities = c(1.0, 1.0)
cm = diag(6)
cm[1, 2] = 0.1
cm[2, 1] = 0.1
cm[1, 3] = 0.1
cm[3, 1] = 0.1
cm[2, 3] = 0.1
cm[3, 2] = 0.1
B = 0.2307692

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

cs = list(Cvoigt(cm), Cvoigt(cm))

inner = function(i0, i1, j0, j1, k0, k1, a, b) {
  # int_0_X x^i0 * x^i1 dx
  ret = 0.0
  if(min(i0, i1, j0, j1, k0, k1) >= 0) {
    ret = (X^(i0 + i1 + 1) / (i0 + i1 + 1)) *
      (Y^(j0 + j1 + 1) / (j0 + j1 + 1)) *
      (b^(k0 + k1 + 1) / (k0 + k1 + 1) - a^(k0 + k1 + 1) / (k0 + k1 + 1))
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
    
    dinp[n0, n1, 1, 1] = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1, k0, k1, 0, Z)
    dinp[n0, n1, 1, 2] = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1, k0, k1, 0, Z)
    dinp[n0, n1, 1, 3] = i0 * k1 * inner(i0 - 1, i1, j0, j1, k0, k1 - 1, 0, Z)
    dinp[n0, n1, 2, 1] = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1, k0, k1, 0, Z)
    dinp[n0, n1, 2, 2] = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1, k0, k1, 0, Z)
    dinp[n0, n1, 2, 3] = j0 * k1 * inner(i0, i1, j0 - 1, j1, k0, k1 - 1, 0, Z)
    dinp[n0, n1, 3, 1] = k0 * i1 * inner(i0, i1 - 1, j0, j1, k0 - 1, k1, 0, Z)
    dinp[n0, n1, 3, 2] = k0 * j1 * inner(i0, i1, j0, j1 - 1, k0 - 1, k1, 0, Z)
    dinp[n0, n1, 3, 3] = k0 * k1 * inner(i0, i1, j0, j1, k0 - 1, k1 - 1, 0, Z)
    inp[n0, n1] = inner(i0, i1, j0, j1, k0, k1, 0, Z)
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
            total = total + cs[[1]][i, j, k, l] * dinp[n0, n1, j, l]
          }
        }
        
        K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] = total
      }
      M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = densities[1] * inp[n0, n1]
    }
  }
}

r = geigen(K, M, TRUE)

print(r$values[7:14])

{
  xs = seq(0.0, X, length = 20)
  y = 0.5
  zs = seq(0.0, Z, length = 20)
  u = matrix(0, nrow = length(xs), ncol = length(zs))
  for(n in 1:N) {
    i = idxs[[n]][1]
    j = idxs[[n]][2]
    k = idxs[[n]][3]
    u = u + outer(xs, zs, function(x, z) { r$vectors[2 * (n - 1) + 1, 1] * (x^i) * (y^j) * (z^k) })
  }
  levelplot(u)
}

# 7.442681e+00  8.283205e+00  8.283205e+00  8.613477e+00  8.750229e+00  8.981048e+00  9.263764e+00 9.934843e+00  9.934843e+00  1.038664e+01  1.090038e+01  1.166401e+01  1.222080e+01  1.235710e+01