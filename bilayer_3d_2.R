library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 5
JN = 5
KN = 11
X = 0.7
Y = 1.2
Z = 1.0
zs = seq(0.0, Z, length = KN)
densities = c(1.0, 2.0)
cm = diag(6)
cm[1, 2] = 0.1
cm[2, 1] = 0.1
cm[1, 3] = 0.1
cm[3, 1] = 0.1
cm[2, 3] = 0.1
cm[3, 2] = 0.1
B = 0.8
m = which.min(abs(B - zs))

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

inner = function(i0, i1, j0, j1) {
  # int_0_X x^i0 * x^i1 dx
  ret = 0.0
  if(min(i0, i1, j0, j1) >= 0) {
    ret = (X^(i0 + i1 + 1) / (i0 + i1 + 1)) *
      (Y^(j0 + j1 + 1) / (j0 + j1 + 1))
  }
  ret
}

idxs = list()
n = 1
for(i in 0:IN) {
  for(j in 0:JN) {
    for(k in 1:KN) {
      if((i + j) <= max(IN, JN)) {
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
    
    if(k0 != k1) {
      next
    }
    
    if(k0 == KN) {
      next
    }

    tmp = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
    dinp[n0, n1, 1, 1] = dinp[n0, n1, 1, 1] + tmp
    dinp[n0 + 1, n1, 1, 1] = dinp[n0 + 1, n1, 1, 1] + tmp / 2.0
    dinp[n0, n1 + 1, 1, 1] = dinp[n0, n1 + 1, 1, 1] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 1, 1] = dinp[n0 + 1, n1 + 1, 1, 1] + tmp

    tmp = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
    dinp[n0, n1, 1, 2] = dinp[n0, n1, 1, 2] + tmp
    dinp[n0 + 1, n1, 1, 2] = dinp[n0 + 1, n1, 1, 2] + tmp / 2.0
    dinp[n0, n1 + 1, 1, 2] = dinp[n0, n1 + 1, 1, 2] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 1, 2] = dinp[n0 + 1, n1 + 1, 1, 2] + tmp

    tmp = i0 * inner(i0 - 1, i1, j0, j1) # f0df1
    dinp[n0, n1, 1, 3] = dinp[n0, n1, 1, 3] - tmp / 2.0
    dinp[n0 + 1, n1, 1, 3] = dinp[n0 + 1, n1, 1, 3] - tmp / 2.0
    dinp[n0, n1 + 1, 1, 3] = dinp[n0, n1 + 1, 1, 3] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 1, 3] = dinp[n0 + 1, n1 + 1, 1, 3] + tmp / 2.0

    tmp = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
    dinp[n0, n1, 2, 1] = dinp[n0, n1, 2, 1] + tmp
    dinp[n0 + 1, n1, 2, 1] = dinp[n0 + 1, n1, 2, 1] + tmp / 2.0
    dinp[n0, n1 + 1, 2, 1] = dinp[n0, n1 + 1, 2, 1] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 2, 1] = dinp[n0 + 1, n1 + 1, 2, 1] + tmp

    tmp = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
    dinp[n0, n1, 2, 2] = dinp[n0, n1, 2, 2] + tmp
    dinp[n0 + 1, n1, 2, 2] = dinp[n0 + 1, n1, 2, 2] + tmp / 2.0
    dinp[n0, n1 + 1, 2, 2] = dinp[n0, n1 + 1, 2, 2] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 2, 2] = dinp[n0 + 1, n1 + 1, 2, 2] + tmp

    tmp = j0 * inner(i0, i1, j0 - 1, j1) # f0df1
    dinp[n0, n1, 2, 3] = dinp[n0, n1, 2, 3] - tmp / 2.0
    dinp[n0 + 1, n1, 2, 3] = dinp[n0 + 1, n1, 2, 3] - tmp / 2.0
    dinp[n0, n1 + 1, 2, 3] = dinp[n0, n1 + 1, 2, 3] + tmp / 2.0
    dinp[n0 + 1, n1 + 1, 2, 3] = dinp[n0 + 1, n1 + 1, 2, 3] + tmp / 2.0

    tmp = i1 * inner(i0, i1 - 1, j0, j1) # df0f1
    dinp[n0, n1, 3, 1] = dinp[n0, n1, 3, 1] - tmp / 2.0
    dinp[n0 + 1, n1, 3, 1] = dinp[n0 + 1, n1, 3, 1] + tmp / 2.0
    dinp[n0, n1 + 1, 3, 1] = dinp[n0, n1 + 1, 3, 1] - tmp / 2.0
    dinp[n0 + 1, n1 + 1, 3, 1] = dinp[n0 + 1, n1 + 1, 3, 1] + tmp / 2.0

    tmp = j1 * inner(i0, i1, j0, j1 - 1) # df0f1
    dinp[n0, n1, 3, 2] = dinp[n0, n1, 3, 2] - tmp / 2.0
    dinp[n0 + 1, n1, 3, 2] = dinp[n0 + 1, n1, 3, 2] + tmp / 2.0
    dinp[n0, n1 + 1, 3, 2] = dinp[n0, n1 + 1, 3, 2] - tmp / 2.0
    dinp[n0 + 1, n1 + 1, 3, 2] = dinp[n0 + 1, n1 + 1, 3, 2] + tmp / 2.0

    tmp = inner(i0, i1, j0, j1) / (zs[k0 + 1] - zs[k0]) # df0df1
    dinp[n0, n1, 3, 3] = dinp[n0, n1, 3, 3] + tmp
    dinp[n0 + 1, n1, 3, 3] = dinp[n0 + 1, n1, 3, 3] - tmp
    dinp[n0, n1 + 1, 3, 3] = dinp[n0, n1 + 1, 3, 3] - tmp
    dinp[n0 + 1, n1 + 1, 3, 3] = dinp[n0 + 1, n1 + 1, 3, 3] + tmp
    
    tmp = inner(i0, i1, j0, j1) * (zs[k0 + 1] - zs[k0]) / 3.0
    inp[n0, n1] = inp[n0, n1] + tmp
    inp[n0 + 1, n1] = inp[n0 + 1, n1] + tmp / 2.0
    inp[n0, n1 + 1] = inp[n0, n1 + 1] + tmp / 2.0
    inp[n0 + 1, n1 + 1] = inp[n0 + 1, n1 + 1] + tmp
  }
}

K = matrix(0, nrow = 3 * N, ncol = 3 * N)
M = matrix(0, nrow = 3 * N, ncol = 3 * N)
for(n0 in 1:N) {
  for(n1 in 1:N) {
    k0 = idxs[[n0]][3]
    k1 = idxs[[n1]][3]
    #ii = if(k0 < m) 1 else 2
    
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
# 3.734092 4.729633 4.734347 6.431031 6.615935 8.074131 8.519344 9.863590

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

