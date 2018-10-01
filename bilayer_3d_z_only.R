library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 10
JN = 10
KN = 5
X = 0.02
Y = 0.03
Z = 0.0001
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

cm1 = buildcm(250, 150, 140)
cm2 = buildcm(250, 150, 140)
#cm2 = buildcm(269.231, 115.385, 76.923)
B = Z / 2.0

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

cs = list(Cvoigt(cm1), Cvoigt(cm1))

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

dinp = array(0, c(2, N, N, 3, 3))
inp = array(0, c(2, N, N))

for(ii in 1:length(densities)) {
  a = if(ii == 1) 0 else B
  b = if(ii == 1) B else Z
  for(n0 in 1:N) {
    for(n1 in 1:N) {
      i0 = idxs[[n0]][1]
      j0 = idxs[[n0]][2]
      k0 = idxs[[n0]][3]
      i1 = idxs[[n1]][1]
      j1 = idxs[[n1]][2]
      k1 = idxs[[n1]][3]
      
      dinp[ii, n0, n1, 1, 1] = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1, k0, k1, a, b)
      dinp[ii, n0, n1, 1, 2] = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1, k0, k1, a, b)
      dinp[ii, n0, n1, 1, 3] = i0 * k1 * inner(i0 - 1, i1, j0, j1, k0, k1 - 1, a, b)
      dinp[ii, n0, n1, 2, 1] = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1, k0, k1, a, b)
      dinp[ii, n0, n1, 2, 2] = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1, k0, k1, a, b)
      dinp[ii, n0, n1, 2, 3] = j0 * k1 * inner(i0, i1, j0 - 1, j1, k0, k1 - 1, a, b)
      dinp[ii, n0, n1, 3, 1] = k0 * i1 * inner(i0, i1 - 1, j0, j1, k0 - 1, k1, a, b)
      dinp[ii, n0, n1, 3, 2] = k0 * j1 * inner(i0, i1, j0, j1 - 1, k0 - 1, k1, a, b)
      dinp[ii, n0, n1, 3, 3] = k0 * k1 * inner(i0, i1, j0, j1, k0 - 1, k1 - 1, a, b)
      inp[ii, n0, n1] = inner(i0, i1, j0, j1, k0, k1, a, b)
    }
  }
}

#K = matrix(0, nrow = 3 * N, ncol = 3 * N)
#M = matrix(0, nrow = 3 * N, ncol = 3 * N)
K = matrix(0, nrow = N, ncol = N)
M = matrix(0, nrow = N, ncol = N)
for(ii in 1:length(densities)) {
  for(n0 in 1:N) {
    for(n1 in 1:N) {
      # for(i in 1:3) {
      #   for(k in 1:3) {
      #     total = 0.0
      #     
      #     for(j in 1:3) {
      #       for(l in 1:3) {
      #         total = total + cs[[ii]][i, j, k, l] * dinp[ii, n0, n1, j, l]
      #       }
      #     }
      #     
      #     K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] = K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] + total
      #   }
      #   M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] + densities[ii] * inp[ii, n0, n1]
      # }
      
      K[n0, n1] = K[n0, n1] + sum(cs[[ii]][3,, 3, ] * dinp[ii, n0, n1,,])
      M[n0, n1] = M[n0, n1] + densities[ii] * inp[ii, n0, n1]
    }
  }
}

r = geigen(K, M, TRUE)

print(r$values[1:14])
print(1e-3 * sqrt(r$values[2:14] * 1e9) / (pi * 2))
# 40623.81 48710.37 52785.86 64962.03 80016.28 85147.94 88924.69 89817.98
# 38714.62 50746.73 51924.30 63328.94 76872.74 82737.92 87248.11 90207.19
{
  xs = seq(0.0, X, length = 20)
  ys = seq(0.0, Y, length = 20)
  z = 0.5
  u = matrix(0, nrow = length(xs), ncol = length(ys))
  for(n in 1:10) {
    i = idxs[[n]][1]
    j = idxs[[n]][2]
    k = idxs[[n]][3]
    u = u + outer(xs, ys, function(x, y) { r$vectors[3 * (n - 1) + 2, 1] * (x^i) * (y^j) * (z^k) })
  }
  levelplot(u)
}

# 7.442681e+00  8.283205e+00  8.283205e+00  8.613477e+00  8.750229e+00  8.981048e+00  9.263764e+00 9.934843e+00  9.934843e+00  1.038664e+01  1.090038e+01  1.166401e+01  1.222080e+01  1.235710e+01

