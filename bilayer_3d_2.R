library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 8
JN = 8
KN = 22
X = 0.02
Y = 0.03
Z = 0.011
zs = seq(0.0, Z, length = KN)
densities = c(8700, 8000)
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
cm2 = buildcm(269.231, 115.385, 76.923)
B = 0.01
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

cs = list(Cvoigt(cm1), Cvoigt(cm2))

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

dinp = array(0, c(2, N, N, 3, 3))
inp = array(0, c(2, N, N))

for(ii in 1:length(densities)) {
  for(n0 in 1:N) {
    for(n1 in 1:N) {
      i0 = idxs[[n0]][1]
      j0 = idxs[[n0]][2]
      k0 = idxs[[n0]][3]
      i1 = idxs[[n1]][1]
      j1 = idxs[[n1]][2]
      k1 = idxs[[n1]][3]
      
      if(ii == 1) {
        if(k0 >= m) {
          next
        }
      } else {
        if(k0 < m) {
          next
        }
      }
      
      if(k0 != k1) {
        next
      }
      
      if(k0 == KN) {
        next
      }
  
      tmp = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
      dinp[ii, n0, n1, 1, 1] = dinp[ii, n0, n1, 1, 1] + tmp
      dinp[ii, n0 + 1, n1, 1, 1] = dinp[ii, n0 + 1, n1, 1, 1] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 1, 1] = dinp[ii, n0, n1 + 1, 1, 1] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 1, 1] = dinp[ii, n0 + 1, n1 + 1, 1, 1] + tmp
  
      tmp = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
      dinp[ii, n0, n1, 1, 2] = dinp[ii, n0, n1, 1, 2] + tmp
      dinp[ii, n0 + 1, n1, 1, 2] = dinp[ii, n0 + 1, n1, 1, 2] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 1, 2] = dinp[ii, n0, n1 + 1, 1, 2] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 1, 2] = dinp[ii, n0 + 1, n1 + 1, 1, 2] + tmp
  
      tmp = i0 * inner(i0 - 1, i1, j0, j1) # f0df1
      dinp[ii, n0, n1, 1, 3] = dinp[ii, n0, n1, 1, 3] - tmp / 2.0
      dinp[ii, n0 + 1, n1, 1, 3] = dinp[ii, n0 + 1, n1, 1, 3] - tmp / 2.0
      dinp[ii, n0, n1 + 1, 1, 3] = dinp[ii, n0, n1 + 1, 1, 3] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 1, 3] = dinp[ii, n0 + 1, n1 + 1, 1, 3] + tmp / 2.0
  
      tmp = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
      dinp[ii, n0, n1, 2, 1] = dinp[ii, n0, n1, 2, 1] + tmp
      dinp[ii, n0 + 1, n1, 2, 1] = dinp[ii, n0 + 1, n1, 2, 1] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 2, 1] = dinp[ii, n0, n1 + 1, 2, 1] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 2, 1] = dinp[ii, n0 + 1, n1 + 1, 2, 1] + tmp
  
      tmp = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1) * (zs[k0 + 1] - zs[k0]) / 3.0 # f0f1
      dinp[ii, n0, n1, 2, 2] = dinp[ii, n0, n1, 2, 2] + tmp
      dinp[ii, n0 + 1, n1, 2, 2] = dinp[ii, n0 + 1, n1, 2, 2] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 2, 2] = dinp[ii, n0, n1 + 1, 2, 2] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 2, 2] = dinp[ii, n0 + 1, n1 + 1, 2, 2] + tmp
  
      tmp = j0 * inner(i0, i1, j0 - 1, j1) # f0df1
      dinp[ii, n0, n1, 2, 3] = dinp[ii, n0, n1, 2, 3] - tmp / 2.0
      dinp[ii, n0 + 1, n1, 2, 3] = dinp[ii, n0 + 1, n1, 2, 3] - tmp / 2.0
      dinp[ii, n0, n1 + 1, 2, 3] = dinp[ii, n0, n1 + 1, 2, 3] + tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 2, 3] = dinp[ii, n0 + 1, n1 + 1, 2, 3] + tmp / 2.0
  
      tmp = i1 * inner(i0, i1 - 1, j0, j1) # df0f1
      dinp[ii, n0, n1, 3, 1] = dinp[ii, n0, n1, 3, 1] - tmp / 2.0
      dinp[ii, n0 + 1, n1, 3, 1] = dinp[ii, n0 + 1, n1, 3, 1] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 3, 1] = dinp[ii, n0, n1 + 1, 3, 1] - tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 3, 1] = dinp[ii, n0 + 1, n1 + 1, 3, 1] + tmp / 2.0
  
      tmp = j1 * inner(i0, i1, j0, j1 - 1) # df0f1
      dinp[ii, n0, n1, 3, 2] = dinp[ii, n0, n1, 3, 2] - tmp / 2.0
      dinp[ii, n0 + 1, n1, 3, 2] = dinp[ii, n0 + 1, n1, 3, 2] + tmp / 2.0
      dinp[ii, n0, n1 + 1, 3, 2] = dinp[ii, n0, n1 + 1, 3, 2] - tmp / 2.0
      dinp[ii, n0 + 1, n1 + 1, 3, 2] = dinp[ii, n0 + 1, n1 + 1, 3, 2] + tmp / 2.0
  
      tmp = inner(i0, i1, j0, j1) / (zs[k0 + 1] - zs[k0]) # df0df1
      dinp[ii, n0, n1, 3, 3] = dinp[ii, n0, n1, 3, 3] + tmp
      dinp[ii, n0 + 1, n1, 3, 3] = dinp[ii, n0 + 1, n1, 3, 3] - tmp
      dinp[ii, n0, n1 + 1, 3, 3] = dinp[ii, n0, n1 + 1, 3, 3] - tmp
      dinp[ii, n0 + 1, n1 + 1, 3, 3] = dinp[ii, n0 + 1, n1 + 1, 3, 3] + tmp
      
      tmp = inner(i0, i1, j0, j1) * (zs[k0 + 1] - zs[k0]) / 3.0
      inp[ii, n0, n1] = inp[ii, n0, n1] + tmp
      inp[ii, n0 + 1, n1] = inp[ii, n0 + 1, n1] + tmp / 2.0
      inp[ii, n0, n1 + 1] = inp[ii, n0, n1 + 1] + tmp / 2.0
      inp[ii, n0 + 1, n1 + 1] = inp[ii, n0 + 1, n1 + 1] + tmp
    }
  }
}

K = matrix(0, nrow = 3 * N, ncol = 3 * N)
M = matrix(0, nrow = 3 * N, ncol = 3 * N)
for(ii in 1:length(densities)) {
  for(n0 in 1:N) {
    cat(n0, "\n");
    for(n1 in 1:N) {
      k0 = idxs[[n0]][3]
      k1 = idxs[[n1]][3]
  
      for(i in 1:3) {
        for(k in 1:3) {
          total = 0.0
  
          for(j in 1:3) {
            for(l in 1:3) {
              total = total + cs[[ii]][i, j, k, l] * dinp[ii, n0, n1, j, l]
            }
          }
  
          K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] = K[3 * (n0 - 1) + i, 3 * (n1 - 1) + k] + total
        }
        M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] + densities[ii] * inp[ii, n0, n1]
      }
    }
  }
}

r = geigen(K, M, TRUE)

print(r$values[1:25])
print(sqrt(r$values[7:14] * 1e9) / (pi * 2))
# 3.734092e+00  4.729633e+00  4.734347e+00  6.431031e+00  6.615935e+00  8.074131e+00
# 3.024603e+00  3.202178e+00  3.542420e+00  5.058626e+00  5.895012e+00  6.216565e+00
# 4.251224e+00  4.271725e+00  5.319324e+00  5.357855e+00  6.569780e+00  7.108487e+00

# 3.692754e+00  4.687945e+00  4.724703e+00  6.199050e+00  6.611294e+00  7.909246e+00
# 2.990619e+00  3.185342e+00  3.509285e+00  4.880287e+00  5.881217e+00  6.094774e+00
# 4.206195e+00  4.308935e+00  5.259553e+00  5.323833e+00  6.564090e+00  7.043638e+00

# 4.894889e+00  5.301620e+00  5.579165e+00  6.930004e+00  8.488756e+00  9.075953e+00

# 4.850074e+00  5.304551e+00  5.724816e+00  7.078923e+00  8.783755e+00  9.096595e+00


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

