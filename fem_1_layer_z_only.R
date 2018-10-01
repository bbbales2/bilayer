library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)

IN = 20
JN = 20
KN = 2
X = 0.03
Y = 0.02
Z = 0.0001
dx = X / (IN - 1)
dy = Y / (JN - 1)
dz = Z / (KN - 1)
ds = c(dx, dy, dz)
densities = c(8700)
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

cs = list(Cvoigt(cm1))

idx = function(i, j, k) {
  (i - 1) * JN * KN + (j - 1) * KN + k
}

idxt = function(ijk) {
  idx(ijk[1], ijk[2], ijk[3])
}

idxs = list()
for(i in 1:(IN - 1)) {
  for(j in 1:(JN - 1)) {
    for(k in 1:(KN - 1)) {
      idxs[[idx(i, j, k)]] = c(i, j, k)
    }
  }
}

N = IN * JN * KN

offsets = list(c(0, 0, 0),
               c(0, 0, 1),
               c(0, 1, 0),
               c(0, 1, 1),
               c(1, 0, 0),
               c(1, 0, 1),
               c(1, 1, 0),
               c(1, 1, 1))

O = length(offsets)
dinp = array(0, c(O, O, 3, 3))
inp = array(0, c(O, O))

for(oi1 in 1:O) {
  for(oi2 in 1:O) {
    o1 = offsets[[oi1]]
    o2 = offsets[[oi2]]
    
    for(i in 1:3) {
      for(j in 1:3) {
        d = 1.0
        if(i == j) {
          for(k in 1:3) {
            if(i == k) {
              if(o1[k] == o2[k]) {
                d = d * 1.0
              } else {
                d = d * -1.0
              }
            } else {
              if(o1[k] == o2[k]) {
                d = d * 1.0 / 3.0
              } else {
                d = d * 1.0 / 6.0
              }
            }
          }
        } else {
          for(k in 1:3) {
            if(i == k) {
              if(o1[k] == 0) {
                d = d * -1.0 / 2.0
              } else {
                d = d * 1.0 / 2.0
              }
            } else if(j == k) {
              if(o2[k] == 0) {
                d = d * -1.0 / 2.0
              } else {
                d = d * 1.0 / 2.0
              }
            } else {
              if(o1[k] == o2[k]) {
                d = d * 1.0 / 3.0
              } else {
                d = d * 1.0 / 6.0
              }
            }
          }
        }
        
        dinp[oi1, oi2, i, j] = d * dx * dy * dz / (ds[i] * ds[j])
      }
    }
    
    m = 1.0
    if(o1[1] == o2[1]) {
      m = m * 1.0 / 3.0
    } else {
      m = m * 1.0 / 6.0
    }

    if(o1[2] == o2[2]) {
      m = m * 1.0 / 3.0
    } else {
      m = m * 1.0 / 6.0
    }

    if(o1[3] == o2[3]) {
      m = m * 1.0 / 3.0
    } else {
      m = m * 1.0 / 6.0
    }

    inp[oi1, oi2] = m * dx * dy * dz
  }
}

#K = matrix(0, nrow = N, ncol = N)
#M = matrix(0, nrow = N, ncol = N)
K = matrix(0, nrow = 3 * N, ncol = 3 * N)
M = matrix(0, nrow = 3 * N, ncol = 3 * N)

for(n in 1:length(idxs)) {
  for(oi1 in 1:length(offsets)) {
    n1 = idxt(idxs[[n]] + offsets[[oi1]])
    for(oi2 in 1:length(offsets)) {
      n2 = idxt(idxs[[n]] + offsets[[oi2]])
      
      #K[n1, n2] = K[n1, n2] + sum(cs[[1]][3,, 3, ] * dinp[oi1, oi2,,])
      #M[n1, n2] = M[n1, n2] + densities[1] * inp[oi1, oi2]
      
      for(i in 1:3) {
        for(k in 1:3) {
          K[3 * (n1 - 1) + i, 3 * (n2 - 1) + k] = K[3 * (n1 - 1) + i, 3 * (n2 - 1) + k] + sum(cs[[1]][i,, k, ] * dinp[oi1, oi2,,])
        }
        M[3 * (n1 - 1) + i, 3 * (n2 - 1) + i] = M[3 * (n1 - 1) + i, 3 * (n2 - 1) + i] + densities[1] * inp[oi1, oi2]
      }
    }
  }
}

isSymmetric(K)
r = geigen(K, M, TRUE)

write.table(M, "m.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(K, "k.csv", sep = ",", row.names = FALSE, col.names = FALSE)

print(r$values[1:25])
print(1e-3 * sqrt(r$values[7:14] * 1e9) / (pi * 2))
# basis functions: [1] 38714.62 50746.73 51924.30 63328.94 76872.74 82737.92 87248.11 90207.19
# 12, 12, 5: 67085.41 100628.12 120939.95 135540.36 168811.16 203310.53 206751.13 214092.56
# 10, 10, 5: 67197.9 100796.8 121142.7 136445.3 169638.8 204667.9 209814.9 215417.1
# 5, 5, 5:   68586.85 102880.27 123646.70 147442.91 179788.10 221164.36 231555.24 239600.60

{
  u = matrix(0, nrow = IN, ncol = JN)
  for(j in 1:JN) {
    for(i in 1:IN) {
      u[i, j] = r$vectors[3 * (idx(i, j, 1) - 1) + 1, 7]
    }
  }
  levelplot(u)
}

#[1] -8.030616e-11  0.000000e+00  0.000000e+00  0.000000e+00  1.916770e-10  7.877108e-10  5.917111e+01  1.016660e+02  1.064391e+02  1.583304e+02  2.332945e+02  2.702520e+02  3.005189e+02
#[14]  3.212492e+02
#> print(sqrt(r$values[7:14] * 1e9) / (pi * 2))
#[1] 38714.62 50746.73 51924.30 63328.94 76872.74 82737.92 87248.11 90207.19
