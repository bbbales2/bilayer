library(lattice)
library(tidyverse)
library(ggplot2)
library(geigen)
library(Ryacas)

IN = 7
JN = 6
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

cm2 = buildcm(250, 150, 140)
cm1 = buildcm(150, 100, 40)

cs = list(Cvoigt(cm1), Cvoigt(cm2))

x = Sym("x")
fs = list("1 - 3 * x + 2 * x * x", "4 * x - 4 * x * x", "-x + 2 * x * x")#list("x", "1 - x")
bp = array(0, c(length(fs), 2, length(fs), 2))
for(i in 1:length(fs)) {
  for(di in 1:2) {
    for(j in 1:length(fs)) {
      for(dj in 1:2) {
        f = Integrate(paste0('(', ifelse(di == 2, "D(x) ", ""), fs[[i]], ') * (', ifelse(dj == 2, "D(x) ", ""), fs[[j]] ,')'), x) %>% as.expression
        bp[i, di, j, dj] = eval(f, list(x = 1.0)) - eval(f, list(x = 0.0))
      }
    }
  }
}

inner = function(i0, i1, j0, j1) {
  # int_0_X x^i0 * x^i1 dx
  ret = 0.0
  if(min(i0, i1, j0, j1) >= 0) {
    ret = (X^(i0 + i1 + 1) / (i0 + i1 + 1)) *
      (Y^(j0 + j1 + 1) / (j0 + j1 + 1))
  }
  ret
}

ntoijk = list()
ntom = list()
nlength = list()
mton = list()
n = 1
m = 1
for(i in 0:IN) {
  for(j in 0:JN) {
    for(k in 1:KN) {
      if((i + j) <= max(IN, JN)) {
        ntoijk[[n]] = c(i, j, k)
        ntom[[n]] = m
        if(k < KN) {
          nlength[[n]] = (length(fs) - 1)
        } else {
          nlength[[n]] = 1
        }
        for(l in 1:nlength[[n]]) {
          mton[[m + l - 1]] = n 
        }
        m = m + nlength[[n]]
        n = n + 1
      }
    }
  }
}

N = length(ntoijk)
M = nlength %>% unlist %>% sum

dinp = array(0, c(2, M, M, 3, 3))
inp = array(0, c(2, M, M))

for(ii in 1:length(densities)) {
  for(n0 in 1:N) {
    for(n1 in 1:N) {
      i0 = ntoijk[[n0]][1]
      j0 = ntoijk[[n0]][2]
      k0 = ntoijk[[n0]][3]
      i1 = ntoijk[[n1]][1]
      j1 = ntoijk[[n1]][2]
      k1 = ntoijk[[n1]][3]
  
      if(ii == 1) {
        if(k0 >= zint) {
          next
        }
      } else {
        if(k0 < zint) {
          next
        }
      }
  
      if(k0 != k1) {
        next
      }
  
      if(k0 == KN) {
        next
      }
  
      dz = (zs[k0 + 1] - zs[k0])
      r0 = ntom[[n0]]:(ntom[[n0]] + nlength[[n0]])
      r1 = ntom[[n1]]:(ntom[[n1]] + nlength[[n1]])
  
      tmp = i0 * i1 * inner(i0 - 1, i1 - 1, j0, j1) # f0f1
      dinp[ii, r0, r1, 1, 1] = dinp[ii, r0, r1, 1, 1] + tmp * bp[, 1,, 1] * dz
  
      tmp = i0 * j1 * inner(i0 - 1, i1, j0, j1 - 1) # f0f1
      dinp[ii, r0, r1, 1, 2] = dinp[ii, r0, r1, 1, 2] + tmp * bp[, 1,, 1] * dz
  
      tmp = i0 * inner(i0 - 1, i1, j0, j1) # f0df1
      dinp[ii, r0, r1, 1, 3] = dinp[ii, r0, r1, 1, 3] + tmp * bp[, 1,, 2]
  
      tmp = j0 * i1 * inner(i0, i1 - 1, j0 - 1, j1) # f0f1
      dinp[ii, r0, r1, 2, 1] = dinp[ii, r0, r1, 2, 1] + tmp * bp[, 1,, 1] * dz
  
      tmp = j0 * j1 * inner(i0, i1, j0 - 1, j1 - 1) # f0f1
      dinp[ii, r0, r1, 2, 2] = dinp[ii, r0, r1, 2, 2] + tmp * bp[, 1,, 1] * dz
  
      tmp = j0 * inner(i0, i1, j0 - 1, j1) # f0df1
      dinp[ii, r0, r1, 2, 3] = dinp[ii, r0, r1, 2, 3] + tmp * bp[, 1,, 2]
  
      tmp = i1 * inner(i0, i1 - 1, j0, j1) # df0f1
      dinp[ii, r0, r1, 3, 1] = dinp[ii, r0, r1, 3, 1] + tmp * bp[, 2,, 1]
  
      tmp = j1 * inner(i0, i1, j0, j1 - 1) # df0f1
      dinp[ii, r0, r1, 3, 2] = dinp[ii, r0, r1, 3, 2] + tmp * bp[, 2,, 1]
  
      tmp = inner(i0, i1, j0, j1) # df0df1
      dinp[ii, r0, r1, 3, 3] = dinp[ii, r0, r1, 3, 3] + tmp * bp[, 2,, 2] / dz
  
      tmp = inner(i0, i1, j0, j1)
      inp[ii, r0, r1] = inp[ii, r0, r1] + tmp * bp[, 1,, 1] * dz
    }
  }
}

K = matrix(0, nrow = 3 * M, ncol = 3 * M)
W = matrix(0, nrow = 3 * M, ncol = 3 * M)
for(ii in 1:length(densities)) {
  for(m0 in 1:M) {
    for(m1 in 1:M) {
      k0 = ntoijk[[mton[[m0]]]][[3]]
      k1 = ntoijk[[mton[[m0]]]][[3]]
      
      if(k0 != k1) {
        next
      }
      
      for(i in 1:3) {
        for(k in 1:3) {
          total = 0.0
          
          for(j in 1:3) {
            for(l in 1:3) {
              total = total + cs[[ii]][i, j, k, l] * dinp[ii, m0, m1, j, l]
            }
          }
          
          K[3 * (m0 - 1) + i, 3 * (m1 - 1) + k] = K[3 * (m0 - 1) + i, 3 * (m1 - 1) + k] + total
        }
        W[3 * (m0 - 1) + i, 3 * (m1 - 1) + i] = W[3 * (m0 - 1) + i, 3 * (m1 - 1) + i] + densities[ii] * inp[ii, m0, m1]
      }
    }
  }
}

r = geigen(K, W, TRUE)

print(r$values[1:25])
print(sqrt(r$values[7:14] * 1e9) / (pi * 2))

#11454.55 13885.22 21559.08 21846.67 25507.50 26803.93 27196.34 32373.27
#11414.40 13781.38 21451.18 21753.09 25504.56 26562.91 26832.72 31874.41
#11408.66 13747.66 21450.40 21721.79 25502.83 26426.39 26725.75 31712.02
