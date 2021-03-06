library(Rvcg)
library(argparse)
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)
library(rgl)
library(tidyverse)
library(ggplot2)
library(rstan)
library(geigen)

dat = read_rdump("fwd.dat")
scale = 10.0

not_defined = setdiff(c("P", "N", "X", "Y", "Z", "density", "c"), names(dat))
if(length(not_defined) > 0) {
  print("Configuration file does not include:", paste(not_defined, collapse = ", "), ". These are all required")
  quit()
}

IN = dat$P
JN = dat$P
KN = dat$P

X = dat$X
Y = dat$Y
Z = dat$Z

cm = dat$c

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

cs = Cvoigt(cm)

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
      M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] = M[3 * (n0 - 1) + i, 3 * (n1 - 1) + i] + dat$density * inp[n0, n1]
    }
  }
}

r = geigen(K, M, TRUE)
print(1e-3 * sqrt(r$values[7:16] * 1e9) / (pi * 2))

frequencies = (1e-3 * sqrt(r$values * (r$values > 0) * 1e9) / (pi * 2))[7:26]
points = list()
clear3d()

for(mode in 1:20) {
  w = mode + 6
  R = 21
  d = 1
  xs = seq(0.0, X, length = R * X / min(X, Y, Z))
  ys = seq(0.0, Y, length = R * Y / min(X, Y, Z))
  zs = seq(0.0, Z, length = R * Z / min(X, Y, Z))
  
  approx = function(N, d, w, x, y, z) {
    inner = function(n) {
      i = idxs[[n]][1]
      j = idxs[[n]][2]
      k = idxs[[n]][3]
      r$vectors[3 * (n - 1) + d, w] * (x^i) * (y^j) * (z^k)
    }
    sapply(1:N, inner) %>% rowSums
  }
  
  tmp = bind_rows(expand.grid(x = xs[1], y = ys, z = zs) %>% as.tibble,
                  expand.grid(x = xs[length(xs)], y = ys, z = zs) %>% as.tibble,
                  expand.grid(x = xs, y = ys[1], z = zs) %>% as.tibble,
                  expand.grid(x = xs, y = ys[length(ys)], z = zs) %>% as.tibble,
                  expand.grid(x = xs, y = ys, z = zs[1]) %>% as.tibble,
                  expand.grid(x = xs, y = ys, z = zs[length(zs)]) %>% as.tibble) %>%
    mutate(ux = approx(N, 1, w, x, y, z),
           uy = approx(N, 2, w, x, y, z),
           uz = approx(N, 3, w, x, y, z)) %>%
    mutate(ux = x + (max(X, Y, Z) / scale) * ux / max(abs(ux), abs(uy), abs(uz)),
           uy = y + (max(X, Y, Z) / scale) * uy / max(abs(ux), abs(uy), abs(uz)),
           uz = z + (max(X, Y, Z) / scale) * uz / max(abs(ux), abs(uy), abs(uz)))
  
  mesh = vcgBallPivoting(cbind(tmp$ux, tmp$uy, tmp$uz), min(X, Y, Z) / 10.0)
  
  points[[mode]] = shade3d(mesh, col = "grey", specular = "black", aspect = c(X, Y, Z))
  
  print(paste0("Mode ", mode, ", ", sprintf("%0.3fKhz", frequencies[mode]), " (of 20)"))
}

aspect3d(X, Y, Z)
axesid = decorate3d()

rglwidget() %>%
  playwidget(start = 0, stop = length(points) - 1, interval = 1,
             subsetControl(value = 1, subsets = points))
