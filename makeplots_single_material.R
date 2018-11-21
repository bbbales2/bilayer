#!/usr/bin/env Rscript
library(argparse)
library(tidyverse)
library(ggplot2)
library(rstan)
library(geigen)

parser = ArgumentParser()

parser$add_argument("--data", default="", help = "File to be sourced that has all config stuff in it (or an output file if --refine enabled)")
parser$add_argument("--output", default="", help = "Place to store the output plots")
parser$add_argument("--noscale", action="store_true", default=FALSE, help = "Do not scale plot size to geometry")

args = parser$parse_args()

if(length(args$config) == 0) {
  print("No configuration file provided")
  quit()
}

dat = read_rdump("fwd.dat")

not_defined = setdiff(c("P", "N", "X", "Y", "Z", "density", "c"), names(dat))
if(length(not_defined) > 0) {
  print("Configuration file does not include:", paste(not_defined, collapse = ", "), ". These are all required")
  quit()
}

IN = dat$P
JN = dat$P
KN = dat$P

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

evals = r$values * (r$values > 0)
#for (w in 1:30) {
w = 1
R = 21
#w = 15
d = 1
xs = seq(0.0, X, length = R)
ys = seq(0.0, Y, length = R)
zs = seq(0.0, Z, length = R)

approx = function(N, d, w, x, y, z) {
  inner = function(n) {
    i = idxs[[n]][1]
    j = idxs[[n]][2]
    k = idxs[[n]][3]
    r$vectors[3 * (n - 1) + d, w] * (x^i) * (y^j) * (z^k)
  }
  sapply(1:N, inner) %>% rowSums
}

tmp = expand.grid(x = xs, y = ys, z = zs) %>% as.tibble %>% mutate(ux = approx(N, 1, w, x, y, z),
                                                                   uy = approx(N, 2, w, x, y, z),
                                                                   uz = approx(N, 3, w, x, y, z)) %>%
  mutate(ux = x + (X / 4.0) * ux / max(abs(ux)),
         uy = y + (Y / 4.0) * uy / max(abs(uy)),
         uz = z + (Z / 4.0) * uz / max(abs(uz)))

points3D(tmp$ux, tmp$uy, tmp$uz, colvar = NULL, pch = '.')
plotdev(phi = 0, theta = 0)
points3D(tmp$ux, tmp$uy, tmp$uz, colvar = NULL, pch = '.')
plotdev(phi = 0, theta = 90)
points3D(tmp$ux, tmp$uy, tmp$uz, colvar = NULL, pch = '.')
plotdev(phi = -90, theta = 0)
points3D(tmp$ux, tmp$uy, tmp$uz, colvar = NULL, pch = '.')
plotdev(phi = -90, theta = 0)
par(mar =rep(1.0, 4))

dtostring = c("x", "y", "z")

#tplot = 
  bind_rows(expand.grid(c1 = ys, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, xs[1], c1, c2), direction = 'along x\nc1 = y\nc2 = z', position = 0.0), 
          expand.grid(c1 = ys, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, xs[(R - 1) / 2], c1, c2), direction = 'along x\nc1 = y\nc2 = z', position = 0.5),
          expand.grid(c1 = ys, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, xs[R], c1, c2), direction = 'along x\nc1 = y\nc2 = z', position = 1.0),
          expand.grid(c1 = xs, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, ys[1], c2), direction = 'along y\nc1 = x\nc2 = z', position = 0.0), 
          expand.grid(c1 = xs, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, ys[(R - 1) / 2], c2), direction = 'along y\nc1 = x\nc2 = z', position = 0.5),
          expand.grid(c1 = xs, c2 = zs) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, ys[R], c2), direction = 'along y\nc1 = x\nc2 = z', position = 1.0),
          expand.grid(c1 = xs, c2 = ys) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, c2, zs[1]), direction = 'along z\nc1 = x\nc2 = y', position = 0.0), 
          expand.grid(c1 = xs, c2 = ys) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, c2, zs[(R - 1) / 2]), direction = 'along z\nc1 = x\nc2 = y', position = 0.5),
          expand.grid(c1 = xs, c2 = ys) %>% as.tibble %>% mutate(u = approx(N, d, w, c1, c2, zs[R]), direction = 'along z\nc1 = x\nc2 = y', position = 1.0)) %>%
  ggplot() +
  geom_raster(aes(c1, c2, fill = u)) +
  scale_fill_gradient2(low = "#0000FF", high = "#FF0000", mid = "white", midpoint = 0.0) +
  facet_grid(direction ~ position, space = "free", scales = "free") +
  ggtitle(sprintf("%s displacement, mode %d: %.3fKhz", dtostring[d], w - 6, sqrt(evals * 1e9) / (1e3 * pi * 2)))

ggsave(filename = sprintf("ux_mode_%03d.png", w - 6), width = 10, height = 8, plot = tplot)
#}
