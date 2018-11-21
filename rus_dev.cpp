#include <Rcpp.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Core>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/Eigen/Sparse>
#include </home/bbales2/R/x86_64-pc-linux-gnu-library/3.4/RSpectra/include/Spectra/SymEigsSolver.h>
#include </home/bbales2/R/x86_64-pc-linux-gnu-library/3.4/RSpectra/include/Spectra/SymEigsShiftSolver.h>
#include </home/bbales2/R/x86_64-pc-linux-gnu-library/3.4/RSpectra/include/Spectra/MatOp/SparseSymShiftSolve.h>
using namespace Rcpp;

/*
 * Input:
 * X: length of sample in X direction
 * Y: length of sample in Y direction
 * Z: vector of discretization points of sample in Z direction
 * NX: Number of basis functions to use in X dimension
 * NY: Number of basis functions to use in Y dimension
 * Cs: Elastic constants of substrate
 * Cc: Elastic constants of coating
 * 
 * Output:
 * dKdci
 * M
 */
Eigen::Tensor<double, 4> Cvoigt(const Eigen::MatrixXd &cm) {
  Eigen::Tensor<double, 4> C(3, 3, 3, 3);
  
  std::vector<std::vector<std::array<int, 2> > > voigt = {
    { { 1, 1 } },
    { { 2, 2 } },
    { { 3, 3 } },
    { { 2, 3 }, { 3, 2 } },
    { { 1, 3 }, { 3, 1 } },
    { { 1, 2 }, { 2, 1 } }
  };
  
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      for(auto&& vi : voigt[i]) {
        int k = vi[0] - 1;
        int l = vi[1] - 1;
        for(auto&& vj : voigt[j]) {
          int n = vj[0] - 1;
          int m = vj[1] - 1;
          C(k, l, n, m) = cm(i, j);
        }
      }
    }
  }
  
  return C;
}

static double bp_[] = { 0.1333333333333333, 0.06666666666666666, -0.03333333333333333,
                        -0.500000000000000, 0.66666666666666666, -0.16666666666666666,
                        0.06666666666666666, 0.53333333333333333, 0.066666666666666666,
                        -0.66666666666666666, 0.000000000000000, 0.66666666666666666,
                        -0.0333333333333333, 0.066666666666666666, 0.133333333333333333,
                        0.16666666666666666, -0.66666666666666666, 0.500000000000000000,
                        -0.5000000000000000, -0.66666666666666666, 0.16666666666666666,
                        2.33333333333333333, -2.66666666666666666, 0.33333333333333333,
                        0.66666666666666666, 0.000000000000000000, -0.66666666666666666,
                        -2.66666666666666666, 5.3333333333333333, -2.66666666666666666,
                        -0.16666666666666666, 0.66666666666666666, 0.500000000000000000,
                        0.33333333333333333, -2.66666666666666666, 2.33333333333333333 };

double inner(double X, double Y, int i0, int i1, int j0, int j1) {
  double ret = 0.0;

  if (std::min(i0, std::min(i1, std::min(j0, j1))) >= 0) {
    ret = (std::pow(X, i0 + i1 + 1) / (i0 + i1 + 1)) *
      (std::pow(Y, j0 + j1 + 1) / (j0 + j1 + 1));
  }

  return ret;
}

// [[Rcpp::export]]
int calc_M(int IN, int JN, int KN) {
  int m = 0;
  for(int i = 0; i <= IN; i++) {
    for(int j = 0; j <= JN; j++) {
      for(int k = 0; k < KN; k++) {
        if((i + j) <= std::max(IN, JN)) {
          if(k < KN - 1) {
            m += 2;
          } else {
            m += 1;
          }
        }
      }
    }
  }
  return m;
}

// [[Rcpp::export]]
NumericVector build_dKdci_M(double X, double Y, NumericVector Z, int IN, int JN, int zint, NumericVector densities) {
  zint = zint - 1;
  
  Eigen::Tensor<double, 4> bp(3, 2, 3, 2);
  for(int i = 0; i < bp.size(); i++)
    bp(i) = bp_[i];
  
  std::vector<std::array<int, 3> > ntoijk;
  std::vector<int> ntom;
  std::vector<int> nlength;
  std::vector<int> mton;
  int n = 0;
  int m = 0;
  for(int i = 0; i <= IN; i++) {
    for(int j = 0; j <= JN; j++) {
      for(int k = 0; k < Z.size(); k++) {
        if((i + j) <= std::max(IN, JN)) {
          ntoijk.push_back({ i, j, k });
          //std::cout << i << ", " << j << ", " << k << std::endl;
          ntom.push_back(m);
          if(k < Z.size() - 1) {
            nlength.push_back(2);
          } else {
            nlength.push_back(1);
          }

          for(int l = 0; l < nlength.back(); l++) {
            mton.push_back(n);
          }

          m = m + nlength.back();
          n = n + 1;
        }
      }
    }
  }

  int N = ntoijk.size();
  int M = m;
  
  Eigen::Tensor<double, 5> dinp(2, M, M, 3, 3);
  dinp.setZero();
  Eigen::Tensor<double, 3> inp(2, M, M);
  inp.setZero();

  //std::cout << "N: " << N << std::endl;
  
  for(int ii = 0; ii < 2; ii++) {
    for(int n0 = 0; n0 < N; n0++) {
      for(int n1 = 0; n1 < N; n1++) {
        //std::cout << "hi" << std::endl;
        int i0 = ntoijk[n0][0];
        int j0 = ntoijk[n0][1];
        int k0 = ntoijk[n0][2];
        int i1 = ntoijk[n1][0];
        int j1 = ntoijk[n1][1];
        int k1 = ntoijk[n1][2];
        
        if(ii == 0) {
          if(k0 >= zint) {
            continue;
          }
        } else {
          if(k0 < zint) {
            continue;
          }
        }
        
        if(k0 != k1) {
          continue;
        }
        
        if(k0 == Z.size() - 1) {
          continue;
        }
        
        double dz = (Z[k0 + 1] - Z[k0]);
        //r0 = ntom[[n0]]:(ntom[[n0]] + nlength[[n0]]);
        //r1 = ntom[[n1]]:(ntom[[n1]] + nlength[[n1]]);
          
        for(int m0 = 0; m0 <= nlength[n0]; m0++) {
          for(int m1 = 0; m1 <= nlength[n1]; m1++) {
            int r0 = ntom[n0] + m0;
            int r1 = ntom[n1] + m1;
            
            //std::cout << "n0, n1: " << n0 << ", " << n1 << "; nlength[n0]: " << nlength[n0] << "; r0, r1: " << r0 << ", " << r1 << std::endl;

            double tmp = i0 * i1 * inner(X, Y, i0 - 1, i1 - 1, j0, j1); //# f0f1
            dinp(ii, r0, r1, 0, 0) = dinp(ii, r0, r1, 0, 0) + tmp * bp(m0, 0, m1, 0) * dz;
            
            tmp = i0 * j1 * inner(X, Y, i0 - 1, i1, j0, j1 - 1); //# f0f1
            dinp(ii, r0, r1, 0, 1) = dinp(ii, r0, r1, 0, 1) + tmp * bp(m0, 0, m1, 0) * dz;
            
            tmp = i0 * inner(X, Y, i0 - 1, i1, j0, j1); //# f0df1
            dinp(ii, r0, r1, 0, 2) = dinp(ii, r0, r1, 0, 2) + tmp * bp(m0, 0, m1, 1);
          
            tmp = j0 * i1 * inner(X, Y, i0, i1 - 1, j0 - 1, j1); //# f0f1
            dinp(ii, r0, r1, 1, 0) = dinp(ii, r0, r1, 1, 0) + tmp * bp(m0, 0, m1, 0) * dz;
            
            tmp = j0 * j1 * inner(X, Y, i0, i1, j0 - 1, j1 - 1); //# f0f1
            dinp(ii, r0, r1, 1, 1) = dinp(ii, r0, r1, 1, 1) + tmp * bp(m0, 0, m1, 0) * dz;
            
            tmp = j0 * inner(X, Y, i0, i1, j0 - 1, j1); //# f0df1
            dinp(ii, r0, r1, 1, 2) = dinp(ii, r0, r1, 1, 2) + tmp * bp(m0, 0, m1, 1);
          
            tmp = i1 * inner(X, Y, i0, i1 - 1, j0, j1); //# df0f1
            dinp(ii, r0, r1, 2, 0) = dinp(ii, r0, r1, 2, 0) + tmp * bp(m0, 1, m1, 0);
          
            tmp = j1 * inner(X, Y, i0, i1, j0, j1 - 1); //# df0f1
            dinp(ii, r0, r1, 2, 1) = dinp(ii, r0, r1, 2, 1) + tmp * bp(m0, 1, m1, 0);
          
            tmp = inner(X, Y, i0, i1, j0, j1); //# df0df1
            dinp(ii, r0, r1, 2, 2) = dinp(ii, r0, r1, 2, 2) + tmp * bp(m0, 1, m1, 1) / dz;
            
            tmp = inner(X, Y, i0, i1, j0, j1);
            //std::cout << "(" << m0 << " " << m1 << "): " << bp(m0, 0, m1, 0) << std::endl;
            inp(ii, r0, r1) = inp(ii, r0, r1) + tmp * bp(m0, 0, m1, 0) * dz;
            //std::cout << "(" << ii << " " << r0 << " " << r1 << "): " << tmp << ", " << bp(m0, 0, m1, 0) << ", " << dz << ", " << inp(ii, r0, r1) << std::endl;
          }
        }
      }
    }
  }
  
  /*NumericVector out(dinp.size() + inp.size());
  std::copy(dinp.data(), dinp.data() + dinp.size(), out.begin());
  std::copy(inp.data(), inp.data() + inp.size(), out.begin() + dinp.size());*/
  
  Eigen::Tensor<double, 5> dKdcij(3 * M, 3 * M, 6, 6, 2);
  dKdcij.setZero();
  Eigen::MatrixXd W(3 * M, 3 * M);
  W.setZero();
  
  for(int ii = 0; ii < 2; ii++) {
    for(int p = 0; p < 6; p++) {
      for(int q = 0; q < 6; q++) {
        Eigen::MatrixXd cm = Eigen::MatrixXd::Zero(6, 6);
        cm(p, q) = 1.0;
        Eigen::Tensor<double, 4> c = Cvoigt(cm);
        for(int m0 = 0; m0 < M; m0++) {
          for(int m1 = 0; m1 < M; m1++) {
            for(int i = 0; i < 3; i++) {
              for(int k = 0; k < 3; k++) {
                double total = 0.0;
                for(int j = 0; j < 3; j++) {
                  for(int l = 0; l < 3; l++) {
                    total += c(i, j, k, l) * dinp(ii, m0, m1, j, l);
                  }
                }
    
                dKdcij(3 * m0 + i, 3 * m1 + k, p, q, ii) = dKdcij(3 * m0 + i, 3 * m1 + k, p, q, ii) + total;
              }
            }
          }
        }
      }
    }
    
    for(int m0 = 0; m0 < M; m0++) {
      for(int m1 = 0; m1 < M; m1++) {
        for(int i = 0; i < 3; i++) {
          W(3 * m0 + i, 3 * m1 + i) = W(3 * m0 + i, 3 * m1 + i) + densities[ii] * inp(ii, m0, m1);
        }
      }
    }
  }
  
  Eigen::LLT<Eigen::MatrixXd> Wc = W.llt();
  
  Eigen::Tensor<double, 5> dKhatdcij(3 * M, 3 * M, 6, 6, 2);
  dKhatdcij.setZero();
  for(int ii = 0; ii < 2; ii++) {
    for(int p = 0; p < 6; p++) {
      for(int q = 0; q < 6; q++) {
        Eigen::Map<Eigen::MatrixXd> dKtmp(&dKdcij(0, 0, p, q, ii), 3 * M, 3 * M);
        Eigen::MatrixXd dKhattmp = Wc.matrixL().solve(Wc.matrixL().solve(dKtmp.transpose()).transpose());
        std::copy(dKhattmp.data(), dKhattmp.data() + 3 * M * 3 * M, &dKhatdcij(0, 0, p, q, ii));
      }
    }
  }
  
  double total = 0.0;
  for(int i = 0; i < dKhatdcij.size(); i++) {
    total += std::abs(dKhatdcij(i));
  }
  std::cout << "total: " << total << std::endl;
  
  NumericVector out(dKhatdcij.size());
  std::copy(dKhatdcij.data(), dKhatdcij.data() + dKhatdcij.size(), out.begin());
  return out;
}

// [[Rcpp::export]]
NumericVector calc_evals(NumericMatrix K) {
  Eigen::MatrixXd K_(K.rows(), K.cols());
  for(int j = 0; j < K.cols(); j++) {
    for(int i = 0; i < K.rows(); i++) {
      K_(i, j) = K(i, j);
    }
  }
  
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_);
  
  NumericVector out(es.eigenvalues().size());
  for(int i = 0; i < es.eigenvalues().size(); i++) {
    out(i) = es.eigenvalues()[i];
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector calc_evals_spectra(NumericMatrix K, int n, double ft) {
  double sigma = std::pow(ft * 2 * 3.14159265359, 2.0) / 1.0e9;
  //const Eigen::Map<Eigen::MatrixXd> K_(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(K));
  Eigen::MatrixXf K_(K.rows(), K.cols());
  for(int j = 0; j < K.cols(); j++) {
    for(int i = 0; i < K.rows(); i++) {
      K_(i, j) = K(i, j);
    }
  }
  
  Spectra::DenseSymShiftSolve<float> op(K_);
  Spectra::SymEigsShiftSolver<float, Spectra::LARGEST_MAGN, Spectra::DenseSymShiftSolve<float> > eigs(&op, n, 2 * n, sigma);

  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  // Retrieve results
  if(eigs.info() == Spectra::SUCCESSFUL) {
    NumericVector out(eigs.eigenvalues().size());
    for(int i = 0; i < eigs.eigenvalues().size(); i++) {
      out(i) = eigs.eigenvalues()[i];
    }
    
    return out;
  } else {
    return NumericVector();
  }
}

// [[Rcpp::export]]
NumericVector calc_evals_sparse(IntegerVector is, IntegerVector js, NumericVector vs, int N, int n, double ft) {
  double sigma = std::pow(ft * 2 * 3.14159265359, 2.0) / 1.0e9;
  
  int nnz = is.size();
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(nnz);
  
  for(int i = 0; i < is.size(); i++) {
    tripletList.push_back(T(is(i), js(i), vs(i)));
  }
  
  Eigen::SparseMatrix<double> mat(N, N);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  
  Spectra::SparseSymShiftSolve<double> op(mat);
  Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigs(&op, n, 2 * n, sigma);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  // Retrieve results
  if(eigs.info() == Spectra::SUCCESSFUL) {
    NumericVector out(eigs.eigenvalues().size());
    for(int i = 0; i < eigs.eigenvalues().size(); i++) {
      out(i) = eigs.eigenvalues()[i];
    }
    
    return out;
  } else {
    return NumericVector();
  }
}

