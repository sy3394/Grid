/*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/iterative/RefinedArnoldi.h

    Refined Arnoldi for interior eigenvalues of a general (non-Hermitian) operator.

    Single-vector Arnoldi builds an orthonormal Krylov basis  M V_m = V_{m+1} Hbar.
    Ritz VALUES are eigenvalues of the m x m Hessenberg; Ritz VECTORS are extracted by
    REFINED Ritz: for each value theta the vector minimises the Euclidean residual
    ||M u - theta u|| over span(V_m), i.e. u = V_m y with y the smallest right singular
    vector of (Hbar - theta [I;0]).  This avoids the spurious interior Ritz vectors of
    plain Rayleigh-Ritz.

    For interior eigenvalues of the Wilson Dirac operator, run on the shift-inverted
    operator  M = (D_W - sigma)^{-1}; a Ritz value mu maps back as lambda = sigma + 1/mu.
    getEvals() returns the mu (caller maps to lambda); convergence is judged on the raw
    operator dW via setRawCheck (residual ||D_W u - lambda u|| / ||u||).

    On the Wilson operator this single-vector + refined scheme is both simpler and more
    effective for interior modes than gamma5-block Lanczos (the gamma5 metric adds an
    oblique-projection penalty and the block seed halves the Krylov depth).
*************************************************************************************/
#pragma once

#include <functional>

NAMESPACE_BEGIN(Grid);

enum RefinedRitzSort {
  RASortAbsImagAscending,  // |Im(mu)| ascending  (near-real-axis modes first)
  RASortAbsAscending,      // |mu|     ascending
  RASortAbsDescending,     // |mu|     descending (shift-invert: nearest sigma first)
  RASortRealAscending      // Re(mu)   ascending
};

template<class Field>
class RefinedArnoldi {
  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;
  typedef std::complex<double> Cd;

  LinearOperatorBase<Field>& M_;        // operator the Krylov space is built on
  GridBase*                  Grid_;
  RealD                      tol_;
  int                        verbose_;

  // raw-convergence operator: judge the residual on dW with lambda = sigma + 1/mu
  LinearOperatorBase<Field>* dW_    = nullptr;
  double                     sigma_ = 0.0;
  bool                       shift_ = false;

  std::vector<Field> V_;                // orthonormal Arnoldi basis (Euclidean)
  CMat               H_;                // (m+1) x m upper Hessenberg
  int                nSteps_ = 0;

  CVec               evals_;            // mu (eigenvalues of M)
  std::vector<Field> evecs_;            // refined Ritz vectors
  std::vector<RealD> residuals_;        // raw Euclidean residual on dW (convergence)

  static Cd cd(const ComplexD& z) { return Cd((double)real(z), (double)imag(z)); }
  Cd lamOf(Cd mu) const { return shift_ ? Cd(sigma_, 0.0) + 1.0/mu : mu; }

  static bool less(const Cd& a, const Cd& b, RefinedRitzSort s) {
    switch (s) {
      case RASortAbsImagAscending: return std::abs(a.imag()) < std::abs(b.imag());
      case RASortAbsAscending:     return std::abs(a)        < std::abs(b);
      case RASortAbsDescending:    return std::abs(a)        > std::abs(b);
      case RASortRealAscending:    return a.real()           < b.real();
    }
    return false;
  }

public:
  RefinedArnoldi(LinearOperatorBase<Field>& M, GridBase* grid, RealD tol = 1e-10, int verbose = 1)
    : M_(M), Grid_(grid), tol_(tol), verbose_(verbose) {}

  // Judge convergence on the raw operator dW (lambda = sigma + 1/mu if shiftInvert).
  void setRawCheck(LinearOperatorBase<Field>* dW, double sigma, bool shiftInvert) {
    dW_ = dW; sigma_ = sigma; shift_ = shiftInvert;
  }

  const CVec&               getEvals()     const { return evals_; }
  const std::vector<Field>& getEvecs()     const { return evecs_; }
  const std::vector<RealD>& getResiduals() const { return residuals_; }
  int                       getNumSteps()  const { return nSteps_; }

  // Build to dimension maxSteps and extract refined Ritz pairs.
  void operator()(const Field& v0, int maxSteps, RefinedRitzSort sort = RASortAbsImagAscending) {
    buildArnoldi(v0, maxSteps);
    extractRitzAt(nSteps_, sort);
  }

  // Single-vector Arnoldi factorisation  M V_m = V_{m+1} Hbar  (reorthogonalised twice).
  void buildArnoldi(const Field& v0, int maxSteps) {
    V_.clear(); V_.reserve(maxSteps);
    Field v(Grid_); v = v0; v = v * (1.0/std::sqrt(norm2(v))); V_.push_back(v);
    H_ = CMat::Zero(maxSteps + 1, maxSteps); nSteps_ = maxSteps;
    Field w(Grid_);
    for (int j = 0; j < maxSteps; j++) {
      M_.Op(V_[j], w);
      for (int pass = 0; pass < 2; pass++)
        for (int i = 0; i <= j; i++) {
          ComplexD h = innerProduct(V_[i], w);
          H_(i, j) += cd(h); w = w - V_[i] * h;
        }
      double hn = std::sqrt(norm2(w));
      H_(j + 1, j) = hn;
      if (hn < 1e-12) { nSteps_ = j + 1; break; }          // invariant subspace
      if (j + 1 < maxSteps) V_.push_back(w * (1.0 / hn));
    }
    if (verbose_)
      std::cout << GridLogMessage << "RefinedArnoldi: Krylov dim " << nSteps_ << std::endl;
  }

  // Refined Ritz extraction from the first m Arnoldi vectors, sorted by `sort`.
  void extractRitzAt(int m, RefinedRitzSort sort) {
    Eigen::ComplexEigenSolver<CMat> es(H_.block(0, 0, m, m), /*computeEigenvectors=*/false);
    CMat Hbar = H_.block(0, 0, m + 1, m);
    std::vector<Cd>  mus(m);
    std::vector<CVec> ys(m);
    for (int j = 0; j < m; j++) {
      Cd th = es.eigenvalues()(j); mus[j] = th;
      CMat C = Hbar; for (int i = 0; i < m; i++) C(i, i) -= th;
      Eigen::JacobiSVD<CMat> svd(C, Eigen::ComputeThinV);
      ys[j] = svd.matrixV().col(m - 1);                    // smallest singular vector
    }
    std::vector<int> idx(m); for (int j = 0; j < m; j++) idx[j] = j;
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return less(mus[a], mus[b], sort); });

    evals_.resize(m); evecs_.clear(); residuals_.clear();
    evecs_.reserve(m); residuals_.reserve(m);
    LinearOperatorBase<Field>* rop = dW_ ? dW_ : &M_;
    Field rw(Grid_), t(Grid_);
    for (int k = 0; k < m; k++) {
      int j = idx[k]; Cd mu = mus[j]; Cd lam = lamOf(mu);
      Field u(Grid_); u = Zero();
      for (int i = 0; i < m; i++) u = u + V_[i] * ys[j](i);
      evals_(k) = mu; evecs_.push_back(u);
      rop->Op(u, rw);
      ComplexD lf(lam.real(), lam.imag());
      t = rw - u * lf;
      residuals_.push_back(std::sqrt(norm2(t) / norm2(u)));
    }
  }

  // Count converged pairs (raw residual < tol).
  int numConverged() const {
    int n = 0; for (auto r : residuals_) if (r < tol_) n++; return n;
  }
};

NAMESPACE_END(Grid);
