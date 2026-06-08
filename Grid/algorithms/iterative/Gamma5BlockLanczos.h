/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/iterative/Gamma5BlockLanczos.h

    gamma5-Block Lanczos for the Wilson Dirac operator D_W.

    D_W is self-adjoint in the indefinite gamma5 inner product (u,v) = u^dag g5 v.
    This makes the block Krylov recurrence three-term (block-tridiagonal T_m) at
    block size 2, building span{ [v,w], M[v,w], M^2[v,w], ... } and targeting
    conjugate eigenvalue pairs together.  M is the operator supplied (D_W, or a
    shift-invert (D_W - sigma)^{-1} for interior eigenvalues).

    Ritz VALUES come from T_m.  Ritz VECTORS are extracted by REFINED Ritz: each
    vector is the EUCLIDEAN-residual minimiser over the Krylov subspace (a small
    generalised eigenproblem), which avoids the suboptimality of the gamma5
    (oblique) Galerkin projection while keeping the cheap gamma5 recurrence and the
    conjugate-pair block structure.  Convergence is judged by the raw Euclidean
    residual ||D_W u - lambda u|| / ||u|| (set via setRawCheck).

    A serious breakdown (a non-zero gamma5-neutral residual direction) is handled
    by look-ahead: the residual block is augmented with M applied to the neutral
    direction, restoring a non-degenerate gamma5-Gram; the block grows and T_m
    stays block-tridiagonal with variable block sizes.

    thickRestart() locks converged conjugate pairs, deflates them, and restarts
    from the live spillover block (never re-applying the raw operator to converged
    vectors, never dividing by a noise-level residual).

*************************************************************************************/
#ifndef GRID_GAMMA5_BLOCK_LANCZOS_H
#define GRID_GAMMA5_BLOCK_LANCZOS_H

#include <functional>
#include <numeric>
#include <vector>

NAMESPACE_BEGIN(Grid);

enum Gamma5RitzSort {
  G5SortAbsImagAscending,  // |Im(lambda)| ascending  (near-real-axis modes first)
  G5SortAbsAscending,      // |lambda|     ascending  (smallest magnitude first)
  G5SortAbsDescending,     // |lambda|     descending (shift-invert: nearest sigma)
  G5SortRealAscending      // Re(lambda)   ascending
};

template<class Field>
class Gamma5BlockLanczos {
public:
  using Gamma5Func = std::function<void(const Field&, Field&)>;

private:
  typedef Eigen::Matrix2cd CMat2;
  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;
  typedef std::complex<double> Cd;

  LinearOperatorBase<Field>& M_;       // operator the Krylov space is built on
  GridBase*                  Grid_;
  Gamma5Func                 g5_;
  RealD                      tol_;
  int                        verbose_;

  // Raw convergence operator: D_W and (for shift-invert) sigma, so the residual
  // is judged on the true operator, lambda = sigma + 1/mu.  Always set in use.
  LinearOperatorBase<Field>* dW_   = nullptr;
  double                     sigma_ = 0.0;
  bool                       shift_ = false;

  // Variable-size coefficient blocks: A_[k] (s_k x s_k), B_[k] (s_{k+1} x s_k),
  // C_[k] (s_{k-1} x s_k), G_[k] (s_k x s_k = +-1 signature).  s_k = Q_[k].size().
  std::vector<CMat> A_, B_, C_, G_;
  std::vector<std::vector<Field>> Q_;  // Q_[k] = columns of block k
  int  nSteps_ = 0;
  long lookaheads_ = 0;

  // Locked converged conjugate pairs (2-D blocks) for restart / deflation.
  std::vector<Field> lockV_;     // 2 columns per locked block
  std::vector<CMat2> lockG_;     // gamma5-Gram of each locked block
  std::vector<Cd>    lockEval_;  // 2 eigenvalues per locked block

  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

  // relative floor for a degenerate / neutral gamma5-Gram eigenvalue
  RealD degenRel_ = 1e-6;

  template<class C> static Cd toStd(const C& z) { return Cd((double)real(z), (double)imag(z)); }

public:
  Gamma5BlockLanczos(LinearOperatorBase<Field>& M, GridBase* grid,
                     Gamma5Func g5, RealD tol = 1e-8, int verbose = 1)
    : M_(M), Grid_(grid), g5_(g5), tol_(tol), verbose_(verbose) {}

  // Judge convergence on the raw operator dW (lambda = sigma + 1/mu if shiftInvert).
  void setRawCheck(LinearOperatorBase<Field>* dW, double sigma, bool shiftInvert) {
    dW_ = dW; sigma_ = sigma; shift_ = shiftInvert;
  }
  void setDegenRel(RealD r) { degenRel_ = r; }

  const CVec&               getEvals()     const { return evals_;     }
  const std::vector<Field>& getEvecs()     const { return evecs_;     }
  const std::vector<RealD>& getResiduals() const { return residuals_; }
  int                       getNumLocked() const { return (int)lockG_.size(); }
  int                       getNumSteps()  const { return nSteps_; }
  long                      getLookaheads()const { return lookaheads_; }

  // Single non-restarted pass (isolation test).
  void operator()(const Field& v0, const Field& v1, int maxSteps,
                  bool reorthog = true, Gamma5RitzSort sort = G5SortAbsImagAscending) {
    reset();
    if (!initStartBlock(v0, v1)) return;
    for (int step = 0; step < maxSteps; step++) {
      if (!lanczosStep(step, reorthog, /*deflate=*/false)) break;
      nSteps_ = step + 1;
      if (B_[step].norm() < tol_) break;
    }
    if (nSteps_ > 0) computeRitzPairs(nSteps_, sort);
  }

  // Re-extract Ritz pairs from the first m completed steps (nested subspaces).
  void extractRitzAt(int m, Gamma5RitzSort sort) {
    if (m >= 1 && m <= nSteps_) computeRitzPairs(m, sort);
  }

  // gamma5-metric thick restart with paired locking + deflation.  Returns the
  // number of converged conjugate pairs; locked pairs are spliced to the front
  // of the output arrays.
  int thickRestart(const Field& v0, const Field& v1,
                   int maxCycles, int cycleSteps, int nWantedPairs,
                   bool reorthog = true, Gamma5RitzSort sort = G5SortAbsImagAscending) {
    lockV_.clear(); lockG_.clear(); lockEval_.clear();
    Field s0(Grid_), s1(Grid_); s0 = v0; s1 = v1;
    for (int cyc = 0; cyc < maxCycles; cyc++) {
      reset();
      if (!initStartBlock(s0, s1)) break;
      for (int step = 0; step < cycleSteps; step++) {
        if (!lanczosStep(step, reorthog, /*deflate=*/true)) break;
        nSteps_ = step + 1;
        if (B_[step].norm() < tol_) break;
      }
      if (nSteps_ == 0) break;
      computeRitzPairs(nSteps_, sort);
      int newly = lockConvergedPairs(nWantedPairs);
      if (verbose_ > 0)
        std::cout << GridLogMessage << "[g5BL] cycle " << cyc << ": +" << newly
                  << " pairs, total " << lockG_.size() << "/" << nWantedPairs << std::endl;
      if ((int)lockG_.size() >= nWantedPairs) break;
      if (nSteps_ >= (int)Q_.size() || Q_[nSteps_].size() < 2) break;
      s0 = Q_[nSteps_][0]; s1 = Q_[nSteps_][1];
      deflateVec(s0); deflateVec(s1);
    }
    spliceLockedToFront();
    return (int)lockG_.size();
  }

private:
  void reset() { Q_.clear(); A_.clear(); B_.clear(); C_.clear(); G_.clear(); nSteps_ = 0; }
  int  sz(int k) const { return (int)Q_[k].size(); }

  std::vector<Field> applyOp(const std::vector<Field>& X) {
    std::vector<Field> Y; Y.reserve(X.size());
    for (auto& x : X) { Field y(Grid_); M_.Op(x, y); Y.push_back(y); }
    return Y;
  }
  // M(i,j) = X[i]^dag g5 Y[j]
  CMat g5Inner(const std::vector<Field>& X, const std::vector<Field>& Y) {
    int m = X.size(), n = Y.size();
    std::vector<Field> gY; gY.reserve(n);
    for (auto& y : Y) { Field gy(Grid_); g5_(y, gy); gY.push_back(gy); }
    CMat O(m, n);
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) O(i, j) = toStd(innerProduct(X[i], gY[j]));
    return O;
  }
  // R[j] -= sum_i X[i] * O(i,j)
  void subtractCombine(std::vector<Field>& R, const std::vector<Field>& X, const CMat& O) {
    for (int j = 0; j < (int)R.size(); j++)
      for (int i = 0; i < (int)X.size(); i++) R[j] = R[j] - X[i] * O(i, j);
  }
  Field combineCol(const std::vector<Field>& X, CVec c) {  // by value: materialises Eigen cols safely
    Field out(Grid_); out = Zero();
    for (int i = 0; i < (int)X.size(); i++) out = out + X[i] * c(i);
    return out;
  }

  bool initStartBlock(const Field& v0, const Field& v1) {
    Field u0(Grid_), u1(Grid_);
    u0 = v0;
    if (!lockG_.empty()) deflateVec(u0);
    RealD n = std::sqrt(norm2(u0));
    if (n < 1e-14) return false;
    u0 = u0 * (1.0 / n);
    u1 = v1;
    if (!lockG_.empty()) deflateVec(u1);
    auto proj = innerProduct(u0, u1);
    u1 = u1 - u0 * proj;
    n = std::sqrt(norm2(u1));
    if (n < 1e-14) return false;
    u1 = u1 * (1.0 / n);
    std::vector<Field> Q0 = {u0, u1};
    CMat G1 = g5Inner(Q0, Q0);
    Eigen::SelfAdjointEigenSolver<CMat> es(G1);
    if (std::abs(es.eigenvalues()(0)) < 1e-13 || std::abs(es.eigenvalues()(1)) < 1e-13) return false;
    G_.push_back(G1); Q_.push_back(Q0);
    return true;
  }

  // One block step: build Q_{k+1} from the residual block, with look-ahead.
  bool lanczosStep(int step, bool reorthog, bool deflate) {
    const std::vector<Field>& Qk = Q_[step];
    int s_k = Qk.size();
    CMat Gk = G_[step];

    std::vector<Field> P = applyOp(Qk);
    CMat A = Gk.inverse() * g5Inner(Qk, P);                  // A_k
    A_.push_back(A);
    CMat C = (step > 0) ? CMat(G_[step-1].inverse() * B_[step-1].adjoint() * Gk)
                        : CMat(CMat::Zero(0, s_k));          // C_k
    C_.push_back(C);

    std::vector<Field> R = P;                               // residual = M Q_k - Q_k A_k - Q_{k-1} C_k
    subtractCombine(R, Qk, A);
    if (step > 0) subtractCombine(R, Q_[step-1], C);
    if (reorthog)
      for (int j = 0; j <= step; j++) {
        CMat H = G_[j].inverse() * g5Inner(Q_[j], R);
        subtractCombine(R, Q_[j], H);
      }
    if (deflate && !lockG_.empty()) for (auto& r : R) deflateVec(r);

    // Look-ahead: while the residual block has a non-zero gamma5-neutral
    // direction, augment with M applied to it (escapes the neutral cone).
    const RealD relEps = degenRel_;
    std::vector<Field> S = R;
    for (int la = 0; la <= 3; la++) {
      Eigen::SelfAdjointEigenSolver<CMat> es(g5Inner(S, S));
      Eigen::VectorXd D = es.eigenvalues(); CMat U = es.eigenvectors();
      RealD dmax = D.cwiseAbs().maxCoeff();
      std::vector<int> neutral;
      for (int i = 0; i < D.size(); i++)
        if (std::abs(D(i)) < relEps * dmax && std::sqrt(norm2(combineCol(S, U.col(i)))) >= tol_)
          neutral.push_back(i);
      if (neutral.empty() || la == 3) break;
      std::vector<Field> add;
      for (int i : neutral) {
        Field ri = combineCol(S, U.col(i));
        ri = ri * (1.0 / std::sqrt(norm2(ri)));
        Field d(Grid_); M_.Op(ri, d);
        for (int j = 0; j <= step; j++) {                  // gamma5-orthogonalise vs history
          int sj = sz(j);
          Field g5d(Grid_); g5_(d, g5d);
          CVec pr(sj); for (int c = 0; c < sj; c++) pr(c) = toStd(innerProduct(Q_[j][c], g5d));
          CVec co = G_[j].inverse() * pr;
          for (int c = 0; c < sj; c++) d = d - Q_[j][c] * co(c);
        }
        RealD dn = std::sqrt(norm2(d));
        if (dn >= tol_) add.push_back(d * (1.0 / dn));
      }
      if (add.empty()) break;
      for (auto& a : add) S.push_back(a);
      lookaheads_++;
    }

    // Build Q_{k+1}: gamma5-orthonormal basis of the (augmented) residual block.
    Eigen::SelfAdjointEigenSolver<CMat> es(g5Inner(S, S));
    Eigen::VectorXd D = es.eigenvalues(); CMat U = es.eigenvectors();
    RealD dmax = D.cwiseAbs().maxCoeff();
    if (dmax < tol_ * tol_) return false;                  // happy breakdown
    RealD floor = std::max(relEps * dmax, tol_ * tol_);
    std::vector<int> keep;
    for (int i = 0; i < D.size(); i++) if (std::abs(D(i)) >= floor) keep.push_back(i);
    if (keep.empty()) return false;

    int s_kp1 = keep.size();
    std::vector<Field> Qkp1; Qkp1.reserve(s_kp1);
    CMat Gkp1 = CMat::Zero(s_kp1, s_kp1);
    for (int a = 0; a < s_kp1; a++) {
      int i = keep[a];
      Qkp1.push_back(combineCol(S, U.col(i)) * (1.0 / std::sqrt(std::abs(D(i)))));
      Gkp1(a, a) = Cd(D(i) > 0 ? 1.0 : -1.0, 0.0);
    }
    CMat Bkp1 = Gkp1.inverse() * g5Inner(Qkp1, R);         // R = Q_{k+1} B_{k+1}
    G_.push_back(Gkp1); B_.push_back(Bkp1); Q_.push_back(Qkp1);
    return true;
  }

  // Ritz values from T_m; Ritz vectors by refined (Euclidean) extraction;
  // residuals by the raw D_W Euclidean residual.
  void computeRitzPairs(int m, Gamma5RitzSort sort) {
    std::vector<int> off(m + 1, 0);
    for (int k = 0; k < m; k++) off[k+1] = off[k] + sz(k);
    int dim = off[m], s_m = sz(m);

    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(off[k], off[k], sz(k), sz(k)) = A_[k];
      if (k < m - 1) {
        Tm.block(off[k+1], off[k],   sz(k+1), sz(k)) = B_[k];
        Tm.block(off[k],   off[k+1], sz(k), sz(k+1)) = C_[k+1];
      }
    }
    Eigen::ComplexEigenSolver<CMat> ces(Tm);
    CVec lam = ces.eigenvalues();
    std::vector<int> idx(dim); std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return ritzLess(lam(a), lam(b), sort); });

    // Euclidean Grams of the augmented basis U = [V_m, Q_m] and of V_m.
    std::vector<const Field*> cols;
    for (int k = 0; k <= m; k++) for (int c = 0; c < sz(k); c++) cols.push_back(&Q_[k][c]);
    int dimA = (int)cols.size();
    CMat GE(dimA, dimA);
    for (int a = 0; a < dimA; a++)
      for (int b = a; b < dimA; b++) { Cd v = toStd(innerProduct(*cols[a], *cols[b])); GE(a,b)=v; GE(b,a)=std::conj(v); }
    CMat VE = GE.topLeftCorner(dim, dim);
    CMat Blink = CMat::Zero(s_m, dim);
    Blink.block(0, off[m-1], s_m, sz(m-1)) = B_[m-1];      // B_{m-1} E^T

    evals_.resize(dim); evecs_.clear(); residuals_.clear();
    evecs_.reserve(dim); residuals_.reserve(dim);
    for (int ji = 0; ji < dim; ji++) {
      Cd mu = lam(idx[ji]);
      // refined vector: min ||(M-mu) V_m z||/||V_m z||  =  smallest gen-eigpair (Krec^dag GE Krec, VE)
      CMat Krec(dim + s_m, dim);
      Krec.topRows(dim)    = Tm - mu * CMat::Identity(dim, dim);
      Krec.bottomRows(s_m) = Blink;
      CMat Pm = Krec.adjoint() * GE * Krec; Pm = 0.5 * (Pm + Pm.adjoint());
      Eigen::GeneralizedSelfAdjointEigenSolver<CMat> ges(Pm, VE);
      CVec z = ges.eigenvectors().col(0);
      Field u(Grid_); u = Zero();
      for (int k = 0; k < m; k++) for (int c = 0; c < sz(k); c++) u = u + Q_[k][c] * z(off[k] + c);
      evals_(ji) = ComplexD(mu.real(), mu.imag());
      evecs_.push_back(u);
      residuals_.push_back(0.0);   // filled by rawResidual below
    }
    rawResidual();
  }

  // residual = ||D_W u - lambda u|| / ||u||  (lambda = sigma + 1/mu in shift-invert)
  void rawResidual() {
    LinearOperatorBase<Field>* op = dW_ ? dW_ : &M_;
    Field w(Grid_);
    for (int i = 0; i < (int)evecs_.size(); i++) {
      Cd mu(real(evals_(i)), imag(evals_(i)));
      Cd lam = (dW_ && shift_) ? (sigma_ + 1.0/mu) : mu;
      op->Op(evecs_[i], w);
      typename Field::scalar_type lf(lam.real(), lam.imag());
      Field t(Grid_); t = w - evecs_[i] * lf;
      residuals_[i] = std::sqrt(norm2(t) / norm2(evecs_[i]));
    }
  }

  // Lock converged conjugate pairs (both members residual < tol).
  int lockConvergedPairs(int nWantedPairs) {
    int dim = (int)evecs_.size(), locked = 0;
    std::vector<bool> used(dim, false);
    for (int i = 0; i < dim; i++) {
      if ((int)lockG_.size() >= nWantedPairs) break;
      if (used[i] || residuals_[i] >= tol_) continue;
      Cd li(real(evals_(i)), imag(evals_(i)));
      int best = -1; double bd = 1e30;
      for (int j = 0; j < dim; j++) {
        if (j == i || used[j] || residuals_[j] >= tol_) continue;
        Cd lj(real(evals_(j)), imag(evals_(j)));
        double d = std::abs(lj - std::conj(li));
        if (d < bd) { bd = d; best = j; }
      }
      if (best < 0 || isAlreadyLocked(li)) { used[i] = (best>=0); continue; }
      std::vector<Field> blk = {evecs_[i], evecs_[best]};
      CMat Gg = g5Inner(blk, blk);
      CMat2 G; G << Gg(0,0), Gg(0,1), Gg(1,0), Gg(1,1);
      if (std::abs(G.determinant()) < 1e-12) { used[i] = used[best] = true; continue; }
      lockV_.push_back(evecs_[i]); lockV_.push_back(evecs_[best]);
      lockG_.push_back(G);
      lockEval_.push_back(li); lockEval_.push_back(Cd(real(evals_(best)), imag(evals_(best))));
      used[i] = used[best] = true; locked++;
    }
    return locked;
  }
  bool isAlreadyLocked(const Cd& lam) const {
    for (auto& e : lockEval_) if (std::abs(e - lam) < 10.0 * tol_) return true;
    return false;
  }
  // x <- x - sum_b V_b G_b^{-1} (V_b^dag g5 x)
  void deflateVec(Field& x) {
    for (size_t b = 0; b < lockG_.size(); b++) {
      const Field& w0 = lockV_[2*b]; const Field& w1 = lockV_[2*b+1];
      Field g5x(Grid_); g5_(x, g5x);
      Eigen::Vector2cd c; c(0) = toStd(innerProduct(w0, g5x)); c(1) = toStd(innerProduct(w1, g5x));
      Eigen::Vector2cd a = lockG_[b].inverse() * c;
      x = x - (w0 * a(0) + w1 * a(1));
    }
  }
  void spliceLockedToFront() {
    if (lockG_.empty()) return;
    int nl = 2 * (int)lockG_.size();
    CVec ev(nl + evals_.size());
    std::vector<Field> ec; std::vector<RealD> rs;
    for (size_t b = 0; b < lockG_.size(); b++) {
      ev(2*b)   = ComplexD(lockEval_[2*b].real(),   lockEval_[2*b].imag());
      ev(2*b+1) = ComplexD(lockEval_[2*b+1].real(), lockEval_[2*b+1].imag());
      ec.push_back(lockV_[2*b]); ec.push_back(lockV_[2*b+1]); rs.push_back(0.0); rs.push_back(0.0);
    }
    for (int i = 0; i < (int)evals_.size(); i++) ev(nl + i) = evals_(i);
    for (auto& v : evecs_) ec.push_back(v);
    for (auto& r : residuals_) rs.push_back(r);
    evals_ = ev; evecs_ = ec; residuals_ = rs;
  }

  static bool ritzLess(const ComplexD& a, const ComplexD& b, Gamma5RitzSort sort) {
    Cd x(real(a), imag(a)), y(real(b), imag(b));
    switch (sort) {
      case G5SortAbsImagAscending: return std::abs(x.imag()) < std::abs(y.imag());
      case G5SortAbsAscending:     return std::abs(x)        < std::abs(y);
      case G5SortAbsDescending:    return std::abs(x)        > std::abs(y);
      case G5SortRealAscending:    return x.real()           < y.real();
    }
    return std::abs(x.imag()) < std::abs(y.imag());
  }
};

NAMESPACE_END(Grid);

#endif
