/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/iterative/Gamma5BlockLanczos.h

    gamma5-Block Lanczos for gamma5-Hermitian operators (Wilson Dirac D_W).

    Reference: S. Yamamoto, "gamma5-Block Krylov (Block Lanczos) Methods for the
    Wilson Dirac Operator" (2026).

    D_W satisfies D_W^dag = g5 D_W g5, so it is self-adjoint in the indefinite
    gamma5-inner product (u,v) = u^dag g5 v.  Block size s = 2.  The three-term
    block recurrence

        Q_{k+1} B_{k+1} = D_W Q_k - Q_k A_k - Q_{k-1} C_k,   Q_0 = 0,

    builds a block-tridiagonal projected matrix T_m whose eigenvalues approximate
    eigenvalues of D_W directly (NOT those of H_W = g5 D_W).

    This file provides:
      operator()  : single non-restarted pass (the "isolation test" path).
      thickRestart: gamma5-metric thick restart with paired locking, which avoids
                    the noise-normalisation singularity and doubler re-amplification
                    that wreck explicit restart (see reference, Sec. 9 / App. B).

    Verified inner kernel (lanczosStep / computeRitzPairs / Gram helpers) follows
    the reference equations directly.  All diagnostic printing is gated by a
    verbosity level and kept minimal.

*************************************************************************************/
#ifndef GRID_GAMMA5_BLOCK_LANCZOS_H
#define GRID_GAMMA5_BLOCK_LANCZOS_H

#include <functional>
#include <numeric>
#include <iomanip>
#include <vector>

NAMESPACE_BEGIN(Grid);

// Ritz selection criterion for sorting the 2m eigenvalues of T_m.
enum Gamma5RitzSort {
  G5SortAbsImagAscending,  // |Im(lambda)| ascending  (near-real-axis / physical modes first)
  G5SortAbsAscending,      // |lambda|     ascending  (smallest magnitude first)
  G5SortAbsDescending,     // |lambda|     descending (largest magnitude first; shift-invert)
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

  LinearOperatorBase<Field>& Linop_;
  GridBase*                  Grid_;
  Gamma5Func                 g5_;
  RealD                      tol_;
  int                        verbose_;

  // Per-step 2x2 coefficient blocks.  Storage index k = paper step (k+1):
  //   A_blk[k] = A_{k+1},  B_blk[k] = B_{k+2},  C_blk[k] = C_{k+1},  G_blk[k] = G_{k+1}.
  std::vector<CMat2> A_blk, B_blk, C_blk, G_blk;

  // Krylov basis: basis_[2k], basis_[2k+1] are the two columns of Q_{k+1}.
  std::vector<Field> basis_;
  int nSteps_;

  // Locked (converged / retained) 2-D blocks for thick restart / deflation.
  // Each locked block is a gamma5-non-degenerate 2-space (a conjugate pair).
  std::vector<Field> lockV_;     // 2 columns per locked block, concatenated
  std::vector<CMat2> lockG_;     // gamma5-Gram of each locked block (2x2, invertible)
  std::vector<Cd>    lockEval_;  // 2 Ritz values per locked block

  // Output (sorted)
  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

  // Accept either ComplexD or ComplexF (Grid innerProduct precision depends on Field).
  template<class C> static inline Cd toStd(const C& z) {
    return Cd((double)real(z), (double)imag(z));
  }

public:
  Gamma5BlockLanczos(LinearOperatorBase<Field>& op, GridBase* grid,
                     Gamma5Func g5, RealD tol = 1e-8, int verbose = 1)
    : Linop_(op), Grid_(grid), g5_(g5), tol_(tol), verbose_(verbose), nSteps_(0) {}

  const CVec&               getEvals()     const { return evals_;     }
  const std::vector<Field>& getEvecs()     const { return evecs_;     }
  const std::vector<RealD>& getResiduals() const { return residuals_; }
  int                       getNumLocked() const { return (int)lockG_.size(); }

  void log(const std::string& s) const {
    if (verbose_ > 0) std::cout << GridLogMessage << "[g5BL] " << s << std::endl;
  }

  // ------------------------------------------------------------------
  // Single non-restarted pass.  This is the ISOLATION TEST path: the
  // verified inner kernel, no restart logic.  Builds up to maxSteps blocks
  // from the starting block Q_1 = [v0', v1'] (gamma5-Gram-Schmidt of v0,v1),
  // then extracts 2*nSteps Ritz pairs.
  // ------------------------------------------------------------------
  void operator()(const Field& v0, const Field& v1, int maxSteps,
                  bool reorthog = false, Gamma5RitzSort sort = G5SortAbsImagAscending)
  {
    reset();
    if (!initStartBlock(v0, v1)) return;

    for (int step = 0; step < maxSteps; step++) {
      bool ok = lanczosStep(step, reorthog, /*deflate=*/false);
      if (!ok) break;
      nSteps_ = step + 1;
      RealD beta = B_blk[step].norm();
      if (verbose_ > 1)
        log("step " + std::to_string(step) + "  ||B|| = " + std::to_string(beta));
      if (beta < tol_) { log("beta < tol; stopping at step " + std::to_string(step)); break; }
    }
    if (nSteps_ == 0) return;
    computeRitzPairs(nSteps_, sort);
  }

  // ------------------------------------------------------------------
  // gamma5-metric thick restart.
  //
  // Each cycle runs `cycleSteps` block-Lanczos steps (deflated against the
  // locked blocks), extracts Ritz pairs, and LOCKS converged conjugate pairs.
  // The next cycle restarts from the live spillover block Q_{m+1} (NOT from a
  // sum of Ritz vectors), gamma5-deflated against the locked space.  This
  // preserves the implicit polynomial filter, never applies raw D_W to a
  // converged Ritz vector, and never divides by a near-zero (noise) residual,
  // so it avoids the explicit-restart failure modes.
  //
  // Returns the number of locked (converged) eigen-pairs found, and leaves the
  // full sorted Ritz set of the LAST cycle in evals_/evecs_/residuals_ (locked
  // pairs are spliced to the front, residuals first).
  // ------------------------------------------------------------------
  int thickRestart(const Field& v0, const Field& v1,
                   int maxCycles, int cycleSteps, int nWantedPairs,
                   bool reorthog = true, Gamma5RitzSort sort = G5SortAbsImagAscending)
  {
    lockV_.clear(); lockG_.clear(); lockEval_.clear();

    Field s0(Grid_), s1(Grid_);
    s0 = v0; s1 = v1;

    for (int cyc = 0; cyc < maxCycles; cyc++) {
      reset();
      if (!initStartBlock(s0, s1)) {
        log("thickRestart: degenerate start block on cycle " + std::to_string(cyc));
        break;
      }

      int done = 0;
      for (int step = 0; step < cycleSteps; step++) {
        bool ok = lanczosStep(step, reorthog, /*deflate=*/true);
        if (!ok) break;
        nSteps_ = step + 1;
        done = step + 1;
        if (B_blk[step].norm() < tol_) break;
      }
      if (nSteps_ == 0) { log("thickRestart: no steps taken; stopping."); break; }

      computeRitzPairs(nSteps_, sort);

      // Lock newly-converged conjugate pairs from this cycle's Ritz set.
      int newlyLocked = lockConvergedPairs(nWantedPairs);

      int nLocked = (int)lockG_.size();
      log("cycle " + std::to_string(cyc) + ": steps=" + std::to_string(done)
          + "  newly locked=" + std::to_string(newlyLocked)
          + "  total locked pairs=" + std::to_string(nLocked)
          + " / " + std::to_string(nWantedPairs));

      if (nLocked >= nWantedPairs) {
        log("thickRestart: converged " + std::to_string(nLocked) + " pairs.");
        break;
      }

      // Next start = live spillover block Q_{m+1}, deflated against locked space.
      // Q_{m+1} columns are basis_[2*nSteps_], basis_[2*nSteps_+1].
      int idx = 2 * nSteps_;
      if (idx + 1 >= (int)basis_.size()) {
        log("thickRestart: no spillover block available; stopping.");
        break;
      }
      s0 = basis_[idx];
      s1 = basis_[idx + 1];
      deflateBlock(s0, s1);
    }

    spliceLockedToFront();
    return (int)lockG_.size();
  }

private:
  // ---- setup ----

  void reset() {
    basis_.clear();
    A_blk.clear(); B_blk.clear(); C_blk.clear(); G_blk.clear();
    nSteps_ = 0;
  }

  // Build Q_1 = [u0,u1] from (v0,v1): normalise u0, L2-orthogonalise u1 vs u0,
  // (optionally deflate both against locked space), check G_1 invertible.
  bool initStartBlock(const Field& v0, const Field& v1) {
    Field u0(Grid_), u1(Grid_);
    u0 = v0;
    if (!lockG_.empty()) deflateVec(u0);
    RealD n = std::sqrt(norm2(u0));
    if (n < 1e-14) { log("init: first start vector vanished"); return false; }
    u0 = u0 * (1.0 / n);

    u1 = v1;
    if (!lockG_.empty()) deflateVec(u1);
    auto proj = innerProduct(u0, u1);
    u1 = u1 - u0 * proj;
    n = std::sqrt(norm2(u1));
    if (n < 1e-14) { log("init: second start vector linearly dependent"); return false; }
    u1 = u1 * (1.0 / n);

    CMat2 G1 = gram(u0, u1);
    Eigen::SelfAdjointEigenSolver<CMat2> es(G1);
    auto ev = es.eigenvalues();
    if (std::abs(ev(0)) < 1e-13 || std::abs(ev(1)) < 1e-13) {
      log("init: degenerate start (G1 eigenvalue ~ 0); pick a different seed");
      return false;
    }
    G_blk.push_back(G1);
    basis_.push_back(u0);
    basis_.push_back(u1);
    return true;
  }

  // ---- core three-term block step (verified kernel) ----

  // One block-Lanczos step.  On success pushes Q_{step+2} and returns true.
  // If deflate==true, the residual block is gamma5-projected against the locked
  // space before normalisation.
  bool lanczosStep(int step, bool reorthog, bool deflate) {
    const Field& q1 = basis_[2*step];
    const Field& q2 = basis_[2*step + 1];
    CMat2 Gk = G_blk[step];

    // (i) matvecs P_k = D_W Q_k
    Field p1(Grid_), p2(Grid_);
    Linop_.Op(q1, p1);
    Linop_.Op(q2, p2);

    // (ii) M_k = Q_k^dag g5 P_k ;  (iii) A_k = G_k^{-1} M_k
    CMat2 Mk = g5Block(q1, q2, p1, p2);
    CMat2 Ak = Gk.inverse() * Mk;
    A_blk.push_back(Ak);

    // (iv) C_k = G_{k-1}^{-1} B_k^dag G_k  (zero at step 0)
    CMat2 Ck = CMat2::Zero();
    if (step > 0) Ck = G_blk[step-1].inverse() * B_blk[step-1].adjoint() * Gk;
    C_blk.push_back(Ck);

    // (v) residual R_k = P_k - Q_k A_k - Q_{k-1} C_k  (columnwise)
    Field r1(Grid_), r2(Grid_);
    r1 = p1 - (q1 * Ak(0,0) + q2 * Ak(1,0));
    r2 = p2 - (q1 * Ak(0,1) + q2 * Ak(1,1));
    if (step > 0) {
      const Field& s1 = basis_[2*(step-1)];
      const Field& s2 = basis_[2*(step-1)+1];
      r1 = r1 - (s1 * Ck(0,0) + s2 * Ck(1,0));
      r2 = r2 - (s1 * Ck(0,1) + s2 * Ck(1,1));
    }

    // optional full gamma5-reorthogonalisation against all previous blocks
    if (reorthog) {
      for (int j = 0; j <= step; j++) {
        CMat2 Mj = g5Block(basis_[2*j], basis_[2*j+1], r1, r2);
        CMat2 Hj = G_blk[j].inverse() * Mj;
        r1 = r1 - (basis_[2*j] * Hj(0,0) + basis_[2*j+1] * Hj(1,0));
        r2 = r2 - (basis_[2*j] * Hj(0,1) + basis_[2*j+1] * Hj(1,1));
      }
    }

    // deflate residual against locked converged blocks
    if (deflate && !lockG_.empty()) deflateBlock(r1, r2);

    // (vi) Gamma_k = R_k^dag g5 R_k ; eigendecompose
    CMat2 Gamma = gram(r1, r2);
    Eigen::SelfAdjointEigenSolver<CMat2> es(Gamma);
    Eigen::Vector2d D = es.eigenvalues();
    CMat2           U = es.eigenvectors();

    // breakdown checks (reference Eq. happy/serious breakdown)
    const RealD breakdownEps = 1e-14;
    for (int j = 0; j < 2; j++) {
      Field rj(Grid_);
      rj = r1 * U(0,j) + r2 * U(1,j);
      RealD r2n = norm2(rj);
      RealD rn  = std::sqrt(r2n);
      if (rn < tol_) {
        if (verbose_ > 1) log("happy breakdown at step " + std::to_string(step)
                              + " dir " + std::to_string(j));
        return false;
      }
      if (std::abs(D(j)) < breakdownEps * r2n) {
        log("serious breakdown at step " + std::to_string(step) + " dir " + std::to_string(j)
            + " (neutral residual; look-ahead not implemented)");
        return false;
      }
    }

    // (vii) indefinite LDL^dag: G_{k+1}=sign(D), B_{k+1}=|D|^{1/2}U^dag, Q_{k+1}=R U|D|^{-1/2}
    CMat2 Gkp1 = CMat2::Zero();
    Gkp1(0,0) = Cd((D(0) > 0.0) ? 1.0 : -1.0, 0.0);
    Gkp1(1,1) = Cd((D(1) > 0.0) ? 1.0 : -1.0, 0.0);
    double sq0 = std::sqrt(std::abs(D(0)));
    double sq1 = std::sqrt(std::abs(D(1)));

    CMat2 Bkp1;
    Bkp1.row(0) = U.col(0).adjoint() * sq0;
    Bkp1.row(1) = U.col(1).adjoint() * sq1;

    Field qn1(Grid_), qn2(Grid_);
    qn1 = (r1 * U(0,0) + r2 * U(1,0)) * (1.0 / sq0);
    qn2 = (r1 * U(0,1) + r2 * U(1,1)) * (1.0 / sq1);

    G_blk.push_back(Gkp1);
    B_blk.push_back(Bkp1);
    basis_.push_back(qn1);
    basis_.push_back(qn2);
    return true;
  }

  // ---- Ritz extraction ----

  void computeRitzPairs(int m, Gamma5RitzSort sort) {
    int dim = 2 * m;
    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(2*k, 2*k, 2, 2) = A_blk[k];
      if (k < m - 1) {
        Tm.block(2*k+2, 2*k,   2, 2) = B_blk[k];
        Tm.block(2*k,   2*k+2, 2, 2) = C_blk[k+1];
      }
    }

    Eigen::ComplexEigenSolver<CMat> ces(Tm);
    CVec lam = ces.eigenvalues();
    CMat Y   = ces.eigenvectors();

    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b){ return ritzLess(lam(a), lam(b), sort); });

    const CMat2& Bm1 = B_blk[m-1];  // B_{m+1}
    evals_.resize(dim);
    evecs_.clear();
    residuals_.clear();
    evecs_.reserve(dim);
    residuals_.reserve(dim);

    for (int ji = 0; ji < dim; ji++) {
      int j = idx[ji];
      evals_(ji) = lam(j);
      CVec yj = Y.col(j);

      // Ritz vector u_j = V_m y_j
      Field uj(Grid_); uj = Zero();
      for (int k = 0; k < m; k++) {
        uj = uj + basis_[2*k]   * yj(2*k);
        uj = uj + basis_[2*k+1] * yj(2*k+1);
      }
      evecs_.push_back(uj);

      // true residual r_j = Q_{m+1} B_{m+1} tau_j, tau_j = last 2 entries of y_j
      Eigen::Vector2cd tau(yj(dim-2), yj(dim-1));
      Eigen::Vector2cd Bt = Bm1 * tau;
      Field rj(Grid_);
      rj = basis_[2*m] * Bt(0) + basis_[2*m+1] * Bt(1);
      residuals_.push_back(std::sqrt(norm2(rj)));
    }
  }

  // ---- locking / deflation for thick restart ----

  // Lock conjugate pairs from the current Ritz set that have converged (both
  // members residual < tol) and that are not already locked.  Returns count.
  int lockConvergedPairs(int nWantedPairs) {
    int dim = (int)evecs_.size();
    int locked = 0;
    std::vector<bool> used(dim, false);

    for (int i = 0; i < dim; i++) {
      if ((int)lockG_.size() >= nWantedPairs) break;
      if (used[i] || residuals_[i] >= tol_) continue;
      // find conjugate partner j with closest conj(lambda_i) and small residual
      Cd li = toStd(ComplexD(real(evals_(i)), imag(evals_(i))));
      int best = -1; double bestd = 1e30;
      for (int j = 0; j < dim; j++) {
        if (j == i || used[j] || residuals_[j] >= tol_) continue;
        Cd lj = toStd(ComplexD(real(evals_(j)), imag(evals_(j))));
        double d = std::abs(lj - std::conj(li));
        if (d < bestd) { bestd = d; best = j; }
      }
      if (best < 0) continue;  // no partner; skip (single real mode handled below)
      if (isAlreadyLocked(li)) { used[i] = used[best] = true; continue; }

      // Build the 2-column locked block from Ritz vectors i and best.
      Field w0(Grid_), w1(Grid_);
      w0 = evecs_[i];
      w1 = evecs_[best];
      CMat2 G = gram(w0, w1);
      // require gamma5-non-degenerate
      if (std::abs(G.determinant()) < 1e-12) { used[i] = used[best] = true; continue; }

      lockV_.push_back(w0);
      lockV_.push_back(w1);
      lockG_.push_back(G);
      lockEval_.push_back(toStd(ComplexD(real(evals_(i)),    imag(evals_(i)))));
      lockEval_.push_back(toStd(ComplexD(real(evals_(best)), imag(evals_(best)))));
      used[i] = used[best] = true;
      locked++;
    }
    return locked;
  }

  bool isAlreadyLocked(const Cd& lam) const {
    for (size_t b = 0; b < lockEval_.size(); b++)
      if (std::abs(lockEval_[b] - lam) < 10.0 * tol_) return true;
    return false;
  }

  // gamma5-oblique projection of a single vector onto the complement of the
  // locked space:  x <- x - sum_b V_b G_b^{-1} (V_b^dag g5 x).
  void deflateVec(Field& x) {
    for (size_t b = 0; b < lockG_.size(); b++) {
      const Field& w0 = lockV_[2*b];
      const Field& w1 = lockV_[2*b+1];
      Field g5x(Grid_); g5_(x, g5x);
      Eigen::Vector2cd c;
      c(0) = toStd(innerProduct(w0, g5x));
      c(1) = toStd(innerProduct(w1, g5x));
      Eigen::Vector2cd a = lockG_[b].inverse() * c;
      x = x - (w0 * a(0) + w1 * a(1));
    }
  }

  void deflateBlock(Field& x0, Field& x1) {
    deflateVec(x0);
    deflateVec(x1);
  }

  // Put locked pairs (residual 0, exact) at the front of the output arrays.
  void spliceLockedToFront() {
    if (lockG_.empty()) return;
    int nl = 2 * (int)lockG_.size();
    CVec ev(nl + evals_.size());
    std::vector<Field> ec; ec.reserve(nl + evecs_.size());
    std::vector<RealD> rs; rs.reserve(nl + residuals_.size());
    for (size_t b = 0; b < lockG_.size(); b++) {
      ev(2*b)   = ComplexD(lockEval_[2*b].real(),   lockEval_[2*b].imag());
      ev(2*b+1) = ComplexD(lockEval_[2*b+1].real(), lockEval_[2*b+1].imag());
      ec.push_back(lockV_[2*b]);
      ec.push_back(lockV_[2*b+1]);
      rs.push_back(0.0);
      rs.push_back(0.0);
    }
    for (int i = 0; i < (int)evals_.size(); i++) ev(nl + i) = evals_(i);
    for (auto& v : evecs_)    ec.push_back(v);
    for (auto& r : residuals_) rs.push_back(r);
    evals_ = ev; evecs_ = ec; residuals_ = rs;
  }

  // ---- gamma5 Gram helpers ----

  // G[i,j] = q^(i)^dag g5 q^(j) for block [q1,q2].
  CMat2 gram(const Field& q1, const Field& q2) {
    Field g1(Grid_), g2(Grid_);
    g5_(q1, g1); g5_(q2, g2);
    CMat2 G;
    G(0,0) = toStd(innerProduct(q1, g1));
    G(0,1) = toStd(innerProduct(q1, g2));
    G(1,0) = toStd(innerProduct(q2, g1));
    G(1,1) = toStd(innerProduct(q2, g2));
    return G;
  }

  // M[i,j] = q^(i)^dag g5 p^(j).
  CMat2 g5Block(const Field& q1, const Field& q2, const Field& p1, const Field& p2) {
    Field g1(Grid_), g2(Grid_);
    g5_(p1, g1); g5_(p2, g2);
    CMat2 M;
    M(0,0) = toStd(innerProduct(q1, g1));
    M(0,1) = toStd(innerProduct(q1, g2));
    M(1,0) = toStd(innerProduct(q2, g1));
    M(1,1) = toStd(innerProduct(q2, g2));
    return M;
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
