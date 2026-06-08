/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/algorithms/iterative/Gamma5BlockLanczos.h

    gamma5-Block Lanczos for gamma5-Hermitian operators (Wilson Dirac D_W).

    Reference: S. Yamamoto, "gamma5-Block Krylov (Block Lanczos) Methods for the
    Wilson Dirac Operator" (2026).

    D_W is self-adjoint in the indefinite gamma5-inner product (u,v) = u^dag g5 v.
    The block recurrence

        Q_{k+1} B_{k+1} = D_W Q_k - Q_k A_k - Q_{k-1} C_k,   Q_0 = 0,

    builds a block-tridiagonal projected matrix T_m whose eigenvalues approximate
    eigenvalues of D_W directly (NOT those of H_W = g5 D_W).

    Block sizes are VARIABLE: nominally 2 (a chiral pair), but a serious breakdown
    (a non-zero gamma5-neutral residual direction) triggers LOOK-AHEAD -- the
    residual block is augmented with D_W applied to the neutral direction, which
    generically escapes the neutral cone, restoring a non-degenerate gamma5-Gram.
    The block then grows; T_m remains block-tridiagonal with variable block sizes.

*************************************************************************************/
#ifndef GRID_GAMMA5_BLOCK_LANCZOS_H
#define GRID_GAMMA5_BLOCK_LANCZOS_H

#include <functional>
#include <numeric>
#include <iomanip>
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

  LinearOperatorBase<Field>& Linop_;
  GridBase*                  Grid_;
  Gamma5Func                 g5_;
  RealD                      tol_;
  int                        verbose_;

  // Variable-size coefficient blocks: A_[k] (s_k x s_k), B_[k] (s_{k+1} x s_k),
  // C_[k] (s_{k-1} x s_k), G_[k] (s_k x s_k = signature).  s_k = Q_[k].size().
  std::vector<CMat> A_, B_, C_, G_;

  // Krylov basis: Q_[k] is the list of columns of block k.
  std::vector<std::vector<Field>> Q_;
  int nSteps_;

  // Locked (converged) 2-D blocks for thick restart / deflation.
  std::vector<Field> lockV_;     // 2 columns per locked block, concatenated
  std::vector<CMat2> lockG_;     // gamma5-Gram of each locked block (2x2)
  std::vector<Cd>    lockEval_;  // 2 Ritz values per locked block

  CVec               evals_;
  std::vector<Field> evecs_;
  std::vector<RealD> residuals_;

  long lookaheadCount_ = 0;      // diagnostic: number of look-ahead expansions used

  // per-step diagnostics (manuscript Sec. 8): oblique-projector conditioning and
  // loss of gamma5-orthogonality.
  std::vector<double> kappaGamma_;   // kappa(Gamma_k) = max|d|/min|d| of residual Gram
  std::vector<double> etaLoss_;      // ||Q_1^dag g5 Q_{k+1}||_F  (should be ~0)
  std::vector<double> cycleBestRes_; // best (smallest) Ritz residual at end of each restart cycle

  template<class C> static inline Cd toStd(const C& z) {
    return Cd((double)real(z), (double)imag(z));
  }

  // relative threshold for a "degenerate / neutral" gamma5-Gram eigenvalue.
  // Set ~ sqrt(machine eps) of the working precision; the example runs single
  // precision (eps ~ 1e-7), so 1e-6 is a safe, precision-aware floor that is
  // also conservative for a double build.  Configurable via setDegenRel().
  RealD degenRel_ = 1e-6;
  RealD degenRel() const { return degenRel_; }

public:
  Gamma5BlockLanczos(LinearOperatorBase<Field>& op, GridBase* grid,
                     Gamma5Func g5, RealD tol = 1e-8, int verbose = 1)
    : Linop_(op), Grid_(grid), g5_(g5), tol_(tol), verbose_(verbose), nSteps_(0) {}

  const CVec&               getEvals()      const { return evals_;     }
  const std::vector<Field>& getEvecs()      const { return evecs_;     }
  const std::vector<RealD>& getResiduals()  const { return residuals_; }
  int                       getNumLocked()  const { return (int)lockG_.size(); }
  long                      getLookaheads() const { return lookaheadCount_; }
  void                      setDegenRel(RealD r)   { degenRel_ = r; }
  const std::vector<double>& getKappaGamma() const { return kappaGamma_; }
  const std::vector<double>& getEtaLoss()    const { return etaLoss_;    }
  const std::vector<double>& getCycleBestRes() const { return cycleBestRes_; }
  int getNumSteps() const { return nSteps_; }
  // Re-extract Ritz pairs using only the first m completed steps (the Krylov
  // subspaces are nested), for residual-vs-Krylov-dimension histories.
  void extractRitzAt(int m, Gamma5RitzSort sort) { if (m >= 1 && m <= nSteps_) computeRitzPairs(m, sort); }

  void log(const std::string& s) const {
    if (verbose_ > 0) std::cout << GridLogMessage << "[g5BL] " << s << std::endl;
  }

  // ---------------- single non-restarted pass (isolation test) ----------------
  void operator()(const Field& v0, const Field& v1, int maxSteps,
                  bool reorthog = false, Gamma5RitzSort sort = G5SortAbsImagAscending)
  {
    reset();
    if (!initStartBlock(v0, v1)) return;
    for (int step = 0; step < maxSteps; step++) {
      if (!lanczosStep(step, reorthog, /*deflate=*/false)) break;
      nSteps_ = step + 1;
      if (B_[step].norm() < tol_) { log("beta<tol; stop at step "+std::to_string(step)); break; }
    }
    if (nSteps_ == 0) return;
    computeRitzPairs(nSteps_, sort);
  }

  // ---------------- gamma5-metric thick restart with paired locking ----------------
  int thickRestart(const Field& v0, const Field& v1,
                   int maxCycles, int cycleSteps, int nWantedPairs,
                   bool reorthog = true, Gamma5RitzSort sort = G5SortAbsImagAscending)
  {
    lockV_.clear(); lockG_.clear(); lockEval_.clear(); cycleBestRes_.clear();
    Field s0(Grid_), s1(Grid_); s0 = v0; s1 = v1;

    for (int cyc = 0; cyc < maxCycles; cyc++) {
      reset();
      if (!initStartBlock(s0, s1)) { log("thickRestart: degenerate start, cycle "+std::to_string(cyc)); break; }

      for (int step = 0; step < cycleSteps; step++) {
        if (!lanczosStep(step, reorthog, /*deflate=*/true)) break;
        nSteps_ = step + 1;
        if (B_[step].norm() < tol_) break;
      }
      if (nSteps_ == 0) { log("thickRestart: no steps; stop."); break; }
      computeRitzPairs(nSteps_, sort);
      { double best = 1e300; for (auto r : residuals_) best = std::min(best, r);
        cycleBestRes_.push_back(best); }
      int newly = lockConvergedPairs(nWantedPairs);
      int nLk = (int)lockG_.size();
      log("cycle "+std::to_string(cyc)+": newly locked "+std::to_string(newly)
          +"  total "+std::to_string(nLk)+"/"+std::to_string(nWantedPairs));
      if (nLk >= nWantedPairs) { log("converged "+std::to_string(nLk)+" pairs."); break; }

      // restart from the live spillover block Q_[nSteps_], deflated against locked
      if (nSteps_ >= (int)Q_.size() || Q_[nSteps_].size() < 2) { log("no spillover; stop."); break; }
      s0 = Q_[nSteps_][0];
      s1 = Q_[nSteps_][1];
      deflateVec(s0); deflateVec(s1);
    }
    spliceLockedToFront();
    return (int)lockG_.size();
  }

private:
  void reset() { Q_.clear(); A_.clear(); B_.clear(); C_.clear(); G_.clear();
                 kappaGamma_.clear(); etaLoss_.clear(); nSteps_ = 0; }

  int sz(int k) const { return (int)Q_[k].size(); }

  // --- linear-algebra helpers on column lists ---

  // D_W applied to each column.
  std::vector<Field> applyOp(const std::vector<Field>& X) {
    std::vector<Field> Y; Y.reserve(X.size());
    for (auto& x : X) { Field y(Grid_); Linop_.Op(x, y); Y.push_back(y); }
    return Y;
  }

  // M(i,j) = X[i]^dag g5 Y[j]   (|X| x |Y|)
  CMat g5Inner(const std::vector<Field>& X, const std::vector<Field>& Y) {
    int m = X.size(), n = Y.size();
    std::vector<Field> gY; gY.reserve(n);
    for (auto& y : Y) { Field gy(Grid_); g5_(y, gy); gY.push_back(gy); }
    CMat M(m, n);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        M(i, j) = toStd(innerProduct(X[i], gY[j]));
    return M;
  }

  // R[j] -= sum_i X[i] * M(i,j)
  void subtractCombine(std::vector<Field>& R, const std::vector<Field>& X, const CMat& M) {
    for (int j = 0; j < (int)R.size(); j++)
      for (int i = 0; i < (int)X.size(); i++)
        R[j] = R[j] - X[i] * M(i, j);
  }

  // single linear combination col = sum_i X[i] * c(i).
  // Take c BY VALUE so an Eigen column expression (U.col(i)) is safely
  // materialised into a VectorXcd at the call (binding a const ref to the
  // Block temporary is unsafe).
  Field combineCol(const std::vector<Field>& X, CVec c) {
    Field out(Grid_); out = Zero();
    for (int i = 0; i < (int)X.size(); i++) out = out + X[i] * c(i);
    return out;
  }

  // --- setup ---
  bool initStartBlock(const Field& v0, const Field& v1) {
    Field u0(Grid_), u1(Grid_);
    u0 = v0;
    if (!lockG_.empty()) deflateVec(u0);
    RealD n = std::sqrt(norm2(u0));
    if (n < 1e-14) { log("init: first vector vanished"); return false; }
    u0 = u0 * (1.0 / n);
    u1 = v1;
    if (!lockG_.empty()) deflateVec(u1);
    auto proj = innerProduct(u0, u1);
    u1 = u1 - u0 * proj;
    n = std::sqrt(norm2(u1));
    if (n < 1e-14) { log("init: second vector dependent"); return false; }
    u1 = u1 * (1.0 / n);

    std::vector<Field> Q0 = {u0, u1};
    CMat G1 = g5Inner(Q0, Q0);
    Eigen::SelfAdjointEigenSolver<CMat> es(G1);
    auto ev = es.eigenvalues();
    if (std::abs(ev(0)) < 1e-13 || std::abs(ev(1)) < 1e-13) {
      log("init: degenerate start (G1 ~ singular)"); return false;
    }
    G_.push_back(G1); Q_.push_back(Q0);
    return true;
  }

  // --- core block step with look-ahead ---
  bool lanczosStep(int step, bool reorthog, bool deflate) {
    const std::vector<Field>& Qk = Q_[step];
    int s_k = Qk.size();
    CMat Gk = G_[step];

    std::vector<Field> P = applyOp(Qk);             // D_W Q_k
    CMat M = g5Inner(Qk, P);                         // s_k x s_k
    CMat A = Gk.inverse() * M;                       // A_k
    A_.push_back(A);

    CMat C;                                          // C_k (s_{k-1} x s_k)
    if (step > 0) C = G_[step-1].inverse() * B_[step-1].adjoint() * Gk;
    else          C = CMat::Zero(0, s_k);
    C_.push_back(C);

    // residual R = D_W Q_k - Q_k A_k - Q_{k-1} C_k  (s_k columns)
    std::vector<Field> R = P;
    subtractCombine(R, Qk, A);
    if (step > 0) subtractCombine(R, Q_[step-1], C);

    if (reorthog)
      for (int j = 0; j <= step; j++) {
        CMat Mj = g5Inner(Q_[j], R);
        CMat Hj = G_[j].inverse() * Mj;
        subtractCombine(R, Q_[j], Hj);
      }
    if (deflate && !lockG_.empty()) for (auto& r : R) deflateVec(r);

    // --- look-ahead: augment S with D_W(neutral dir) until no serious breakdown ---
    const RealD relEps = degenRel();
    std::vector<Field> S = R;                        // working augmented column set
    const int maxLA = 3;
    for (int la = 0; ; la++) {
      CMat Gamma = g5Inner(S, S);
      Eigen::SelfAdjointEigenSolver<CMat> es(Gamma);
      Eigen::VectorXd D = es.eigenvalues();
      CMat U = es.eigenvectors();
      RealD dmax = D.cwiseAbs().maxCoeff();
      std::vector<int> serious;
      for (int i = 0; i < D.size(); i++) {
        if (std::abs(D(i)) >= relEps * dmax) continue;   // non-degenerate
        Field ri = combineCol(S, U.col(i));
        if (std::sqrt(norm2(ri)) < tol_) continue;        // happy (zero) -> will be dropped
        serious.push_back(i);                             // neutral but non-zero -> serious
      }
      if (serious.empty()) break;                         // no serious breakdown
      if (la >= maxLA) { log("look-ahead exhausted at step "+std::to_string(step)); break; }
      // augment with D_W of each serious-neutral direction (escapes the neutral cone)
      // Collect the new directions from the CURRENT S (U matches this S); append
      // only AFTER the loop -- growing S mid-loop would desync U/combineCol sizes.
      std::vector<Field> newdirs;
      for (int i : serious) {
        Field ri = combineCol(S, U.col(i));
        RealD rin = std::sqrt(norm2(ri));
        if (rin < tol_) continue;
        ri = ri * (1.0 / rin);
        Field dri(Grid_); Linop_.Op(ri, dri);
        // gamma5-orthogonalise the new direction against all previous blocks
        for (int j = 0; j <= step; j++) {
          int sj = (int)Q_[j].size();
          Field g5d(Grid_); g5_(dri, g5d);
          CVec proj(sj);
          for (int c = 0; c < sj; c++) proj(c) = toStd(innerProduct(Q_[j][c], g5d));
          CVec coef = G_[j].inverse() * proj;
          for (int c = 0; c < sj; c++) dri = dri - Q_[j][c] * coef(c);
        }
        RealD dn = std::sqrt(norm2(dri));
        if (dn < tol_) continue;                 // D_W(neutral) already in span -> skip
        newdirs.push_back(dri * (1.0 / dn));     // normalised
      }
      int added = (int)newdirs.size();
      for (auto& d : newdirs) S.push_back(d);
      if (added == 0) break;                      // nothing new to add -> stop expanding
      lookaheadCount_++;
      log("look-ahead at step "+std::to_string(step)+": block grown to "+std::to_string((int)S.size()));
    }

    // --- build Q_{k+1} from the (possibly augmented) set S ---
    CMat GammaS = g5Inner(S, S);
    Eigen::SelfAdjointEigenSolver<CMat> es(GammaS);
    Eigen::VectorXd D = es.eigenvalues();
    CMat U = es.eigenvectors();
    RealD dmax = D.cwiseAbs().maxCoeff();
    if (dmax < tol_ * tol_) { log("happy breakdown at step "+std::to_string(step)); return false; }
    // keep non-degenerate directions; absolute floor avoids dividing by ~0
    const RealD dfloor = std::max(relEps * dmax, tol_ * tol_);
    std::vector<int> keep;
    for (int i = 0; i < D.size(); i++) if (std::abs(D(i)) >= dfloor) keep.push_back(i);
    if (keep.empty()) { log("happy breakdown at step "+std::to_string(step)); return false; }

    // diagnostic: residual-Gram condition number kappa(Gamma_k) (oblique proj. blow-up)
    { double dmx = 0, dmn = 1e300;
      for (int i : keep) { double a = std::abs(D(i)); dmx = std::max(dmx,a); dmn = std::min(dmn,a); }
      kappaGamma_.push_back(dmn > 0 ? dmx/dmn : 1e300); }

    int s_kp1 = keep.size();
    std::vector<Field> Qkp1; Qkp1.reserve(s_kp1);
    CMat Gkp1 = CMat::Zero(s_kp1, s_kp1);
    for (int a = 0; a < s_kp1; a++) {
      int i = keep[a];
      Field q = combineCol(S, U.col(i)) * (1.0 / std::sqrt(std::abs(D(i))));
      Qkp1.push_back(q);
      Gkp1(a, a) = Cd(D(i) > 0 ? 1.0 : -1.0, 0.0);
    }
    // B_{k+1} = G_{k+1}^{-1} (Q_{k+1}^dag g5 R)   (s_{k+1} x s_k);  R = Q_{k+1} B_{k+1}
    CMat QtR = g5Inner(Qkp1, R);
    CMat Bkp1 = Gkp1.inverse() * QtR;

    G_.push_back(Gkp1); B_.push_back(Bkp1); Q_.push_back(Qkp1);

    // diagnostic: loss of gamma5-orthogonality, ||Q_1^dag g5 Q_{k+1}||_F (~0 ideally)
    { CMat e = g5Inner(Q_[0], Qkp1); etaLoss_.push_back(e.norm()); }
    return true;
  }

  // --- Ritz extraction (variable block-tridiagonal T_m) ---
  void computeRitzPairs(int m, Gamma5RitzSort sort) {
    std::vector<int> off(m + 1, 0);
    for (int k = 0; k < m; k++) off[k+1] = off[k] + sz(k);
    int dim = off[m];

    CMat Tm = CMat::Zero(dim, dim);
    for (int k = 0; k < m; k++) {
      Tm.block(off[k], off[k], sz(k), sz(k)) = A_[k];
      if (k < m - 1) {
        Tm.block(off[k+1], off[k],   sz(k+1), sz(k)) = B_[k];     // sub
        Tm.block(off[k],   off[k+1], sz(k), sz(k+1)) = C_[k+1];   // super
      }
    }

    Eigen::ComplexEigenSolver<CMat> ces(Tm);
    CVec lam = ces.eigenvalues();
    CMat Y   = ces.eigenvectors();

    std::vector<int> idx(dim);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b){ return ritzLess(lam(a), lam(b), sort); });

    const CMat& Bm1 = B_[m-1];      // s_m x s_{m-1}
    int s_last = sz(m-1);
    evals_.resize(dim); evecs_.clear(); residuals_.clear();
    evecs_.reserve(dim); residuals_.reserve(dim);

    for (int ji = 0; ji < dim; ji++) {
      int j = idx[ji];
      evals_(ji) = lam(j);
      CVec yj = Y.col(j);

      Field uj(Grid_); uj = Zero();
      for (int k = 0; k < m; k++)
        for (int c = 0; c < sz(k); c++)
          uj = uj + Q_[k][c] * yj(off[k] + c);
      evecs_.push_back(uj);

      // residual r_j = Q_{m} B_{m-1} tau_j,  tau_j = last-block entries of y_j
      CVec tau(s_last);
      for (int c = 0; c < s_last; c++) tau(c) = yj(off[m-1] + c);
      CVec Bt = Bm1 * tau;                  // length s_m
      Field rj(Grid_); rj = Zero();
      for (int c = 0; c < (int)Q_[m].size(); c++) rj = rj + Q_[m][c] * Bt(c);
      residuals_.push_back(std::sqrt(norm2(rj)));
    }
  }

  // --- locking / deflation (converged conjugate pairs) ---
  int lockConvergedPairs(int nWantedPairs) {
    int dim = (int)evecs_.size();
    int locked = 0;
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
      if (best < 0) continue;
      if (isAlreadyLocked(li)) { used[i] = used[best] = true; continue; }
      Field w0 = evecs_[i], w1 = evecs_[best];
      CMat2 G; std::vector<Field> blk = {w0, w1};
      CMat Gg = g5Inner(blk, blk);
      G << Gg(0,0), Gg(0,1), Gg(1,0), Gg(1,1);
      if (std::abs(G.determinant()) < 1e-12) { used[i] = used[best] = true; continue; }
      lockV_.push_back(w0); lockV_.push_back(w1);
      lockG_.push_back(G);
      lockEval_.push_back(Cd(real(evals_(i)),    imag(evals_(i))));
      lockEval_.push_back(Cd(real(evals_(best)), imag(evals_(best))));
      used[i] = used[best] = true; locked++;
    }
    return locked;
  }

  bool isAlreadyLocked(const Cd& lam) const {
    for (size_t b = 0; b < lockEval_.size(); b++)
      if (std::abs(lockEval_[b] - lam) < 10.0 * tol_) return true;
    return false;
  }

  // x <- x - sum_b V_b G_b^{-1} (V_b^dag g5 x)
  void deflateVec(Field& x) {
    for (size_t b = 0; b < lockG_.size(); b++) {
      const Field& w0 = lockV_[2*b]; const Field& w1 = lockV_[2*b+1];
      Field g5x(Grid_); g5_(x, g5x);
      Eigen::Vector2cd c;
      c(0) = toStd(innerProduct(w0, g5x));
      c(1) = toStd(innerProduct(w1, g5x));
      Eigen::Vector2cd a = lockG_[b].inverse() * c;
      x = x - (w0 * a(0) + w1 * a(1));
    }
  }

  void spliceLockedToFront() {
    if (lockG_.empty()) return;
    int nl = 2 * (int)lockG_.size();
    CVec ev(nl + evals_.size());
    std::vector<Field> ec; ec.reserve(nl + evecs_.size());
    std::vector<RealD> rs; rs.reserve(nl + residuals_.size());
    for (size_t b = 0; b < lockG_.size(); b++) {
      ev(2*b)   = ComplexD(lockEval_[2*b].real(),   lockEval_[2*b].imag());
      ev(2*b+1) = ComplexD(lockEval_[2*b+1].real(), lockEval_[2*b+1].imag());
      ec.push_back(lockV_[2*b]); ec.push_back(lockV_[2*b+1]);
      rs.push_back(0.0); rs.push_back(0.0);
    }
    for (int i = 0; i < (int)evals_.size(); i++) ev(nl + i) = evals_(i);
    for (auto& v : evecs_)    ec.push_back(v);
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
