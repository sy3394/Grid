/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_gamma5_block_lanczos.cc

    gamma5-Block Lanczos driver for the Wilson Dirac operator D_W.

    Computes complex eigenvalues of D_W directly (not H_W = g5 D_W) via the
    gamma5-Block Lanczos algorithm (Grid/algorithms/iterative/Gamma5BlockLanczos.h):
    a gamma5-metric block recurrence (cheap, conjugate-pair targeting) with refined
    Euclidean Ritz-vector extraction and raw-residual convergence.

    Modes:
      DIRECT (default)             : Lanczos on D_W (peripheral eigenvalues).
      SHIFT-INVERT (--shift sigma) : Lanczos on (D_W - sigma)^{-1}; converges the
                                     eigenvalues of D_W nearest real sigma.  The
                                     inner solve uses normal-equations CG; a Ritz
                                     value theta maps back as lambda = sigma + 1/theta.
      SWEEP (--shift-sweep lo:hi:n): n shifts in [lo,hi] in one job; accumulate and
                                     de-duplicate the converged modes of each window.

    Diagnostics / comparison:
      --compare : head-to-head against plain Euclidean Arnoldi on the SAME
                  shift-inverted operator (same inner CG), equal Krylov dim.
      --history : residual & #converged vs Krylov dim for g5bl and Arnoldi
                  (-> <out>.g5bl.hist, <out>.arnoldi.hist), Euclidean D_W residual.
      --check   : verify each output eigenpair against the raw operator.
      --dense   : exact dense diagonalisation of D_W (small lattices).

    Test gauges:  --cold (free field) | --weak eps (slight perturbation) |
                  --config PATH (NERSC) | (none -> hot random).

    Cost reporting: Krylov dim = number of (D_W - sigma)^{-1} applications;
    total CG iters = inner-solve work (the absolute D_W-matvec cost).

*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/Gamma5BlockLanczos.h>
#include <Grid/algorithms/iterative/RefinedArnoldi.h>
#include <functional>

using namespace std;
using namespace Grid;

typedef WilsonFermionD                         WilsonOp;
typedef typename WilsonFermionD::FermionField  FermionField;

// ---- shift-invert operator:  Op(in) = (D_W - sigma)^{-1} in  via normal-eq CG ----
template<class Matrix, class Field>
class ShiftInvertNE : public LinearOperatorBase<Field> {
  Matrix&                           D_;
  MdagMLinearOperator<Matrix,Field> MdagM_;
  ConjugateGradient<Field>&         cg_;
public:
  long nApply = 0;   // number of (D_W-sigma)^{-1} applications (= Krylov steps)
  long nCG    = 0;   // total inner CG iterations (absolute cost)
  ShiftInvertNE(Matrix& D, ConjugateGradient<Field>& cg) : D_(D), MdagM_(D), cg_(cg) {}
  void Op(const Field& in, Field& out) {
    Field Mdb(in.Grid()); D_.Mdag(in, Mdb); out = Zero();
    cg_(MdagM_, Mdb, out); nApply++; nCG += cg_.IterationsToComplete;
  }
  void AdjOp(const Field&, Field&)                       { assert(0); }
  void OpDiag(const Field&, Field&)                      { assert(0); }
  void OpDir(const Field&, Field&, int, int)             { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)       { assert(0); }
  void HermOp(const Field&, Field&)                      { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&){ assert(0); }
};

// ---- exact dense diagonalisation of D_W (small lattices only) ----
template<class Field>
static void denseDiag(LinearOperatorBase<Field>& Op, GridCartesian* grid,
                      const std::string& fname, const std::string& tag) {
  int V = 1; for (int d = 0; d < Nd; d++) V *= grid->FullDimensions()[d];
  int N = Ns * Nc * V;
  std::cout << GridLogMessage << "denseDiag: " << N << "x" << N << " (" << N << " matvecs)" << std::endl;
  Eigen::MatrixXcd M(N, N);
  Field ek(grid), Dek(grid);
  typedef typename Field::scalar_object sobj;
  std::vector<sobj> col;
  for (int k = 0; k < N; k++) {
    int site = k / (Ns*Nc), sc = k % (Ns*Nc), s = sc / Nc, c = sc % Nc;
    Coordinate coor(Nd); Lexicographic::CoorFromIndex(coor, site, grid->FullDimensions());
    ek = Zero(); sobj o; o = Zero(); o()(s)(c) = Complex(1.0, 0.0); pokeSite(o, ek, coor);
    Op.Op(ek, Dek); unvectorizeToLexOrdArray(col, Dek);
    for (int j = 0; j < V; j++) for (int sj = 0; sj < Ns; sj++) for (int cc = 0; cc < Nc; cc++)
      M(j*(Ns*Nc) + sj*Nc + cc, k) = static_cast<std::complex<double>>(col[j]()(sj)(cc));
  }
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(M, false);
  Eigen::VectorXcd lam = es.eigenvalues();
  std::vector<std::complex<double>> v(lam.data(), lam.data() + lam.size());
  std::sort(v.begin(), v.end(), [](auto a, auto b){ return std::abs(a.imag()) < std::abs(b.imag()); });
  std::ofstream f(fname); f << std::setprecision(10);
  for (auto& z : v) f << tag << "  " << z.real() << "  " << z.imag() << "\n";
  std::cout << GridLogMessage << "denseDiag -> " << fname << std::endl;
}

// ---- plain Euclidean Arnoldi on Op (baseline); returns V and H to budget m ----
template<class Field>
static void arnoldi(LinearOperatorBase<Field>& Op, GridBase* grid, const Field& v0,
                    int m, std::vector<Field>& V, Eigen::MatrixXcd& H, int& mdone) {
  Field v(grid); v = v0; v = v * (1.0/std::sqrt(norm2(v))); V.push_back(v);
  H = Eigen::MatrixXcd::Zero(m + 1, m); mdone = m;
  for (int j = 0; j < m; j++) {
    Field w(grid); Op.Op(V[j], w);
    for (int i = 0; i <= j; i++) {
      auto h = innerProduct(V[i], w);
      H(i, j) = std::complex<double>((double)real(h), (double)imag(h));
      w = w - V[i] * h;
    }
    double hn = std::sqrt(norm2(w)); H(j + 1, j) = hn;
    if (hn < 1e-12) { mdone = j + 1; break; }
    if (j + 1 < m) V.push_back(w * (1.0/hn));
  }
}

// ---- Block Arnoldi (block size 2) seeded with [v, g5 v]; Euclidean-orthonormal ----
// Spans the SAME block Krylov subspace as g5bl, but in the Euclidean metric (so the
// basis is perfectly conditioned and the Ritz extraction has no oblique penalty).
// Returns the flat column basis V (p*(nblk+1) columns) and the block upper-Hessenberg
// H of size p*(nblk+1) x p*nblk.  Orthogonalisation is done twice for stability.
template<class Field>
static void blockArnoldiG5(LinearOperatorBase<Field>& Op, GridBase* grid, const Field& v0,
                           std::function<void(const Field&, Field&)> g5, int nblk,
                           std::vector<Field>& V, Eigen::MatrixXcd& H) {
  const int p = 2;
  auto cd = [](ComplexD z){ return std::complex<double>((double)real(z), (double)imag(z)); };
  GridParallelRNG rng(grid); rng.SeedFixedIntegers({9,8,7,6});
  V.clear();
  Field a(grid), b(grid), w(grid);
  a = v0; a = a * (1.0/std::sqrt(norm2(a)));
  g5(v0, b); { ComplexD h = innerProduct(a, b); b = b - a*h; }   // [v, g5 v], orthonormal
  b = b * (1.0/std::sqrt(norm2(b)));
  V.push_back(a); V.push_back(b);
  H = Eigen::MatrixXcd::Zero(p*(nblk+1), p*nblk);
  for (int j = 0; j < nblk; j++) {
    for (int c = 0; c < p; c++) {
      Op.Op(V[p*j + c], w);
      int prev = (int)V.size();                  // all columns already orthonormal
      for (int pass = 0; pass < 2; pass++)
        for (int i = 0; i < prev; i++) {
          ComplexD h = innerProduct(V[i], w);
          H(i, p*j + c) += cd(h); w = w - V[i]*h;
        }
      double nrm = std::sqrt(norm2(w));
      if (nrm < 1e-10) {                          // breakdown: continue with a fresh dir
        gaussian(rng, w);
        for (int i = 0; i < prev; i++) { ComplexD h = innerProduct(V[i], w); w = w - V[i]*h; }
        nrm = std::sqrt(norm2(w));
      }
      H(prev, p*j + c) = nrm;
      V.push_back(w * (1.0/nrm));
    }
  }
}

// ---- Refined Ritz extraction from an (block) Arnoldi factorisation A V_m = V_{m+1} Hbar.
// For each Ritz value theta of H[0:pm,0:pm], the refined vector minimises the Euclidean
// residual ||A u - theta u|| over span(V_m): the smallest right singular vector of
// (Hbar - theta [I;0]).  lamOf maps theta back to the D_W eigenvalue. ----
template<class Field>
static void refinedRitz(const std::vector<Field>& V, const Eigen::MatrixXcd& H, int p, int m,
                        GridBase* grid, std::function<std::complex<double>(std::complex<double>)> lamOf,
                        std::vector<std::complex<double>>& L, std::vector<Field>& U) {
  int pm = p*m;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, pm, pm), false);
  Eigen::MatrixXcd Hbar = H.block(0, 0, p*(m+1), pm);
  for (int j = 0; j < pm; j++) {
    std::complex<double> th = es.eigenvalues()(j);
    Eigen::MatrixXcd C = Hbar; for (int i = 0; i < pm; i++) C(i, i) -= th;
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(C, Eigen::ComputeThinV);
    Eigen::VectorXcd y = svd.matrixV().col(pm - 1);
    Field u(grid); u = Zero(); for (int k = 0; k < pm; k++) u = u + V[k]*y(k);
    L.push_back(lamOf(th)); U.push_back(u);
  }
}

static std::string getOpt(int c, char** v, const std::string& k, const std::string& d="") {
  for (int i = 1; i < c-1; i++) if (k == v[i]) return v[i+1]; return d;
}
static bool hasOpt(int c, char** v, const std::string& k) {
  for (int i = 1; i < c; i++) if (k == v[i]) return true; return false;
}
static std::vector<double> parseSweep(const std::string& s) {
  std::vector<double> r; double lo, hi; int n;
  if (sscanf(s.c_str(), "%lf:%lf:%d", &lo, &hi, &n) == 3 && n >= 1) {
    if (n == 1) r.push_back(lo);
    else for (int i = 0; i < n; i++) r.push_back(lo + (hi - lo) * i / (n - 1));
  }
  return r;
}

typedef std::pair<std::complex<double>, double> EvalRes;  // (lambda, residual)
static void addDedup(std::vector<EvalRes>& C, std::complex<double> l, double r, double eps) {
  for (auto& e : C) if (std::abs(e.first - l) < eps) { if (r < e.second) e = {l, r}; return; }
  C.push_back({l, r});
}
static void writeEvals(const std::string& fn, const std::string& tag, std::vector<EvalRes> v) {
  std::sort(v.begin(), v.end(), [](const EvalRes& a, const EvalRes& b){
    return std::abs(a.first.imag()) < std::abs(b.first.imag()); });
  std::ofstream f(fn); f << std::setprecision(10);
  for (auto& e : v) f << tag << "  " << e.first.real() << "  " << e.first.imag() << "\n";
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  RealD mass   = std::stod(getOpt(argc, argv, "--mass",   "-0.5"));
  int   steps  = std::stoi(getOpt(argc, argv, "--steps",  "40"));
  int   wanted = std::stoi(getOpt(argc, argv, "--wanted", "10"));
  int   cycles = std::stoi(getOpt(argc, argv, "--cycles", "8"));
  RealD tol    = std::stod(getOpt(argc, argv, "--tol",    "1e-8"));
  RealD stol   = std::stod(getOpt(argc, argv, "--stol",   "1e-9"));
  int   siter  = std::stoi(getOpt(argc, argv, "--siter",  "20000"));
  RealD accept = std::stod(getOpt(argc, argv, "--accept", "1e-3"));
  RealD dedupe = std::stod(getOpt(argc, argv, "--dedupe", "1e-4"));
  double weak  = hasOpt(argc, argv, "--weak") ? std::stod(getOpt(argc, argv, "--weak", "0.1")) : 0.0;
  std::string cfg = getOpt(argc, argv, "--config", "");
  std::string out = getOpt(argc, argv, "--out",    "g5bl_evals.dat");
  std::string tag = getOpt(argc, argv, "--tag",    std::to_string(mass));
  bool reorth  = !hasOpt(argc, argv, "--noreorth");
  bool isoOnly = hasOpt(argc, argv, "--isolation-only");
  bool cold    = hasOpt(argc, argv, "--cold");
  bool check   = hasOpt(argc, argv, "--check");

  std::vector<double> sigmas;
  if      (hasOpt(argc, argv, "--shift-sweep")) sigmas = parseSweep(getOpt(argc, argv, "--shift-sweep", ""));
  else if (hasOpt(argc, argv, "--shift"))       sigmas = { std::stod(getOpt(argc, argv, "--shift", "0")) };
  bool shiftMode = !sigmas.empty();

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::cout << GridLogMessage << "g5bl: lattice=" << GridDefaultLatt() << " mass=" << mass
            << (shiftMode ? "  shift-invert" : "  direct") << std::endl;

  // ---- gauge ----
  LatticeGaugeField Umu(UGrid);
  if (cold) { SU<Nc>::ColdConfiguration(Umu); }
  else if (weak > 0.0) {
    GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4});
    LatticeColourMatrix Ul(UGrid);
    for (int mu = 0; mu < Nd; mu++) { SU<Nc>::LieRandomize(pRNG, Ul, weak); PokeIndex<LorentzIndex>(Umu, Ul, mu); }
  } else if (cfg.empty()) {
    GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4}); SU<Nc>::HotConfiguration(pRNG, Umu);
  } else {
    FieldMetaData h; NerscIO::readConfiguration(Umu, h, cfg);
    std::cout << GridLogMessage << "config " << cfg << " plaquette=" << h.plaquette << std::endl;
  }

  Gamma G5(Gamma::Algebra::Gamma5);
  auto gamma5 = [&G5](const FermionField& in, FermionField& out){ out = G5 * in; };
  std::vector<Complex> bc = {1,1,1,-1};
  WilsonOp::ImplParams wpar(bc);
  WilsonOp Dw(Umu, *UGrid, *UrbGrid, mass, wpar);
  NonHermitianLinearOperator<WilsonOp, FermionField> DW(Dw);

  GridParallelRNG RNG(UGrid); RNG.SeedFixedIntegers({5,6,7,8});
  FermionField v0(UGrid), v1(UGrid); random(RNG, v0); random(RNG, v1);

  if (hasOpt(argc, argv, "--dense")) denseDiag(DW, UGrid, out + ".dense", tag);

  // raw residual of (u,lambda) against D_W
  FermionField wbuf(UGrid);
  auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
    DW.Op(u, wbuf); ComplexD lf(lam.real(), lam.imag());
    FermionField t(UGrid); t = wbuf - u * lf; return std::sqrt(norm2(t)/norm2(u)); };

  // ================= DIRECT mode (eigenvalue computation) =================
  if (!shiftMode && !hasOpt(argc, argv, "--history")) {
    auto report = [&](Gamma5BlockLanczos<FermionField>& g, const std::string& fn){
      const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
      std::vector<EvalRes> v;
      for (int i = 0; i < (int)ev.size(); i++) v.push_back({ev(i), rs[i]});
      std::sort(v.begin(), v.end(), [](const EvalRes& a, const EvalRes& b){
        return std::abs(a.first.imag()) < std::abs(b.first.imag()); });
      for (int i = 0; i < std::min((int)v.size(), 2*wanted); i++)
        std::cout << GridLogMessage << "  [" << std::setw(3) << i << "] lambda=("
                  << v[i].first.real() << "," << v[i].first.imag() << ")  res=" << v[i].second << std::endl;
      writeEvals(fn, tag, v);
    };
    Gamma5BlockLanczos<FermionField> g(DW, UGrid, gamma5, tol, 1);
    g.setRawCheck(&DW, 0.0, /*shiftInvert=*/false);
    if (isoOnly) { std::cout << GridLogMessage << "-- isolation --" << std::endl; g(v0, v1, steps, reorth); report(g, out); }
    else         { std::cout << GridLogMessage << "-- thick restart --" << std::endl;
                   g.thickRestart(v0, v1, cycles, steps, wanted, reorth); report(g, out); }
    std::cout << GridLogMessage << "Done." << std::endl; Grid_finalize(); return 0;
  }

  // ===== --history: internal residual vs RAW Euclidean residual, vs Krylov dim =====
  // Works in DIRECT (M=D_W, lambda=mu) and SHIFT-INVERT (M=(D_W-sigma)^-1,
  // lambda=sigma+1/mu) regimes.  Columns:
  //   krylov_dim  min_internal_res  min_raw_res  n_raw<1e-2  n_raw<1e-4
  // internal_res = the residual the method reports natively (g5bl: gamma5-Galerkin
  // ||Q_{m+1}B_{m+1}tau||; Arnoldi: Hessenberg estimate).  raw_res = honest
  // ||D_W u - lambda u||/||u||.  A fast-dropping internal residual with a lagging
  // raw residual = optimistic convergence.
  if (hasOpt(argc, argv, "--history")) {
    bool shift = shiftMode; double sg = shift ? sigmas[0] : 0.0;
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, shift ? mass - sg : mass, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    LinearOperatorBase<FermionField>& M = shift
        ? static_cast<LinearOperatorBase<FermionField>&>(SI)
        : static_cast<LinearOperatorBase<FermionField>&>(DW);
    auto lamOf = [&](std::complex<double> mu){ return shift ? (sg + 1.0/mu) : mu; };
    auto rawStats = [&](const std::vector<std::complex<double>>& L, const std::vector<FermionField>& U){
      double mn = 1e30; int b2 = 0, b4 = 0;
      for (int i = 0; i < (int)L.size(); i++) { double r = rawRes(U[i], L[i]); mn = std::min(mn,r); if(r<1e-2)b2++; if(r<1e-4)b4++; }
      return std::array<double,3>{mn,(double)b2,(double)b4}; };
    Gamma5RitzSort sort = shift ? G5SortAbsDescending : G5SortAbsImagAscending;

    Gamma5BlockLanczos<FermionField> g(M, UGrid, gamma5, tol, 0);
    if (shift) g.setRawCheck(&DW, sg, true);
    g(v0, v1, steps, reorth, sort);
    // Columns: galerkin_res = ||Q_{m+1}B_{m+1}tau|| of the standard V_m y vector (the
    // metric the original code reported); galerkin_raw = the RAW residual of that
    // SAME standard vector; refined_raw = the RAW residual of the refined vector.
    std::ofstream fg(out + ".g5bl.hist");
    fg << "# krylov_dim min_galerkin_res min_galerkin_raw min_refined_raw n_refined_1e-2 n_refined_1e-4\n";
    for (int m = 1; m <= g.getNumSteps(); m++) {
      g.extractRitzAt(m, sort);
      const auto& ev = g.getEvals(); const auto& uv = g.getEvecs();
      const auto& gr = g.getGalerkinResiduals(); const auto& grr = g.getGalerkinRawResiduals();
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      double gmin = 1e30, grmin = 1e30;
      for (int i = 0; i < (int)ev.size(); i++) {
        L.push_back(lamOf(ev(i))); U.push_back(uv[i]);
        if (i < (int)gr.size())  gmin  = std::min(gmin,  (double)gr[i]);
        if (i < (int)grr.size()) grmin = std::min(grmin, (double)grr[i]);
      }
      auto s = rawStats(L, U);
      fg << (int)ev.size() << " " << gmin << " " << grmin << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    std::cout << GridLogMessage << "g5bl: " << g.getNumSteps() << " steps" << std::endl;

    std::vector<FermionField> Va; Eigen::MatrixXcd H; int adone;
    arnoldi(M, UGrid, v0, 2*steps, Va, H, adone);
    std::ofstream fa(out + ".arnoldi.hist");
    fa << "# krylov_dim min_internal_res min_raw_res n_raw_1e-2 n_raw_1e-4\n";
    for (int m = 1; m <= adone; m++) {
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, m, m));
      auto lam = es.eigenvalues(); auto Y = es.eigenvectors();
      double hsub = (m < (int)H.rows()) ? std::abs(H(m, m-1)) : 0.0, gmin = 1e30;
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      for (int j = 0; j < m; j++) {
        L.push_back(lamOf(es.eigenvalues()(j)));
        FermionField u(UGrid); u = Zero();
        for (int k = 0; k < m && k < (int)Va.size(); k++) u = u + Va[k] * Y(k, j);
        U.push_back(u);
        gmin = std::min(gmin, hsub * std::abs(Y(m-1, j)));   // Arnoldi internal estimate
      }
      auto s = rawStats(L, U);
      fa << m << " " << gmin << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    std::cout << GridLogMessage << "arnoldi (single-vector): " << adone << " steps" << std::endl;

    std::function<std::complex<double>(std::complex<double>)> lamFn =
        [&](std::complex<double> mu){ return lamOf(mu); };

    // single-vector Arnoldi with REFINED extraction: isolates the refined-extraction
    // effect from the block-seeding effect (same basis as arnoldi.hist above).
    std::ofstream fr(out + ".arnoldiref.hist");
    fr << "# krylov_dim min_refined_raw n_refined_1e-2 n_refined_1e-4\n";
    for (int m = 1; m <= adone; m++) {
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      refinedRitz<FermionField>(Va, H, 1, m, UGrid, lamFn, L, U);
      auto s = rawStats(L, U);
      fr << m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }

    // block Arnoldi seeded with [v, g5 v], refined extraction -- the real competitor.
    std::vector<FermionField> Vb; Eigen::MatrixXcd Hb;
    blockArnoldiG5<FermionField>(M, UGrid, v0, gamma5, steps, Vb, Hb);
    std::ofstream fb(out + ".blockarnoldi.hist");
    fb << "# krylov_dim min_refined_raw n_refined_1e-2 n_refined_1e-4\n";
    for (int m = 1; m <= steps; m++) {
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      refinedRitz<FermionField>(Vb, Hb, 2, m, UGrid, lamFn, L, U);
      auto s = rawStats(L, U);
      fb << 2*m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    std::cout << GridLogMessage << "histories -> " << out
              << ".{g5bl,arnoldi,arnoldiref,blockarnoldi}.hist" << std::endl;
    Grid_finalize(); return 0;
  }

  // ================= --compare (head-to-head at one Krylov budget) =================
  if (hasOpt(argc, argv, "--compare")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    int budget = 2 * steps;
    // count DISTINCT converged modes (dedup nearby eigenvalues so duplicate refined
    // Ritz vectors don't inflate the breadth count).
    auto nconv = [&](const std::vector<std::complex<double>>& L, const std::vector<FermionField>& U){
      std::vector<EvalRes> C;
      for (int i = 0; i < (int)L.size(); i++) { double r = rawRes(U[i], L[i]); if (r < accept) addDedup(C, L[i], r, dedupe); }
      return (int)C.size(); };

    SI.nApply = SI.nCG = 0;
    Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 0);
    g.setRawCheck(&DW, sg, true);
    GridStopWatch t1; t1.Start(); g(v0, v1, steps, reorth, G5SortAbsDescending); t1.Stop();
    std::vector<std::complex<double>> Lg; std::vector<FermionField> Ug;
    { const auto& ev = g.getEvals(); const auto& uv = g.getEvecs();
      for (int i = 0; i < (int)ev.size(); i++) { Lg.push_back(sg + 1.0/ev(i)); Ug.push_back(uv[i]); } }
    long g_app = SI.nApply, g_cg = SI.nCG; int g_c = nconv(Lg, Ug);

    SI.nApply = SI.nCG = 0;
    std::vector<FermionField> Va; Eigen::MatrixXcd H; int adone;
    GridStopWatch t2; t2.Start(); arnoldi(SI, UGrid, v0, budget, Va, H, adone); t2.Stop();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, adone, adone));
    std::vector<std::complex<double>> La; std::vector<FermionField> Ua;
    for (int j = 0; j < adone; j++) {
      La.push_back(sg + 1.0/es.eigenvalues()(j));
      FermionField u(UGrid); u = Zero();
      for (int k = 0; k < adone && k < (int)Va.size(); k++) u = u + Va[k] * es.eigenvectors()(k, j);
      Ua.push_back(u);
    }
    long a_app = SI.nApply, a_cg = SI.nCG; int a_c = nconv(La, Ua);

    std::function<std::complex<double>(std::complex<double>)> lamFn =
        [&](std::complex<double> mu){ return sg + 1.0/mu; };
    std::vector<std::complex<double>> Lar; std::vector<FermionField> Uar;
    refinedRitz<FermionField>(Va, H, 1, adone, UGrid, lamFn, Lar, Uar);   // single-vec REFINED
    int ar_c = nconv(Lar, Uar);

    SI.nApply = SI.nCG = 0;
    std::vector<FermionField> Vb; Eigen::MatrixXcd Hb;
    GridStopWatch t3; t3.Start();
    blockArnoldiG5<FermionField>(SI, UGrid, v0, gamma5, steps, Vb, Hb);
    std::vector<std::complex<double>> Lb; std::vector<FermionField> Ub;
    refinedRitz<FermionField>(Vb, Hb, 2, steps, UGrid, lamFn, Lb, Ub);
    t3.Stop();
    long b_app = SI.nApply, b_cg = SI.nCG; int b_c = nconv(Lb, Ub);

    std::cout << GridLogMessage << "HEAD-TO-HEAD sigma=" << sg << " accept=" << accept
              << " (Krylov dim = 2*steps = " << budget << ")" << std::endl;
    std::cout << GridLogMessage << std::setw(20) << "method" << std::setw(12) << "krylov_dim"
              << std::setw(12) << "CG_iters" << std::setw(10) << "time(s)" << std::setw(8) << "conv" << std::endl;
    auto row = [&](const std::string& nm, long app, long cg, double tm, int cv){
      std::cout << GridLogMessage << std::setw(20) << nm << std::setw(12) << app << std::setw(12) << cg
                << std::setw(10) << tm << std::setw(8) << cv << std::endl; };
    std::cout << GridLogMessage << "(conv = DISTINCT modes with raw residual < accept)" << std::endl;
    row("g5bl (refined)",        g_app, g_cg, t1.useconds()*1e-6, g_c);
    row("arnoldi (1-vec std)",   a_app, a_cg, t2.useconds()*1e-6, a_c);
    row("arnoldi (1-vec refed)", a_app, a_cg, t2.useconds()*1e-6, ar_c);
    row("blockArnoldi [v,g5v]",  b_app, b_cg, t3.useconds()*1e-6, b_c);
    Grid_finalize(); return 0;
  }

  // ================= SHIFT-INVERT (single or sweep): the eigenvalue computation =================
  // Default solver: RefinedArnoldi (single-vector Arnoldi + refined extraction) -- best
  // for interior D_W modes.  --g5bl opts back into gamma5-block Lanczos (thick restart).
  bool useG5bl = hasOpt(argc, argv, "--g5bl");
  int  kdim    = std::stoi(getOpt(argc, argv, "--krylov", std::to_string(2*steps)));
  std::cout << GridLogMessage << "solver: " << (useG5bl ? "gamma5-block Lanczos" : "RefinedArnoldi")
            << (useG5bl ? "" : ("  (Krylov dim " + std::to_string(kdim) + ")")) << std::endl;
  std::vector<EvalRes> collected;
  long totApply = 0, totCG = 0; double tSolve = 0;
  for (size_t is = 0; is < sigmas.size(); is++) {
    double sg = sigmas[is];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    GridStopWatch sw; sw.Start();
    int nc = 0;
    if (useG5bl) {
      Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 1);
      g.setRawCheck(&DW, sg, true);
      if (isoOnly) g(v0, v1, steps, reorth, G5SortAbsDescending);
      else         g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsDescending);
      const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
      for (int i = 0; i < (int)ev.size(); i++)
        if (rs[i] < accept) { addDedup(collected, sg + 1.0/ev(i), rs[i], dedupe); nc++; }
    } else {
      RefinedArnoldi<FermionField> a(SI, UGrid, tol, 1);
      a.setRawCheck(&DW, sg, true);
      a(v0, kdim, RASortAbsDescending);                 // nearest-sigma modes first
      const auto& ev = a.getEvals(); const auto& rs = a.getResiduals();
      for (int i = 0; i < (int)ev.size(); i++)
        if (rs[i] < accept) { addDedup(collected, sg + 1.0/ev(i), rs[i], dedupe); nc++; }
    }
    sw.Stop(); tSolve += sw.useconds()*1e-6; totApply += SI.nApply; totCG += SI.nCG;
    std::cout << GridLogMessage << "sigma=" << sg << ": " << nc << " converged; total " << collected.size() << std::endl;
  }
  std::cout << GridLogMessage << "Krylov dim (applications): " << totApply
            << "   total CG iters: " << totCG << "   time: " << tSolve << " s" << std::endl;
  std::cout << GridLogMessage << "collected " << collected.size() << " distinct eigenvalues -> " << out << std::endl;
  (void)check;  // residuals in `collected` are already the raw D_W Euclidean residual
  writeEvals(out, tag, collected);

  std::cout << GridLogMessage << "Done." << std::endl;
  Grid_finalize();
  return 0;
}
