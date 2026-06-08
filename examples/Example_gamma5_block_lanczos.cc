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

  // ================= DIRECT mode =================
  if (!shiftMode) {
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

  // ================= --history (g5bl vs Arnoldi, residual & #conv vs Krylov dim) =================
  if (hasOpt(argc, argv, "--history")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    auto stats = [&](const std::vector<std::complex<double>>& L, const std::vector<FermionField>& U){
      double mn = 1e30; int b2 = 0, b4 = 0;
      for (int i = 0; i < (int)L.size(); i++) { double r = rawRes(U[i], L[i]); mn = std::min(mn, r); if (r<1e-2) b2++; if (r<1e-4) b4++; }
      return std::array<double,3>{mn, (double)b2, (double)b4}; };

    Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 0);
    g.setRawCheck(&DW, sg, true);
    g(v0, v1, steps, reorth, G5SortAbsDescending);
    std::ofstream fg(out + ".g5bl.hist"); fg << "# krylov_dim min_raw_res n_below_1e-2 n_below_1e-4\n";
    for (int m = 1; m <= g.getNumSteps(); m++) {
      g.extractRitzAt(m, G5SortAbsDescending);
      const auto& ev = g.getEvals(); const auto& uv = g.getEvecs();
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      for (int i = 0; i < (int)ev.size(); i++) { L.push_back(sg + 1.0/ev(i)); U.push_back(uv[i]); }
      auto s = stats(L, U);
      fg << (int)ev.size() << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    std::cout << GridLogMessage << "g5bl: " << SI.nApply << " applications, " << SI.nCG << " total CG iters" << std::endl;

    SI.nApply = SI.nCG = 0;
    std::vector<FermionField> Va; Eigen::MatrixXcd H; int adone;
    arnoldi(SI, UGrid, v0, 2*steps, Va, H, adone);
    std::ofstream fa(out + ".arnoldi.hist"); fa << "# krylov_dim min_raw_res n_below_1e-2 n_below_1e-4\n";
    for (int m = 1; m <= adone; m++) {
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, m, m));
      auto lam = es.eigenvalues(); auto Y = es.eigenvectors();
      std::vector<std::complex<double>> L; std::vector<FermionField> U;
      for (int j = 0; j < m; j++) {
        L.push_back(sg + 1.0/lam(j));
        FermionField u(UGrid); u = Zero();
        for (int k = 0; k < m && k < (int)Va.size(); k++) u = u + Va[k] * Y(k, j);
        U.push_back(u);
      }
      auto s = stats(L, U);
      fa << m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    std::cout << GridLogMessage << "arnoldi: " << SI.nApply << " applications, " << SI.nCG << " total CG iters" << std::endl;
    std::cout << GridLogMessage << "histories -> " << out << ".{g5bl,arnoldi}.hist" << std::endl;
    Grid_finalize(); return 0;
  }

  // ================= --compare (head-to-head at one Krylov budget) =================
  if (hasOpt(argc, argv, "--compare")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    int budget = 2 * steps;
    auto nconv = [&](const std::vector<std::complex<double>>& L, const std::vector<FermionField>& U){
      int n = 0; for (int i = 0; i < (int)L.size(); i++) if (rawRes(U[i], L[i]) < accept) n++; return n; };

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

    std::cout << GridLogMessage << "HEAD-TO-HEAD sigma=" << sg << " accept=" << accept << std::endl;
    std::cout << GridLogMessage << std::setw(18) << "method" << std::setw(12) << "krylov_dim"
              << std::setw(12) << "CG_iters" << std::setw(10) << "time(s)" << std::setw(10) << "conv" << std::endl;
    std::cout << GridLogMessage << std::setw(18) << "g5bl" << std::setw(12) << g_app << std::setw(12) << g_cg
              << std::setw(10) << t1.useconds()*1e-6 << std::setw(10) << g_c << std::endl;
    std::cout << GridLogMessage << std::setw(18) << "arnoldi" << std::setw(12) << a_app << std::setw(12) << a_cg
              << std::setw(10) << t2.useconds()*1e-6 << std::setw(10) << a_c << std::endl;
    Grid_finalize(); return 0;
  }

  // ================= SHIFT-INVERT (single or sweep): the eigenvalue computation =================
  std::vector<EvalRes> collected;
  long totApply = 0, totCG = 0; double tSolve = 0;
  for (size_t is = 0; is < sigmas.size(); is++) {
    double sg = sigmas[is];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 1);
    g.setRawCheck(&DW, sg, true);
    GridStopWatch sw; sw.Start();
    if (isoOnly) g(v0, v1, steps, reorth, G5SortAbsDescending);
    else         g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsDescending);
    sw.Stop(); tSolve += sw.useconds()*1e-6; totApply += SI.nApply; totCG += SI.nCG;

    const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
    int nc = 0;
    for (int i = 0; i < (int)ev.size(); i++) {
      if (rs[i] >= accept) continue;
      addDedup(collected, sg + 1.0/ev(i), rs[i], dedupe); nc++;
    }
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
