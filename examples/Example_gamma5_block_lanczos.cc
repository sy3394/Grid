/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_gamma5_block_lanczos.cc

    Driver: interior eigenvalues of the Wilson Dirac operator D_W.

    Demonstrates two Grid eigensolvers on D_W (both compute the COMPLEX eigenvalues
    of D_W directly, not H_W = g5 D_W):
      RefinedArnoldi      (Grid/algorithms/iterative/RefinedArnoldi.h)      -- default
      Gamma5BlockLanczos  (Grid/algorithms/iterative/Gamma5BlockLanczos.h)  -- --g5bl

    Interior modes are reached by SHIFT-INVERT: run the solver on (D_W - sigma)^{-1}
    (assembled here from Grid's ConjugateGradient via normal equations); a Ritz value
    mu maps back as lambda = sigma + 1/mu.  Real sigma keeps (D_W - sigma) g5-Hermitian.

    Modes:
      --shift sigma         single shift-invert about real sigma
      --shift-sweep lo:hi:n  n shifts in [lo,hi] in one job (accumulate + de-duplicate)
      (neither)             direct Lanczos/Arnoldi on D_W (peripheral eigenvalues)
    Diagnostics:
      --compare             RefinedArnoldi vs Gamma5BlockLanczos at equal Krylov dim
      --history             residual & #converged vs Krylov dim, both solvers
      --dense               exact dense diagonalisation of D_W (small lattices)
    Gauge:  --cold | --weak eps | --config PATH (NERSC) | (none -> hot random)

*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/Gamma5BlockLanczos.h>
#include <Grid/algorithms/iterative/RefinedArnoldi.h>

using namespace Grid;

typedef WilsonFermionD                         WilsonOp;
typedef typename WilsonFermionD::FermionField  FermionField;

// ---- shift-invert operator  Op(in) = (D_W - sigma)^{-1} in  (normal-equations CG) ----
// Shows how to compose Grid's ConjugateGradient + MdagMLinearOperator into an operator.
template<class Matrix, class Field>
class ShiftInvertNE : public LinearOperatorBase<Field> {
  Matrix&                           D_;
  MdagMLinearOperator<Matrix,Field> MdagM_;
  ConjugateGradient<Field>&         cg_;
public:
  long nApply = 0;   // (D_W - sigma)^{-1} applications (= Krylov dim)
  long nCG    = 0;   // total inner CG iterations (absolute cost)
  ShiftInvertNE(Matrix& D, ConjugateGradient<Field>& cg) : D_(D), MdagM_(D), cg_(cg) {}
  void Op(const Field& in, Field& out) {
    Field Mdb(in.Grid()); D_.Mdag(in, Mdb); out = Zero();
    cg_(MdagM_, Mdb, out); nApply++; nCG += cg_.IterationsToComplete;
  }
  void AdjOp(const Field&, Field&)                        { assert(0); }
  void OpDiag(const Field&, Field&)                       { assert(0); }
  void OpDir(const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)        { assert(0); }
  void HermOp(const Field&, Field&)                       { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&){ assert(0); }
};

// ---- exact dense diagonalisation of D_W (small lattices; validation only) ----
template<class Field>
static void denseDiag(LinearOperatorBase<Field>& Op, GridCartesian* grid,
                      const std::string& fname, const std::string& tag) {
  int V = 1; for (int d = 0; d < Nd; d++) V *= grid->FullDimensions()[d];
  int N = Ns * Nc * V;
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
  std::vector<std::complex<double>> v(es.eigenvalues().data(), es.eigenvalues().data() + N);
  std::sort(v.begin(), v.end(), [](auto a, auto b){ return std::abs(a.imag()) < std::abs(b.imag()); });
  std::ofstream f(fname); f << std::setprecision(12);
  for (auto& z : v) f << tag << "  " << z.real() << "  " << z.imag() << "\n";
  std::cout << GridLogMessage << "denseDiag -> " << fname << std::endl;
}

// ---- CLI helpers ----
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
  std::ofstream f(fn); f << std::setprecision(12);
  for (auto& e : v) f << tag << "  " << e.first.real() << "  " << e.first.imag() << "\n";
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  RealD mass   = std::stod(getOpt(argc, argv, "--mass",   "0"));
  int   steps  = std::stoi(getOpt(argc, argv, "--steps",  "30"));
  int   kdim   = std::stoi(getOpt(argc, argv, "--krylov", std::to_string(2*steps)));
  int   wanted = std::stoi(getOpt(argc, argv, "--wanted", "10"));
  int   cycles = std::stoi(getOpt(argc, argv, "--cycles", "8"));
  RealD tol    = std::stod(getOpt(argc, argv, "--tol",    "1e-10"));  // double precision
  RealD stol   = std::stod(getOpt(argc, argv, "--stol",   "1e-11"));  // inner CG < accept
  int   siter  = std::stoi(getOpt(argc, argv, "--siter",  "30000"));
  RealD accept = std::stod(getOpt(argc, argv, "--accept", "1e-8"));   // converged if raw res < accept
  RealD dedupe = std::stod(getOpt(argc, argv, "--dedupe", "1e-6"));
  double weak  = hasOpt(argc, argv, "--weak") ? std::stod(getOpt(argc, argv, "--weak", "0.1")) : 0.0;
  std::string cfg = getOpt(argc, argv, "--config", "");
  std::string out = getOpt(argc, argv, "--out",    "g5bl_evals.dat");
  std::string tag = getOpt(argc, argv, "--tag",    std::to_string(mass));
  bool reorth  = !hasOpt(argc, argv, "--noreorth");
  bool iso     = hasOpt(argc, argv, "--isolation-only");
  bool cold    = hasOpt(argc, argv, "--cold");
  bool useG5bl = hasOpt(argc, argv, "--g5bl");

  std::vector<double> sigmas;
  if      (hasOpt(argc, argv, "--shift-sweep")) sigmas = parseSweep(getOpt(argc, argv, "--shift-sweep", ""));
  else if (hasOpt(argc, argv, "--shift"))       sigmas = { std::stod(getOpt(argc, argv, "--shift", "0")) };
  bool shiftMode = !sigmas.empty();

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::cout << GridLogMessage << "lattice=" << GridDefaultLatt() << " mass=" << mass
            << (shiftMode ? "  shift-invert" : "  direct")
            << "  solver=" << (useG5bl ? "Gamma5BlockLanczos" : "RefinedArnoldi") << std::endl;

  // ---- gauge field ----
  LatticeGaugeField Umu(UGrid);
  if (cold) SU<Nc>::ColdConfiguration(Umu);
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

  if (hasOpt(argc, argv, "--dense")) { denseDiag(DW, UGrid, out + ".dense", tag); Grid_finalize(); return 0; }

  auto g5sort = shiftMode ? G5SortAbsDescending : G5SortAbsImagAscending;  // nearest-sigma / near-real
  auto rasort = shiftMode ? RASortAbsDescending : RASortAbsImagAscending;
  auto lamOf  = [&](std::complex<double> mu, double sg){ return shiftMode ? sg + 1.0/mu : mu; };

  // raw Euclidean residual of (u, lambda) against D_W
  FermionField wbuf(UGrid);
  auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
    DW.Op(u, wbuf); ComplexD lf(lam.real(), lam.imag());
    FermionField t(UGrid); t = wbuf - u * lf; return std::sqrt(norm2(t)/norm2(u)); };

  // ============================ --compare (equal Krylov dim) ============================
  if (hasOpt(argc, argv, "--compare")) {
    double sg = shiftMode ? sigmas[0] : 0.0;
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, shiftMode ? mass - sg : mass, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    LinearOperatorBase<FermionField>& M = shiftMode
        ? (LinearOperatorBase<FermionField>&)SI : (LinearOperatorBase<FermionField>&)DW;
    auto distinct = [&](const std::vector<std::complex<double>>& L, const std::vector<FermionField>& U){
      std::vector<EvalRes> C;
      for (int i = 0; i < (int)L.size(); i++) { double r = rawRes(U[i], L[i]); if (r < accept) addDedup(C, L[i], r, dedupe); }
      return (int)C.size(); };

    SI.nApply = SI.nCG = 0; GridStopWatch ta; ta.Start();
    RefinedArnoldi<FermionField> a(M, UGrid, tol, 0); a.setRawCheck(&DW, sg, shiftMode);
    a(v0, kdim, rasort); ta.Stop();
    std::vector<std::complex<double>> La; std::vector<FermionField> Ua;
    for (int i = 0; i < (int)a.getEvals().size(); i++) { La.push_back(lamOf(a.getEvals()(i), sg)); Ua.push_back(a.getEvecs()[i]); }
    long a_cg = SI.nCG; int a_c = distinct(La, Ua);

    SI.nApply = SI.nCG = 0; GridStopWatch tg; tg.Start();
    Gamma5BlockLanczos<FermionField> g(M, UGrid, gamma5, tol, 0); g.setRawCheck(&DW, sg, shiftMode);
    g(v0, v1, steps, reorth, g5sort);   // single pass: equal Krylov dim (2*steps) to RefinedArnoldi
    tg.Stop();
    std::vector<std::complex<double>> Lg; std::vector<FermionField> Ug;
    for (int i = 0; i < (int)g.getEvals().size(); i++) { Lg.push_back(lamOf(g.getEvals()(i), sg)); Ug.push_back(g.getEvecs()[i]); }
    long g_cg = SI.nCG; int g_c = distinct(Lg, Ug);

    std::cout << GridLogMessage << "HEAD-TO-HEAD sigma=" << sg << " accept=" << accept
              << "  (distinct converged modes)" << std::endl;
    std::cout << GridLogMessage << std::setw(20) << "method" << std::setw(12) << "CG_iters"
              << std::setw(10) << "time(s)" << std::setw(8) << "conv" << std::endl;
    std::cout << GridLogMessage << std::setw(20) << "RefinedArnoldi" << std::setw(12) << a_cg
              << std::setw(10) << ta.useconds()*1e-6 << std::setw(8) << a_c << std::endl;
    std::cout << GridLogMessage << std::setw(20) << "Gamma5BlockLanczos" << std::setw(12) << g_cg
              << std::setw(10) << tg.useconds()*1e-6 << std::setw(8) << g_c << std::endl;
    Grid_finalize(); return 0;
  }

  // ===================== --history (residual & #converged vs Krylov dim) =====================
  if (hasOpt(argc, argv, "--history")) {
    double sg = shiftMode ? sigmas[0] : 0.0;
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, shiftMode ? mass - sg : mass, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
    LinearOperatorBase<FermionField>& M = shiftMode
        ? (LinearOperatorBase<FermionField>&)SI : (LinearOperatorBase<FermionField>&)DW;
    auto stats = [&](const std::vector<RealD>& R){               // min, #<1e-6, #<1e-10
      double mn = 1e30; int b6 = 0, b10 = 0;
      for (double r : R) { mn = std::min(mn, r); if (r<1e-6) b6++; if (r<1e-10) b10++; }
      return std::array<double,3>{mn, (double)b6, (double)b10}; };

    RefinedArnoldi<FermionField> a(M, UGrid, tol, 0); a.setRawCheck(&DW, sg, shiftMode);
    a.buildArnoldi(v0, kdim);
    std::ofstream fa(out + ".refinedarnoldi.hist");
    fa << "# krylov_dim min_raw_res n_below_1e-6 n_below_1e-10\n";
    for (int m = 1; m <= a.getNumSteps(); m++) {
      a.extractRitzAt(m, rasort); auto s = stats(a.getResiduals());
      fa << m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }

    Gamma5BlockLanczos<FermionField> g(M, UGrid, gamma5, tol, 0); g.setRawCheck(&DW, sg, shiftMode);
    g(v0, v1, steps, reorth, g5sort);
    std::ofstream fg(out + ".g5bl.hist");
    fg << "# krylov_dim min_raw_res n_below_1e-6 n_below_1e-10\n";
    for (int m = 1; m <= g.getNumSteps(); m++) {
      g.extractRitzAt(m, g5sort); auto s = stats(g.getResiduals());
      fg << 2*m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";   // block size 2
    }
    std::cout << GridLogMessage << "histories -> " << out << ".{refinedarnoldi,g5bl}.hist" << std::endl;
    Grid_finalize(); return 0;
  }

  // ===================== eigenvalue computation (direct / shift-invert / sweep) =====================
  // Run the chosen solver on operator M (D_W direct, or (D_W-sigma)^-1) and keep the
  // raw-residual-converged modes, de-duplicated across shift windows.
  auto solveOn = [&](LinearOperatorBase<FermionField>& M, double sg,
                     std::vector<EvalRes>& out_modes)->int {
    std::vector<std::complex<double>> L; std::vector<FermionField> U; std::vector<RealD> R;
    if (useG5bl) {
      Gamma5BlockLanczos<FermionField> g(M, UGrid, gamma5, tol, 1); g.setRawCheck(&DW, sg, shiftMode);
      if (iso) g(v0, v1, steps, reorth, g5sort); else g.thickRestart(v0, v1, cycles, steps, wanted, reorth, g5sort);
      for (int i = 0; i < (int)g.getEvals().size(); i++) { L.push_back(lamOf(g.getEvals()(i), sg)); R.push_back(g.getResiduals()[i]); U.push_back(g.getEvecs()[i]); }
    } else {
      RefinedArnoldi<FermionField> a(M, UGrid, tol, 1); a.setRawCheck(&DW, sg, shiftMode);
      a(v0, kdim, rasort);
      for (int i = 0; i < (int)a.getEvals().size(); i++) { L.push_back(lamOf(a.getEvals()(i), sg)); R.push_back(a.getResiduals()[i]); U.push_back(a.getEvecs()[i]); }
    }
    int nc = 0;
    for (int i = 0; i < (int)L.size(); i++) if (R[i] < accept) { addDedup(out_modes, L[i], R[i], dedupe); nc++; }
    return nc;
  };

  std::vector<EvalRes> collected;
  GridStopWatch sw; sw.Start();
  if (!shiftMode) {
    int nc = solveOn(DW, 0.0, collected);
    std::cout << GridLogMessage << "direct: " << nc << " converged" << std::endl;
  } else {
    for (double sg : sigmas) {
      WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
      ConjugateGradient<FermionField> cg(stol, siter, false);
      ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);
      int nc = solveOn(SI, sg, collected);
      std::cout << GridLogMessage << "sigma=" << sg << ": " << nc << " converged; total " << collected.size()
                << "  (" << SI.nApply << " apps, " << SI.nCG << " CG its)" << std::endl;
    }
  }
  sw.Stop();
  std::cout << GridLogMessage << "collected " << collected.size() << " distinct eigenvalues in "
            << sw.useconds()*1e-6 << " s -> " << out << std::endl;
  writeEvals(out, tag, collected);

  Grid_finalize();
  return 0;
}
