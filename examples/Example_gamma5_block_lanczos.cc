/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_gamma5_block_lanczos.cc

    gamma5-Block Lanczos driver for the Wilson Dirac operator D_W.

    Computes eigenvalues of D_W DIRECTLY (not H_W = g5 D_W) using the
    gamma5-Block Lanczos algorithm (Grid/algorithms/iterative/Gamma5BlockLanczos.h).

    Modes:
      DIRECT (default)             : Lanczos on D_W; converges peripheral
                                     eigenvalues (doublers, spectral edges).
      SHIFT-INVERT (--shift sigma) : Lanczos on (D_W - sigma)^{-1}; converges the
                                     eigenvalues of D_W CLOSEST to real sigma (the
                                     interior near-real cluster around Re ~ sigma).
      SWEEP (--shift-sweep lo:hi:n): loop n real shifts from lo to hi in ONE job,
                                     accumulate + de-duplicate the converged modes
                                     of each window.  Resolves a whole interior
                                     band of the near-real D_W spectrum without
                                     multiple submissions.

    Shift-invert (sigma real): D_W(m) - sigma = D_W(m - sigma) (additive Wilson
    mass).  D_W(m - sigma) is inverted by normal-equations CG (HPD operator
    D^dag D), robust where BiCGSTAB on the indefinite operator fails.  A Ritz
    value theta maps back as  lambda = sigma + 1/theta.

    Output (3-column "<tag> Re Im"):
      DIRECT / single-shift : <out>.iso (isolation) and <out> (thick restart)
      SWEEP                 : <out>     (all converged, de-duplicated modes)

    Usage:
      --config PATH    NERSC gauge config       (else --cold, or hot random)
      --cold           unit gauge (free field; analytic spectrum)
      --grid X.Y.Z.T   lattice dims
      --mass m         Wilson bare mass          (default -0.5)
      --shift sigma    single shift-invert about real sigma
      --shift-sweep lo:hi:n   sweep n shifts in [lo,hi] (single submission)
      --accept eps     keep modes with theta-residual < eps as converged (default 1e-3)
      --dedupe eps     merge modes within eps in the complex plane (default 1e-4)
      --stol eps       inner CG tolerance        (default 1e-8)
      --siter N        inner CG max iterations    (default 20000)
      --steps N        Lanczos steps per pass    (default 40)
      --wanted N       wanted conjugate pairs    (default 10)
      --cycles N       thick-restart cycles      (default 8)
      --tol eps        Lanczos residual tol      (default 1e-8)
      --tag s          column-1 label in output  (default = mass)
      --out PATH       output .dat path          (default g5bl_evals.dat)
      --noreorth       disable gamma5-reorthogonalisation
      --isolation-only single pass only (direct / single-shift)

*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/Gamma5BlockLanczos.h>

using namespace std;
using namespace Grid;

typedef WilsonFermionF                         WilsonOp;
typedef typename WilsonFermionF::FermionField  FermionField;

// ---- shift-invert operator: Op(in) = (D_W - sigma)^{-1} in via normal-eq CG ----
template<class Matrix, class Field>
class ShiftInvertNE : public LinearOperatorBase<Field> {
  Matrix&                           D_;       // D_W(mass - sigma)
  MdagMLinearOperator<Matrix,Field> MdagM_;
  ConjugateGradient<Field>&         cg_;
public:
  ShiftInvertNE(Matrix& D, ConjugateGradient<Field>& cg) : D_(D), MdagM_(D), cg_(cg) {}
  void Op(const Field& in, Field& out) {
    Field Mdagb(in.Grid());
    D_.Mdag(in, Mdagb);          // D^dag b
    out = Zero();
    cg_(MdagM_, Mdagb, out);     // (D^dag D)^{-1} D^dag b = D^{-1} b
  }
  void AdjOp(const Field& in, Field& out)                            { assert(0); }
  void OpDiag(const Field& in, Field& out)                           { assert(0); }
  void OpDir(const Field& in, Field& out, int dir, int disp)         { assert(0); }
  void OpDirAll(const Field& in, std::vector<Field>& out)            { assert(0); }
  void HermOp(const Field& in, Field& out)                           { assert(0); }
  void HermOpAndNorm(const Field& in, Field& out, RealD& a, RealD& b){ assert(0); }
};

// ---- command-line helpers ----
static std::string getOpt(int argc, char** argv, const std::string& key,
                          const std::string& def="") {
  for (int i = 1; i < argc-1; i++) if (key == argv[i]) return argv[i+1];
  return def;
}
static bool hasOpt(int argc, char** argv, const std::string& key) {
  for (int i = 1; i < argc; i++) if (key == argv[i]) return true;
  return false;
}

// parse "lo:hi:n" -> n evenly spaced shifts in [lo,hi]
static std::vector<double> parseSweep(const std::string& s) {
  std::vector<double> v;
  double lo, hi; int n;
  if (sscanf(s.c_str(), "%lf:%lf:%d", &lo, &hi, &n) == 3 && n >= 1) {
    if (n == 1) { v.push_back(lo); return v; }
    for (int i = 0; i < n; i++) v.push_back(lo + (hi - lo) * i / (n - 1));
  }
  return v;
}

typedef std::pair<std::complex<double>, double> EvalRes;  // (lambda, theta-residual)

// add lambda to the collection unless a near-duplicate already exists (keep smaller res)
static void addDedup(std::vector<EvalRes>& coll, std::complex<double> lam,
                     double res, double dedupe) {
  for (auto& e : coll)
    if (std::abs(e.first - lam) < dedupe) { if (res < e.second) e = {lam, res}; return; }
  coll.push_back({lam, res});
}

static void writeEvals(const std::string& fname, const std::string& tag,
                       std::vector<EvalRes> v) {
  std::sort(v.begin(), v.end(), [](const EvalRes& a, const EvalRes& b){
    return std::abs(a.first.imag()) < std::abs(b.first.imag());
  });
  std::ofstream f(fname);
  f << std::setprecision(10);
  for (auto& e : v) f << tag << "  " << e.first.real() << "  " << e.first.imag() << "\n";
  f.close();
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  RealD mass   = std::stod(getOpt(argc, argv, "--mass",   "-0.5"));
  int   steps  = std::stoi(getOpt(argc, argv, "--steps",  "40"));
  int   wanted = std::stoi(getOpt(argc, argv, "--wanted", "10"));
  int   cycles = std::stoi(getOpt(argc, argv, "--cycles", "8"));
  RealD tol    = std::stod(getOpt(argc, argv, "--tol",    "1e-8"));
  RealD stol   = std::stod(getOpt(argc, argv, "--stol",   "1e-8"));
  int   siter  = std::stoi(getOpt(argc, argv, "--siter",  "20000"));
  RealD accept = std::stod(getOpt(argc, argv, "--accept", "1e-3"));
  RealD dedupe = std::stod(getOpt(argc, argv, "--dedupe", "1e-4"));
  std::string cfg = getOpt(argc, argv, "--config", "");
  std::string out = getOpt(argc, argv, "--out",    "g5bl_evals.dat");
  std::string tag = getOpt(argc, argv, "--tag",    std::to_string(mass));
  bool reorth  = !hasOpt(argc, argv, "--noreorth");
  bool isoOnly = hasOpt(argc, argv, "--isolation-only");
  bool cold    = hasOpt(argc, argv, "--cold");

  // shift list: --shift-sweep lo:hi:n  (sweep)  |  --shift sigma  (single)  |  none (direct)
  std::vector<double> sigmas;
  bool sweepMode  = hasOpt(argc, argv, "--shift-sweep");
  bool singleShift= hasOpt(argc, argv, "--shift");
  if (sweepMode)        sigmas = parseSweep(getOpt(argc, argv, "--shift-sweep", ""));
  else if (singleShift) sigmas = { std::stod(getOpt(argc, argv, "--shift", "0")) };
  bool shiftMode = !sigmas.empty();

  GridCartesian* UGridD = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::cout << GridLogMessage << "=====================================================" << std::endl;
  std::cout << GridLogMessage << " gamma5-Block Lanczos for D_W   lattice=" << GridDefaultLatt()
            << "  mass=" << mass << std::endl;
  if (sweepMode) {
    std::cout << GridLogMessage << "   MODE = SHIFT-INVERT SWEEP over " << sigmas.size()
              << " shifts: [";
    for (size_t i = 0; i < sigmas.size(); i++) std::cout << sigmas[i] << (i+1<sigmas.size()?", ":"");
    std::cout << "]" << std::endl;
  } else if (singleShift) {
    std::cout << GridLogMessage << "   MODE = SHIFT-INVERT about sigma=" << sigmas[0] << std::endl;
  } else {
    std::cout << GridLogMessage << "   MODE = DIRECT" << std::endl;
  }
  std::cout << GridLogMessage << "   steps=" << steps << " wanted=" << wanted
            << " cycles=" << cycles << " tol=" << tol
            << " accept=" << accept << " stol=" << stol << std::endl;
  std::cout << GridLogMessage << "=====================================================" << std::endl;

  // ---- gauge ----
  LatticeGaugeField  UmuD(UGridD);
  LatticeGaugeFieldF Umu(UGrid);
  if (cold) {
    std::cout << GridLogMessage << "COLD (unit) gauge: free Wilson operator." << std::endl;
    SU<Nc>::ColdConfiguration(UmuD);
  } else if (cfg.empty()) {
    std::cout << GridLogMessage << "HOT random gauge (smoke test only)." << std::endl;
    GridParallelRNG pRNG(UGridD); pRNG.SeedFixedIntegers({1,2,3,4});
    SU<Nc>::HotConfiguration(pRNG, UmuD);
  } else {
    FieldMetaData header;
    NerscIO::readConfiguration(UmuD, header, cfg);
    std::cout << GridLogMessage << "Loaded NERSC config: " << cfg
              << "  plaquette=" << header.plaquette << std::endl;
  }
  precisionChange(Umu, UmuD);

  Gamma G5(Gamma::Algebra::Gamma5);
  auto gamma5 = [&G5](const FermionField& in, FermionField& out){ out = G5 * in; };

  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonOp::ImplParams wpar(boundary);

  GridParallelRNG RNG(UGrid); RNG.SeedFixedIntegers({5,6,7,8});
  FermionField v0(UGrid), v1(UGrid);
  random(RNG, v0); random(RNG, v1);

  // ================= DIRECT mode =================
  if (!shiftMode) {
    WilsonOp Dw(Umu, *UGrid, *UrbGrid, mass, wpar);
    NonHermitianLinearOperator<WilsonOp, FermionField> DLinOp(Dw);

    auto collect = [&](Gamma5BlockLanczos<FermionField>& g, const std::string& fn){
      const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
      std::vector<EvalRes> v(ev.size());
      for (int i = 0; i < (int)ev.size(); i++) v[i] = { ev(i), rs[i] };
      std::sort(v.begin(), v.end(), [](const EvalRes& a, const EvalRes& b){
        return std::abs(a.first.imag()) < std::abs(b.first.imag()); });
      for (int i = 0; i < std::min((int)v.size(), 2*wanted); i++)
        std::cout << GridLogMessage << "   ["<<std::setw(3)<<i<<"]  lambda=("
                  << v[i].first.real()<<", "<<v[i].first.imag()<<")  |res|="<<v[i].second<<std::endl;
      writeEvals(fn, tag, v);
    };

    std::cout << GridLogMessage << "\n--- ISOLATION TEST (single pass) ---" << std::endl;
    { Gamma5BlockLanczos<FermionField> g(DLinOp, UGrid, gamma5, tol, 1);
      g(v0, v1, steps, reorth, G5SortAbsImagAscending); collect(g, out + ".iso"); }
    if (!isoOnly) {
      std::cout << GridLogMessage << "\n--- THICK RESTART ---" << std::endl;
      Gamma5BlockLanczos<FermionField> g(DLinOp, UGrid, gamma5, tol, 1);
      g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsImagAscending);
      collect(g, out);
    }
    std::cout << GridLogMessage << "Done." << std::endl;
    Grid_finalize(); return 0;
  }

  // ================= SHIFT-INVERT (single or sweep) =================
  std::vector<EvalRes> collected;
  for (size_t is = 0; is < sigmas.size(); is++) {
    double sg = sigmas[is];
    std::cout << GridLogMessage << "\n--- SHIFT-INVERT sigma=" << sg
              << "  (" << is+1 << "/" << sigmas.size() << ", inverting D_W(m="
              << mass - sg << ")) ---" << std::endl;

    WilsonOp Dshift(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, /*err_on_no_conv=*/false);
    ShiftInvertNE<WilsonOp, FermionField> SIop(Dshift, cg);
    Gamma5BlockLanczos<FermionField> g(SIop, UGrid, gamma5, tol, 1);

    if (isoOnly) g(v0, v1, steps, reorth, G5SortAbsDescending);
    else         g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsDescending);

    const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
    int nconv = 0;
    for (int i = 0; i < (int)ev.size(); i++) {
      if (rs[i] >= accept) continue;                 // keep only converged windows modes
      std::complex<double> lam = sg + 1.0/ev(i);     // map back to D_W eigenvalue
      addDedup(collected, lam, rs[i], dedupe);
      nconv++;
    }
    std::cout << GridLogMessage << "   sigma=" << sg << ": " << nconv
              << " converged modes (res<" << accept << "); running total "
              << collected.size() << std::endl;
  }

  std::cout << GridLogMessage << "\n=== collected " << collected.size()
            << " distinct converged D_W eigenvalues ===" << std::endl;
  {
    std::vector<EvalRes> v = collected;
    std::sort(v.begin(), v.end(), [](const EvalRes& a, const EvalRes& b){
      return std::abs(a.first.imag()) < std::abs(b.first.imag()); });
    for (int i = 0; i < (int)v.size(); i++)
      std::cout << GridLogMessage << "   ["<<std::setw(3)<<i<<"]  lambda=("
                << v[i].first.real()<<", "<<v[i].first.imag()<<")  |res|="<<v[i].second<<std::endl;
  }
  writeEvals(out, tag, collected);
  std::cout << GridLogMessage << "written to " << out << std::endl;

  std::cout << GridLogMessage << "Done." << std::endl;
  Grid_finalize();
  return 0;
}
