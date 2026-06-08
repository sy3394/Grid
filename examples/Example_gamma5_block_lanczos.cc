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
      --refined        refined (Euclidean) Ritz extraction -- fixes the oblique
                       gamma5-Galerkin penalty (12x more converged modes on
                       interacting fields at equal cost); still g5bl
      --check          verify/accept each eigenpair by the RAW Euclidean residual
      --weak eps       slight interacting perturbation of the free field
      --dense          exact dense diagonalisation reference (small lattices)
      --compare        g5bl vs plain Arnoldi head-to-head (same operator)
      --diag           bottleneck diagnostics (kappa(Gamma), eta, ghosts, cycles)
      --history        Euclidean residual vs Krylov-dim history (+ --refined)
      --degen eps      override the serious-breakdown / look-ahead threshold

*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/Gamma5BlockLanczos.h>

using namespace std;
using namespace Grid;

// Double precision (run on a cluster with adequate memory; no precision-change
// dance, and the breakdown thresholds sit well above the double noise floor).
typedef WilsonFermionD                         WilsonOp;
typedef typename WilsonFermionD::FermionField  FermionField;

// ---- shift-invert operator: Op(in) = (D_W - sigma)^{-1} in via normal-eq CG ----
template<class Matrix, class Field>
class ShiftInvertNE : public LinearOperatorBase<Field> {
  Matrix&                           D_;       // D_W(mass - sigma)
  MdagMLinearOperator<Matrix,Field> MdagM_;
  ConjugateGradient<Field>&         cg_;
public:
  // cost counters: each Op() is one "inverter solve"; nCG accumulates CG iterations
  long nOp = 0;
  long nCG = 0;
  ShiftInvertNE(Matrix& D, ConjugateGradient<Field>& cg) : D_(D), MdagM_(D), cg_(cg) {}
  void Op(const Field& in, Field& out) {
    Field Mdagb(in.Grid());
    D_.Mdag(in, Mdagb);          // D^dag b
    out = Zero();
    cg_(MdagM_, Mdagb, out);     // (D^dag D)^{-1} D^dag b = D^{-1} b
    nOp++;
    nCG += cg_.IterationsToComplete;
  }
  void AdjOp(const Field& in, Field& out)                            { assert(0); }
  void OpDiag(const Field& in, Field& out)                           { assert(0); }
  void OpDir(const Field& in, Field& out, int dir, int disp)         { assert(0); }
  void OpDirAll(const Field& in, std::vector<Field>& out)            { assert(0); }
  void HermOp(const Field& in, Field& out)                           { assert(0); }
  void HermOpAndNorm(const Field& in, Field& out, RealD& a, RealD& b){ assert(0); }
};

// ---- exact dense diagonalisation of D_W (small lattices only) ----
// Builds the full N x N matrix (N = 12*Volume) by applying the operator to each
// unit basis vector, then diagonalises with Eigen.  Reference for validation on
// an interacting background where no analytic spectrum exists.
template<class Field>
static void denseDiag(LinearOperatorBase<Field>& Op, GridCartesian* grid,
                      const std::string& fname, const std::string& tag) {
  int V = 1; for (int d = 0; d < Nd; d++) V *= grid->FullDimensions()[d];
  int N = Ns * Nc * V;
  std::cout << GridLogMessage << "denseDiag: building " << N << " x " << N
            << " matrix (" << N << " matvecs)..." << std::endl;
  Eigen::MatrixXcd M(N, N);
  Field ek(grid), Dek(grid);
  typedef typename Field::scalar_object sobj;
  std::vector<sobj> col;                    // V entries, lexicographic order
  for (int k = 0; k < N; k++) {
    int site = k / (Ns*Nc), sc = k % (Ns*Nc), s = sc / Nc, c = sc % Nc;
    Coordinate coor(Nd); Lexicographic::CoorFromIndex(coor, site, grid->FullDimensions());
    ek = Zero();
    sobj o; o = Zero(); o()(s)(c) = Complex(1.0, 0.0);
    pokeSite(o, ek, coor);
    Op.Op(ek, Dek);                         // column k of D_W
    unvectorizeToLexOrdArray(col, Dek);     // bulk read of the whole column
    for (int j = 0; j < V; j++)
      for (int sj = 0; sj < Ns; sj++)
        for (int cc = 0; cc < Nc; cc++)
          M(j*(Ns*Nc) + sj*Nc + cc, k) = static_cast<std::complex<double>>(col[j]()(sj)(cc));
    if ((k % 512) == 0)
      std::cout << GridLogMessage << "denseDiag: column " << k << "/" << N << std::endl;
  }
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(M, false);
  Eigen::VectorXcd lam = es.eigenvalues();
  std::vector<std::complex<double>> v(lam.data(), lam.data() + lam.size());
  std::sort(v.begin(), v.end(), [](auto a, auto b){
    return std::abs(a.imag()) < std::abs(b.imag()); });
  std::ofstream f(fname); f << std::setprecision(10);
  for (auto& z : v) f << tag << "  " << z.real() << "  " << z.imag() << "\n";
  f.close();
  std::cout << GridLogMessage << "denseDiag: " << N << " exact eigenvalues written to "
            << fname << std::endl;
}

// ---- plain Euclidean Arnoldi on a given operator (apples-to-apples baseline) ----
// Same shift-invert operator Op = (D_W - sigma)^{-1} as g5bl, but the OUTER
// eigensolver is standard Arnoldi (full Euclidean orthogonalisation, upper
// Hessenberg).  Returns (lambda = sigma + 1/theta, residual) per Ritz value and
// the Arnoldi residual estimate.  matvecs = m (one Op per step).
template<class Field>
static void arnoldiShiftInvert(LinearOperatorBase<Field>& Op, GridBase* grid,
                               const Field& v0, int m, double sigma,
                               std::vector<std::complex<double>>& lambdas,
                               std::vector<Field>& evecs) {
  std::vector<Field> V;
  Field v(grid); v = v0;
  v = v * (1.0 / std::sqrt(norm2(v)));
  V.push_back(v);
  Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(m + 1, m);
  int mdone = m;
  for (int j = 0; j < m; j++) {
    Field w(grid); Op.Op(V[j], w);
    for (int i = 0; i <= j; i++) {
      auto h = innerProduct(V[i], w);
      H(i, j) = std::complex<double>((double)real(h), (double)imag(h));
      w = w - V[i] * h;
    }
    double hn = std::sqrt(norm2(w));
    if (j + 1 <= m) H(j + 1, j) = hn;
    if (hn < 1e-12) { mdone = j + 1; break; }
    if (j + 1 < m) V.push_back(w * (1.0 / hn));
  }
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, mdone, mdone));
  auto lam = es.eigenvalues(); auto Y = es.eigenvectors();
  for (int j = 0; j < mdone; j++) {
    lambdas.push_back(sigma + 1.0 / lam(j));
    Field uj(grid); uj = Zero();
    for (int k = 0; k < mdone && k < (int)V.size(); k++) uj = uj + V[k] * Y(k, j);
    evecs.push_back(uj);
  }
}

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
  bool check   = hasOpt(argc, argv, "--check");   // verify each eigenpair against raw D_W
  bool refined = hasOpt(argc, argv, "--refined"); // refined (Euclidean) Ritz extraction
  RealD degen  = std::stod(getOpt(argc, argv, "--degen", "-1")); // override breakdown threshold

  // shift list: --shift-sweep lo:hi:n  (sweep)  |  --shift sigma  (single)  |  none (direct)
  std::vector<double> sigmas;
  bool sweepMode  = hasOpt(argc, argv, "--shift-sweep");
  bool singleShift= hasOpt(argc, argv, "--shift");
  if (sweepMode)        sigmas = parseSweep(getOpt(argc, argv, "--shift-sweep", ""));
  else if (singleShift) sigmas = { std::stod(getOpt(argc, argv, "--shift", "0")) };
  bool shiftMode = !sigmas.empty();

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
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
  double weak = hasOpt(argc, argv, "--weak") ? std::stod(getOpt(argc, argv, "--weak", "0.1")) : 0.0;
  LatticeGaugeField Umu(UGrid);
  if (cold) {
    std::cout << GridLogMessage << "COLD (unit) gauge: free Wilson operator." << std::endl;
    SU<Nc>::ColdConfiguration(Umu);
  } else if (weak > 0.0) {
    std::cout << GridLogMessage << "WEAK gauge: U_mu = exp(i * " << weak
              << " * random algebra) -- slight perturbation of the free field."
              << std::endl;
    GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4});
    LatticeColourMatrix Ulink(UGrid);
    for (int mu = 0; mu < Nd; mu++) {
      SU<Nc>::LieRandomize(pRNG, Ulink, weak);
      PokeIndex<LorentzIndex>(Umu, Ulink, mu);
    }
  } else if (cfg.empty()) {
    std::cout << GridLogMessage << "HOT random gauge (smoke test only)." << std::endl;
    GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4});
    SU<Nc>::HotConfiguration(pRNG, Umu);
  } else {
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, cfg);
    std::cout << GridLogMessage << "Loaded NERSC config: " << cfg
              << "  plaquette=" << header.plaquette << std::endl;
  }

  Gamma G5(Gamma::Algebra::Gamma5);
  auto gamma5 = [&G5](const FermionField& in, FermionField& out){ out = G5 * in; };

  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonOp::ImplParams wpar(boundary);

  GridParallelRNG RNG(UGrid); RNG.SeedFixedIntegers({5,6,7,8});
  FermionField v0(UGrid), v1(UGrid);
  random(RNG, v0); random(RNG, v1);

  // exact dense reference (small lattices): writes ALL D_W eigenvalues
  if (hasOpt(argc, argv, "--dense")) {
    WilsonOp Dwd(Umu, *UGrid, *UrbGrid, mass, wpar);
    NonHermitianLinearOperator<WilsonOp, FermionField> Ld(Dwd);
    denseDiag(Ld, UGrid, out + ".dense", tag);
  }

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
      if (degen>0) g.setDegenRel(degen);
    if (refined) g.setRefined(true);
      g(v0, v1, steps, reorth, G5SortAbsImagAscending); collect(g, out + ".iso"); }
    if (!isoOnly) {
      std::cout << GridLogMessage << "\n--- THICK RESTART ---" << std::endl;
      Gamma5BlockLanczos<FermionField> g(DLinOp, UGrid, gamma5, tol, 1);
      if (degen>0) g.setDegenRel(degen);
    if (refined) g.setRefined(true);
      GridStopWatch sw; sw.Start();
      g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsImagAscending);
      sw.Stop();
      collect(g, out);
      std::cout << GridLogMessage << "DIRECT thick-restart solve time: "
                << sw.useconds()*1e-6 << " s" << std::endl;
    }
    std::cout << GridLogMessage << "Done." << std::endl;
    Grid_finalize(); return 0;
  }

  // ================= SHIFT-INVERT (single or sweep) =================
  // Direct operator D_W(mass), for the optional raw-operator eigenpair check.
  WilsonOp Dw_direct(Umu, *UGrid, *UrbGrid, mass, wpar);
  NonHermitianLinearOperator<WilsonOp, FermionField> DLinDirect(Dw_direct);

  // --- common-inverter baseline: one standard Wilson solve D_W(mass) x = b ---
  // (normal-equations CG on D^dag D, the usual robust Wilson inverter) so the
  // eigensolve cost can be quoted in units of a "common inverter" solve.
  long baseCG = 0;
  {
    MdagMLinearOperator<WilsonOp, FermionField> MdagM(Dw_direct);
    ConjugateGradient<FermionField> cgb(stol, siter, false);
    FermionField b(UGrid), Mdb(UGrid), x(UGrid);
    b = v0; Dw_direct.Mdag(b, Mdb); x = Zero();
    GridStopWatch sw; sw.Start(); cgb(MdagM, Mdb, x); sw.Stop();
    baseCG = cgb.IterationsToComplete;
    std::cout << GridLogMessage << "[baseline] one common inverter solve D_W(m="
              << mass << ") x=b : " << baseCG << " CG iters, "
              << sw.useconds()*1e-6 << " s" << std::endl;
  }

  // ===== apples-to-apples head-to-head: g5bl vs plain Arnoldi, same operator =====
  // Both invert the SAME (D_W - sigma)^{-1} (same inner CG) and target the same
  // complex eigenvalues; only the outer eigensolver differs.  Equal matvec budget.
  // ===== residual-vs-Krylov-dimension history (Euclidean), g5bl vs Arnoldi =====
  // Writes <out>.g5bl.hist and <out>.arnoldi.hist with columns:
  //   matvecs  krylov_dim  min_raw_res  n_below_1e-2  n_below_1e-4
  // where min_raw_res = smallest ||D_W u - lambda u||/||u|| over all Ritz pairs.
  // matvecs is the cost axis (= inner CG solves); at equal matvecs both methods
  // span the same Krylov dimension, so it is a fair head-to-head.
  if (hasOpt(argc, argv, "--history")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cgh(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cgh);
    FermionField w(UGrid);
    auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
      DLinDirect.Op(u, w); ComplexD lf(lam.real(), lam.imag());
      FermionField t(UGrid); t = w - u * lf; return std::sqrt(norm2(t)/norm2(u)); };
    auto stats = [&](const std::vector<std::complex<double>>& lams,
                     const std::vector<FermionField>& vecs)->std::array<double,3> {
      double mn = 1e30; int b2 = 0, b4 = 0;
      for (int i = 0; i < (int)lams.size(); i++) {
        double r = rawRes(vecs[i], lams[i]);
        mn = std::min(mn, r); if (r < 1e-2) b2++; if (r < 1e-4) b4++;
      }
      return {mn, (double)b2, (double)b4}; };

    // --- g5bl history ---
    Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 0);
    if (degen > 0) g.setDegenRel(degen);
    if (refined) g.setRefined(true);
    g(v0, v1, steps, reorth, G5SortAbsDescending);
    std::ofstream fg(out + ".g5bl.hist");
    fg << "# matvecs krylov_dim min_raw_res n_below_1e-2 n_below_1e-4\n";
    for (int m = 1; m <= g.getNumSteps(); m++) {
      g.extractRitzAt(m, G5SortAbsDescending);
      const auto& ev = g.getEvals(); const auto& uv = g.getEvecs();
      std::vector<std::complex<double>> L; std::vector<FermionField> Vv;
      for (int i = 0; i < (int)ev.size(); i++) { L.push_back(sg + 1.0/ev(i)); Vv.push_back(uv[i]); }
      auto s = stats(L, Vv);
      fg << 2*m << " " << (int)ev.size() << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    fg.close();
    std::cout << GridLogMessage << "g5bl history -> " << out << ".g5bl.hist" << std::endl;

    // --- g5bl with REFINED (Euclidean) extraction ---
    if (refined) {
      Gamma5BlockLanczos<FermionField> gr(SI, UGrid, gamma5, tol, 0);
      if (degen > 0) gr.setDegenRel(degen);
      gr.setRefined(true);
      gr(v0, v1, steps, reorth, G5SortAbsDescending);
      std::ofstream fr(out + ".g5bl_refined.hist");
      fr << "# matvecs krylov_dim min_raw_res n_below_1e-2 n_below_1e-4\n";
      for (int mm = 1; mm <= gr.getNumSteps(); mm++) {
        gr.extractRitzAt(mm, G5SortAbsDescending);
        const auto& ev = gr.getEvals(); const auto& uv = gr.getEvecs();
        std::vector<std::complex<double>> L; std::vector<FermionField> Vv;
        for (int i = 0; i < (int)ev.size(); i++) { L.push_back(sg + 1.0/ev(i)); Vv.push_back(uv[i]); }
        auto s = stats(L, Vv);
        fr << 2*mm << " " << (int)ev.size() << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
      }
      fr.close();
      std::cout << GridLogMessage << "g5bl_refined history -> " << out << ".g5bl_refined.hist" << std::endl;
    }

    // --- Arnoldi history (build V,H once to 2*steps, extract at each m) ---
    int budget = 2 * steps;
    std::vector<FermionField> Va;
    { FermionField v(UGrid); v = v0; v = v * (1.0/std::sqrt(norm2(v))); Va.push_back(v); }
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(budget + 1, budget);
    int adone = budget;
    for (int j = 0; j < budget; j++) {
      FermionField wj(UGrid); SI.Op(Va[j], wj);
      for (int i = 0; i <= j; i++) {
        auto h = innerProduct(Va[i], wj);
        H(i, j) = std::complex<double>((double)real(h), (double)imag(h));
        wj = wj - Va[i] * h;
      }
      double hn = std::sqrt(norm2(wj));
      if (j + 1 <= budget) H(j + 1, j) = hn;
      if (hn < 1e-12) { adone = j + 1; break; }
      if (j + 1 < budget) Va.push_back(wj * (1.0/hn));
    }
    std::ofstream fa(out + ".arnoldi.hist");
    fa << "# matvecs krylov_dim min_raw_res n_below_1e-2 n_below_1e-4\n";
    for (int m = 1; m <= adone; m++) {
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.block(0, 0, m, m));
      auto lam = es.eigenvalues(); auto Y = es.eigenvectors();
      std::vector<std::complex<double>> L; std::vector<FermionField> Vv;
      for (int j = 0; j < m; j++) {
        L.push_back(sg + 1.0/lam(j));
        FermionField uj(UGrid); uj = Zero();
        for (int k = 0; k < m && k < (int)Va.size(); k++) uj = uj + Va[k] * Y(k, j);
        Vv.push_back(uj);
      }
      auto s = stats(L, Vv);
      fa << m << " " << m << " " << s[0] << " " << (int)s[1] << " " << (int)s[2] << "\n";
    }
    fa.close();
    std::cout << GridLogMessage << "Arnoldi history -> " << out << ".arnoldi.hist" << std::endl;
    Grid_finalize(); return 0;
  }

  // ===== diagnostic: pinpoint the g5bl bottleneck (manuscript Sec. 8) =====
  // Tracks per-step kappa(Gamma_k) (oblique-projector conditioning) and
  // eta=||Q_1^dag g5 Q_k|| (loss of gamma5-orthogonality); per-Ritz compares the
  // gamma5-space residual g5bl uses to declare convergence against the RAW
  // Euclidean residual (false convergence => oblique projection is the culprit);
  // checks conjugate-pair symmetry (ghosts); and tracks the best residual per
  // restart cycle (does it worsen each restart?).
  if (hasOpt(argc, argv, "--diag")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cgd(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cgd);
    FermionField w(UGrid);
    auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
      DLinDirect.Op(u, w); ComplexD lf(lam.real(), lam.imag());
      FermionField t(UGrid); t = w - u * lf; return std::sqrt(norm2(t)/norm2(u)); };

    std::cout << GridLogMessage << "\n===== g5bl DIAGNOSTIC (sigma="<<sg<<") =====" << std::endl;
    Gamma5BlockLanczos<FermionField> g(SI, UGrid, gamma5, tol, 0);
    if (degen > 0) g.setDegenRel(degen);
    if (refined) g.setRefined(true);
    g(v0, v1, steps, reorth, G5SortAbsDescending);

    // (A) per-step trajectories
    auto kap = g.getKappaGamma(); auto eta = g.getEtaLoss();
    std::cout << GridLogMessage << "[A] per-step diagnostics:" << std::endl;
    std::cout << GridLogMessage << "   step   kappa(Gamma_k)   eta=||Q1^dag g5 Qk||" << std::endl;
    for (int k = 0; k < (int)kap.size(); k++)
      std::cout << GridLogMessage << "   " << std::setw(4) << k
                << std::setw(16) << kap[k] << std::setw(22) << (k<(int)eta.size()?eta[k]:0.0) << std::endl;

    // (B) gamma5-residual (what g5bl uses) vs raw Euclidean residual, per Ritz pair
    const auto& ev = g.getEvals(); const auto& rs = g.getResiduals(); const auto& uv = g.getEvecs();
    std::cout << GridLogMessage << "[B] per-Ritz: gamma5-residual vs RAW residual "
              << "(false convergence => oblique projection):" << std::endl;
    std::cout << GridLogMessage << std::setw(34) << "lambda"
              << std::setw(13) << "res_g5" << std::setw(13) << "res_raw"
              << std::setw(10) << "raw/g5" << std::setw(8) << "ghost?" << std::endl;
    int nfalse = 0, nshow = std::min((int)ev.size(), 4*std::max(1,wanted));
    for (int i = 0; i < nshow; i++) {
      std::complex<double> lam = sg + 1.0/ev(i);
      double rg = rs[i], rr = rawRes(uv[i], lam);
      bool ghost = (rg < accept && rr >= accept);     // g5 says converged, raw says no
      if (ghost) nfalse++;
      std::cout << GridLogMessage << "  ("<<std::setw(9)<<lam.real()<<","<<std::setw(11)<<lam.imag()<<")"
                << std::setw(13) << rg << std::setw(13) << rr
                << std::setw(10) << rr/std::max(rg,1e-30) << std::setw(8) << (ghost?"YES":"-") << std::endl;
    }
    std::cout << GridLogMessage << "   false-convergence count (res_g5<acc but res_raw>=acc): "
              << nfalse << " / " << nshow << std::endl;

    // (C) conjugate-pair symmetry of the spectrum (broken => spurious modes)
    int unpaired = 0;
    for (int i = 0; i < (int)ev.size(); i++) {
      std::complex<double> li = sg + 1.0/ev(i);
      double best = 1e30;
      for (int j = 0; j < (int)ev.size(); j++) if (j!=i) {
        std::complex<double> lj = sg + 1.0/ev(j);
        best = std::min(best, std::abs(lj - std::conj(li)));
      }
      if (best > 1e-3) unpaired++;
    }
    std::cout << GridLogMessage << "[C] conjugate-pair symmetry: " << unpaired
              << " / " << ev.size() << " eigenvalues lack a conjugate partner (within 1e-3)" << std::endl;

    // (D) does the residual worsen across restart cycles?
    Gamma5BlockLanczos<FermionField> gr(SI, UGrid, gamma5, tol, 0);
    if (degen > 0) gr.setDegenRel(degen);
    gr.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsDescending);
    auto cb = gr.getCycleBestRes();
    std::cout << GridLogMessage << "[D] best gamma5-residual per restart cycle "
              << "(increasing => restart degrades):" << std::endl;
    for (int c = 0; c < (int)cb.size(); c++)
      std::cout << GridLogMessage << "   cycle " << std::setw(3) << c
                << "   best res_g5 = " << cb[c] << std::endl;
    std::cout << GridLogMessage << "   locked pairs: " << gr.getNumLocked() << std::endl;
    Grid_finalize(); return 0;
  }

  if (hasOpt(argc, argv, "--compare")) {
    double sg = sigmas[0];
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cgc(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cgc);
    FermionField w(UGrid);
    auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
      DLinDirect.Op(u, w);
      ComplexD lf(lam.real(), lam.imag());
      FermionField t(UGrid); t = w - u * lf;
      return std::sqrt(norm2(t) / norm2(u));
    };
    int budget = 2 * steps;   // g5bl block-2 does 2 matvecs/step; match Arnoldi steps

    // --- g5bl (gamma5-block Lanczos) ---
    SI.nOp = SI.nCG = 0;
    Gamma5BlockLanczos<FermionField> gg(SI, UGrid, gamma5, tol, 0);
    if (degen > 0) gg.setDegenRel(degen);
    if (refined) gg.setRefined(true);
    GridStopWatch sw1; sw1.Start();
    gg(v0, v1, steps, reorth, G5SortAbsDescending);
    sw1.Stop();
    long g_nOp = SI.nOp, g_nCG = SI.nCG;
    int g_conv = 0;
    { const auto& ev = gg.getEvals(); const auto& uv = gg.getEvecs();
      for (int i = 0; i < (int)ev.size() && i < (int)uv.size(); i++) {
        std::complex<double> lam = sg + 1.0 / ev(i);   // theta -> D_W eigenvalue
        if (rawRes(uv[i], lam) < accept) g_conv++;
      } }

    // --- plain Arnoldi, equal matvec budget ---
    SI.nOp = SI.nCG = 0;
    std::vector<std::complex<double>> a_lam; std::vector<FermionField> a_vec;
    GridStopWatch sw2; sw2.Start();
    arnoldiShiftInvert(SI, UGrid, v0, budget, sg, a_lam, a_vec);
    sw2.Stop();
    long a_nOp = SI.nOp, a_nCG = SI.nCG;
    int a_conv = 0;
    for (int i = 0; i < (int)a_lam.size(); i++)
      if (rawRes(a_vec[i], a_lam[i]) < accept) a_conv++;

    std::cout << GridLogMessage << "\n===== HEAD-TO-HEAD (same D_W, same shift-invert, sigma="
              << sg << ", accept=" << accept << ") =====" << std::endl;
    std::cout << GridLogMessage << std::setw(22) << "method"
              << std::setw(12) << "matvecs" << std::setw(12) << "CG iters"
              << std::setw(10) << "time(s)" << std::setw(14) << "conv(raw<acc)" << std::endl;
    std::cout << GridLogMessage << std::setw(22) << "g5-block Lanczos"
              << std::setw(12) << g_nOp << std::setw(12) << g_nCG
              << std::setw(10) << sw1.useconds()*1e-6 << std::setw(14) << g_conv << std::endl;
    std::cout << GridLogMessage << std::setw(22) << "plain Arnoldi"
              << std::setw(12) << a_nOp << std::setw(12) << a_nCG
              << std::setw(10) << sw2.useconds()*1e-6 << std::setw(14) << a_conv << std::endl;
    Grid_finalize(); return 0;
  }

  std::vector<EvalRes> collected;
  double tSolve = 0.0;
  long totOp = 0, totCG = 0;
  for (size_t is = 0; is < sigmas.size(); is++) {
    double sg = sigmas[is];
    std::cout << GridLogMessage << "\n--- SHIFT-INVERT sigma=" << sg
              << "  (" << is+1 << "/" << sigmas.size() << ", inverting D_W(m="
              << mass - sg << ")) ---" << std::endl;

    WilsonOp Dshift(Umu, *UGrid, *UrbGrid, mass - sg, wpar);
    ConjugateGradient<FermionField> cg(stol, siter, /*err_on_no_conv=*/false);
    ShiftInvertNE<WilsonOp, FermionField> SIop(Dshift, cg);
    Gamma5BlockLanczos<FermionField> g(SIop, UGrid, gamma5, tol, 1);
    if (degen>0) g.setDegenRel(degen);
    if (refined) g.setRefined(true);

    GridStopWatch sw; sw.Start();
    if (isoOnly) g(v0, v1, steps, reorth, G5SortAbsDescending);
    else         g.thickRestart(v0, v1, cycles, steps, wanted, reorth, G5SortAbsDescending);
    sw.Stop();
    tSolve += sw.useconds()*1e-6;

    const auto& ev = g.getEvals(); const auto& rs = g.getResiduals();
    const auto& uv = g.getEvecs();
    int nconv = 0;
    FermionField w(UGrid);
    for (int i = 0; i < (int)ev.size(); i++) {
      std::complex<double> lam = sg + 1.0/ev(i);     // map back to D_W eigenvalue
      double resAccept = rs[i];                       // default: Lanczos theta-residual
      // raw-operator residual ||D_W u - lambda u||/||u|| -- the honest convergence
      // measure (rejects shift-invert ghosts whose theta-residual looks small).
      if (check && i < (int)uv.size()) {
        DLinDirect.Op(uv[i], w);
        ComplexD lamf(lam.real(), lam.imag());
        w = w - uv[i] * lamf;
        resAccept = std::sqrt(norm2(w) / norm2(uv[i]));   // use RAW residual to accept
      }
      if (resAccept >= accept) continue;
      addDedup(collected, lam, resAccept, dedupe);
      nconv++;
    }
    totOp += SIop.nOp; totCG += SIop.nCG;
    std::cout << GridLogMessage << "   sigma=" << sg << ": " << nconv
              << " converged modes (res<" << accept << ")"
              << "   [outer " << (isoOnly?steps:steps*cycles) << " Lanczos steps, "
              << SIop.nOp << " inverter solves, " << SIop.nCG << " CG iters, "
              << sw.useconds()*1e-6 << " s]   running total " << collected.size() << std::endl;
  }
  // --- performance summary, quoted against the common inverter ---
  std::cout << GridLogMessage << "\n=== performance (shift-invert) ===" << std::endl;
  std::cout << GridLogMessage << "  shifts            : " << sigmas.size() << std::endl;
  std::cout << GridLogMessage << "  inverter solves   : " << totOp
            << "   (each ~ one common inverter solve)" << std::endl;
  std::cout << GridLogMessage << "  total CG iters    : " << totCG
            << "   (baseline 1 inverter = " << baseCG << " iters)" << std::endl;
  std::cout << GridLogMessage << "  ~inverter-equiv   : " << (baseCG>0 ? (double)totCG/baseCG : 0.0)
            << " common-inverter solves" << std::endl;
  std::cout << GridLogMessage << "  collected modes   : " << collected.size() << std::endl;
  std::cout << GridLogMessage << "  total solve time  : " << tSolve << " s"
            << "   (" << (collected.size()? tSolve/collected.size():0.0) << " s/mode)" << std::endl;

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
