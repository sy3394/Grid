/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Wilson_DW_spectrum.cc

    Compute the interior eigenvalues of the Wilson Dirac operator D_W in a window of
    the real axis, Re(lambda) in [lo, hi], for a NERSC gauge configuration.

    Method: a shift-invert SWEEP over real sigma tiling [lo, hi].  For each sigma the
    solver runs on M = (D_W - sigma)^{-1} (assembled from Grid's ConjugateGradient via
    the normal equations) and converges the eigenvalues of D_W nearest sigma; a Ritz
    value mu maps back as lambda = sigma + 1/mu.  Real sigma keeps (D_W - sigma)
    gamma5-Hermitian.  Converged modes (raw residual ||D_W u - lambda u||/||u|| < accept)
    are accumulated and de-duplicated across windows, then written out (those with
    Re(lambda) in [lo, hi]).

    Solver: RefinedArnoldi (Grid/algorithms/iterative/RefinedArnoldi.h).

    Usage (32^4 example):
      Wilson_DW_spectrum --grid 32.32.32.32 --mpi 1.2.2.2 --config <NERSC> \
          --mass 0 --window 0:2 --nshift 21 --krylov 80 \
          --stol 1e-11 --accept 1e-8 --tag <mdtime> --out evals_DW.dat
    Output rows:  <tag>  Re(lambda)  Im(lambda)   (sorted by |Im|).

*************************************************************************************/

#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/RefinedArnoldi.h>

using namespace Grid;

typedef WilsonFermionD                         WilsonOp;
typedef typename WilsonFermionD::FermionField  FermionField;

// (D_W - sigma)^{-1} via normal-equations CG (compose Grid's CG + MdagMLinearOperator).
template<class Matrix, class Field>
class ShiftInvertNE : public LinearOperatorBase<Field> {
  Matrix&                           D_;
  MdagMLinearOperator<Matrix,Field> MdagM_;
  ConjugateGradient<Field>&         cg_;
public:
  long nApply = 0, nCG = 0;
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

static std::string argOr(int argc, char** argv, const std::string& k, const std::string& d) {
  return GridCmdOptionExists(argv, argv + argc, k) ? GridCmdOptionPayload(argv, argv + argc, k) : d;
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  RealD mass   = std::stod(argOr(argc, argv, "--mass",   "0"));
  int   nshift = std::stoi(argOr(argc, argv, "--nshift", "21"));
  int   kdim   = std::stoi(argOr(argc, argv, "--krylov", "80"));
  RealD stol   = std::stod(argOr(argc, argv, "--stol",   "1e-11"));
  int   siter  = std::stoi(argOr(argc, argv, "--siter",  "30000"));
  RealD accept = std::stod(argOr(argc, argv, "--accept", "1e-8"));
  RealD dedupe = std::stod(argOr(argc, argv, "--dedupe", "1e-6"));
  std::string cfg = argOr(argc, argv, "--config", "");
  std::string out = argOr(argc, argv, "--out",    "evals_DW.dat");
  std::string tag = argOr(argc, argv, "--tag",    "0");
  double lo = 0.0, hi = 2.0;
  sscanf(argOr(argc, argv, "--window", "0:2").c_str(), "%lf:%lf", &lo, &hi);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::cout << GridLogMessage << "D_W spectrum on " << GridDefaultLatt() << " mass=" << mass
            << "  window Re in [" << lo << "," << hi << "]  nshift=" << nshift
            << "  Krylov dim=" << kdim << std::endl;

  // gauge field
  LatticeGaugeField Umu(UGrid);
  if (cfg.empty()) {
    std::cout << GridLogMessage << "no --config: cold (free field)" << std::endl;
    SU<Nc>::ColdConfiguration(Umu);
  } else {
    FieldMetaData h; NerscIO::readConfiguration(Umu, h, cfg);
    std::cout << GridLogMessage << "config " << cfg << " plaquette=" << h.plaquette << std::endl;
  }

  std::vector<Complex> bc = {1,1,1,-1};
  WilsonOp::ImplParams wpar(bc);
  WilsonOp Dw(Umu, *UGrid, *UrbGrid, mass, wpar);
  NonHermitianLinearOperator<WilsonOp, FermionField> DW(Dw);

  GridParallelRNG RNG(UGrid); RNG.SeedFixedIntegers({5,6,7,8});
  FermionField v0(UGrid); random(RNG, v0);
  FermionField wbuf(UGrid);
  auto rawRes = [&](const FermionField& u, std::complex<double> lam)->double {
    DW.Op(u, wbuf); ComplexD lf(lam.real(), lam.imag());
    FermionField t(UGrid); t = wbuf - u * lf; return std::sqrt(norm2(t)/norm2(u)); };

  // collected (lambda, residual), de-duplicated
  typedef std::pair<std::complex<double>, double> EvalRes;
  std::vector<EvalRes> modes;
  auto addDedup = [&](std::complex<double> l, double r){
    for (auto& e : modes) if (std::abs(e.first - l) < dedupe) { if (r < e.second) e = {l, r}; return; }
    modes.push_back({l, r}); };

  long totApply = 0, totCG = 0;
  GridStopWatch sw; sw.Start();
  for (int is = 0; is < nshift; is++) {
    double sg = (nshift == 1) ? 0.5*(lo+hi) : lo + (hi - lo) * is / (nshift - 1);
    WilsonOp Dsh(Umu, *UGrid, *UrbGrid, mass - sg, wpar);   // D_W(m) - sigma = D_W(m - sigma)
    ConjugateGradient<FermionField> cg(stol, siter, false);
    ShiftInvertNE<WilsonOp, FermionField> SI(Dsh, cg);

    RefinedArnoldi<FermionField> a(SI, UGrid, accept, 1);
    a.setRawCheck(&DW, sg, /*shiftInvert=*/true);
    a(v0, kdim, RASortAbsDescending);                       // nearest-sigma modes first

    int nc = 0;
    for (int i = 0; i < (int)a.getEvals().size(); i++) {
      std::complex<double> lam = sg + 1.0/a.getEvals()(i);
      if (a.getResiduals()[i] < accept) { addDedup(lam, a.getResiduals()[i]); nc++; }
    }
    totApply += SI.nApply; totCG += SI.nCG;
    std::cout << GridLogMessage << "sigma=" << sg << ": " << nc << " converged; total "
              << modes.size() << "  (" << SI.nApply << " apps, " << SI.nCG << " CG its)" << std::endl;
  }
  sw.Stop();

  // keep only modes inside the requested window, sort by |Im|
  std::vector<EvalRes> win;
  for (auto& e : modes) if (e.first.real() >= lo && e.first.real() <= hi) win.push_back(e);
  std::sort(win.begin(), win.end(), [](const EvalRes& a, const EvalRes& b){
    return std::abs(a.first.imag()) < std::abs(b.first.imag()); });

  std::ofstream f(out); f << std::setprecision(12);
  for (auto& e : win) f << tag << "  " << e.first.real() << "  " << e.first.imag() << "\n";

  std::cout << GridLogMessage << win.size() << " D_W eigenvalues with Re in [" << lo << "," << hi
            << "]  (" << totApply << " applications, " << totCG << " CG iters, "
            << sw.useconds()*1e-6 << " s)  -> " << out << std::endl;

  Grid_finalize();
  return 0;
}
