/*
  Lever-1 consistency check: optimised padded Rs staple vs unoptimised reference.

  Compares Rect_Stout::RectStaplePaddedRs (depth-2 PaddedCell + GeneralLocalStencil)
  against WilsonLoops::RectStapleUnoptimisedRs, in the BaseSmear flow_kernel=2
  convention Cmu = adj(rho * Stap), for all four directions.

  Criterion: bitwise match expected (per-site products are evaluated in the
  adjoint-reversed order, which is exact in FP); tolerance fallback 1e-24 on
  |diff|^2 / |ref|^2 in case compiler FMA contraction differs between paths.
*/
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  typedef PeriodicGimplD Gimpl;
  typedef typename Gimpl::GaugeField GaugeField;
  typedef typename Gimpl::GaugeLinkField GaugeLinkField;

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
							GridDefaultSimd(Nd, vComplexD::Nsimd()),
							GridDefaultMpi());
  std::vector<int> seeds({1, 2, 3, 5});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds);

  GaugeField U(UGrid);
  SU<Nc>::HotConfiguration(RNG4, U); // never unit gauge: it hides path-convention errors

  RealD rho = 0.124;

  PaddedCell Ghost(2, UGrid);
  GaugeField gU = Ghost.ExchangePeriodic(U);
  GridBase *ggrid = Ghost.grids[Nd - 1];

  WilsonLoops<Gimpl> WL;
  int fail = 0;
  for (int mu = 0; mu < Nd; mu++) {
    GeneralLocalStencil gStencil = Rect_Stout<Gimpl>::RectStapleStencilRs(ggrid, mu);

    GaugeLinkField gC(ggrid);
    Rect_Stout<Gimpl>::RectStaplePaddedRs(gC, gU, gStencil, mu, rho);
    GaugeLinkField Copt(UGrid);
    Copt = Ghost.Extract(gC);

    GaugeLinkField Cref(UGrid);
    WL.RectStapleUnoptimisedRs(Cref, U, mu);
    Cref = adj(rho * Cref);

    GaugeLinkField diff(UGrid);
    diff = Copt - Cref;
    RealD nd = norm2(diff);
    RealD nr = norm2(Cref);
    bool ok = (nd <= 1e-24 * nr);
    std::cout << GridLogMessage << "mu=" << mu
	      << "  |Cref|^2 = " << nr
	      << "  |Copt-Cref|^2 = " << nd
	      << (nd == 0.0 ? "  BITWISE MATCH" : (ok ? "  tolerance match" : "  MISMATCH"))
	      << std::endl;
    if (!ok) fail++;
  }

  if (fail) std::cout << GridLogMessage << "FAIL: " << fail << " direction(s) mismatch" << std::endl;
  else      std::cout << GridLogMessage << "PASS" << std::endl;

  Grid_finalize();
  return fail ? 1 : 0;
}
