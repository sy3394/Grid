/*
  Driver 2: per-level force AND Jacobian (ln det) consistency, default
  (optimised) vs old implementations (anchor A comparison, see README.md).

  Naming scheme on feature/rect_fthmc-optimise: the optimised routines carry
  the original names and are the default; the old routines take a leading
  `int old` argument and are kept only until consistency is confirmed.

  One fixed hot configuration; one SmearedConfigurationRect with schedule
  mask_types={2,1} (one rect step, one plaquette step; Nsmear=16); both
  overloads run on the same smeared tower, level by level, then the full
  chain rule / level sum end-to-end.

  PASS criteria (last-bit differences expected — fused kernels reorder
  per-site arithmetic; see check_log.md):
    force:  |F_def - F_old|^2 / |F_old|^2  < 1e-24  per level and total
    ln det: |S_def - S_old| / |S_old|      < 1e-12  per level and total
*/
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/GaugeConfigurationRect.h>

using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  typedef PeriodicGimplR Gimpl;

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd, vComplex::Nsimd()),
                                                        GridDefaultMpi());

  std::vector<int> seeds({1, 2, 3, 5});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds);

  LatticeGaugeField U(UGrid);
  SU<Nc>::HotConfiguration(RNG4, U);

  RealD rho = 0.1;
  Smear_Stout<Gimpl> SmearerS(rho);
  Rect_Stout<Gimpl>  SmearerR(0.0, rho, 0.0);
  std::vector<Smear_Stout<Gimpl> *> Stouts = {&SmearerR, &SmearerS};
  std::vector<int> mask_types = {2, 1};
  const int Nstep_per_type = 2 * Nd;
  const int Nsmear = Nstep_per_type * (int)mask_types.size();

  SmearedConfigurationRect<Gimpl> SmartConfig(UGrid, Nsmear, Stouts, mask_types);
  SmartConfig.set_Field(U);

  const RealD tolF = 1.0e-24; // on |diff|^2/|ref|^2
  const RealD tolS = 1.0e-12; // on |diff|/|ref|
  int fail = 0;

  LatticeGaugeField Fold(UGrid), Fdef(UGrid), Fdiff(UGrid);

  ////////////////////////////////////////////////////////////////////
  // Per-level comparison: force and ln det
  ////////////////////////////////////////////////////////////////////
  for (int smr = Nsmear - 1; smr >= 0; smr--) {
    const LatticeGaugeField &Uin = (smr > 0) ? SmartConfig.get_smeared_conf(smr - 1) : U;
    int knl = mask_types[smr / Nstep_per_type];

    Fold = Zero();
    SmartConfig.logDetJacobianForceLevel(0, Uin, Fold, smr);
    Fdef = Zero();
    SmartConfig.logDetJacobianForceLevel(Uin, Fdef, smr);

    Fdiff = Fdef - Fold;
    RealD n2d = norm2(Fdiff), n2r = norm2(Fold);
    bool okF = (n2d <= tolF * n2r);
    if (!okF) fail++;
    std::cout << GridLogMessage << "LEVELCHECK force smr=" << smr << " kernel=" << knl
              << " |Fold|^2 = " << n2r << "  |Fdef-Fold|^2 = " << n2d
              << "  rel = " << ((n2r > 0.0) ? n2d / n2r : 0.0)
              << (okF ? "  ok" : "  FAIL") << std::endl;

    RealD Sold = SmartConfig.logDetJacobianLevel(0, Uin, smr);
    RealD Sdef = SmartConfig.logDetJacobianLevel(Uin, smr);
    RealD dS = fabs(Sdef - Sold), aS = fabs(Sold);
    bool okS = (dS <= tolS * aS);
    if (!okS) fail++;
    std::cout << GridLogMessage << "LEVELCHECK lndet smr=" << smr << " kernel=" << knl
              << " lndet_old = " << Sold << "  lndet_def = " << Sdef
              << "  rel = " << ((aS > 0.0) ? dS / aS : 0.0)
              << (okS ? "  ok" : "  FAIL") << std::endl;
  }

  ////////////////////////////////////////////////////////////////////
  // End-to-end: full chain rule and level sum
  ////////////////////////////////////////////////////////////////////
  SmartConfig.logDetJacobianForce(0, Fold);
  SmartConfig.logDetJacobianForce(Fdef);
  Fdiff = Fdef - Fold;
  RealD n2d = norm2(Fdiff), n2r = norm2(Fold);
  bool okF = (n2d <= tolF * n2r);
  if (!okF) fail++;
  std::cout << GridLogMessage << "TOTALCHECK force |Fold|^2 = " << n2r
            << "  |Fdef-Fold|^2 = " << n2d
            << "  rel = " << ((n2r > 0.0) ? n2d / n2r : 0.0)
            << (okF ? "  ok" : "  FAIL") << std::endl;

  RealD Sold = SmartConfig.logDetJacobian(0);
  RealD Sdef = SmartConfig.logDetJacobian();
  RealD dS = fabs(Sdef - Sold), aS = fabs(Sold);
  bool okS = (dS <= tolS * aS);
  if (!okS) fail++;
  std::cout << GridLogMessage << "TOTALCHECK lndet old = " << Sold
            << "  def = " << Sdef
            << "  rel = " << ((aS > 0.0) ? dS / aS : 0.0)
            << (okS ? "  ok" : "  FAIL") << std::endl;

  if (fail) std::cout << GridLogMessage << "FAIL (" << fail << " checks)" << std::endl;
  else      std::cout << GridLogMessage << "PASS" << std::endl;

  Grid_finalize();
  return fail ? 1 : 0;
}
