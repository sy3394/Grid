/*
  Driver 3: cross-class check for the plq kernel (anchor B, see README.md).

  SmearedConfigurationMasked and SmearedConfigurationRect with mask_types={1}
  implement the same plq FTHMC with different internal normalisations
  (Masked: ta=2i, dJdX -1.0, final -0.5; Rect: ta=i, dJdX -0.5, final -1.0).
  The output force and lndet are algebra-basis independent, so they must
  agree numerically — this validates the whole convention ledger and the
  rect class's plumbing against the independently written Masked class.

  Compares per level: smeared configs, force, lndet; and the totals.
*/
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/GaugeConfigurationRect.h>

using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  std::cout << std::setprecision(14);

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
  const int Nsmear = 2 * Nd;

  SmearedConfigurationMasked<Gimpl> CfgM(UGrid, Nsmear, SmearerS);
  std::vector<Smear_Stout<Gimpl> *> Stouts = {&SmearerS};
  std::vector<int> mask_types = {1};
  SmearedConfigurationRect<Gimpl> CfgR(UGrid, Nsmear, Stouts, mask_types);

  CfgM.set_Field(U);
  CfgR.set_Field(U);

  const RealD tolC = 1.0e-20; // smeared configs: identical maps expected
  const RealD tolF = 1.0e-18; // force/lndet: different internal arithmetic
  int fail = 0;

  LatticeGaugeField FM(UGrid), FR(UGrid), D(UGrid);

  for (int smr = Nsmear - 1; smr >= 0; smr--) {
    // smeared configuration at this level
    D = CfgM.get_smeared_conf(smr) - CfgR.get_smeared_conf(smr);
    RealD nC = norm2(D), rC = norm2(CfgM.get_smeared_conf(smr));
    bool okC = (nC <= tolC * rC);
    if (!okC) fail++;

    const LatticeGaugeField &UM = (smr > 0) ? CfgM.get_smeared_conf(smr - 1) : U;
    const LatticeGaugeField &UR = (smr > 0) ? CfgR.get_smeared_conf(smr - 1) : U;

    FM = Zero(); CfgM.logDetJacobianForceLevel(UM, FM, smr);
    FR = Zero(); CfgR.logDetJacobianForceLevel(UR, FR, smr);
    D = FM - FR;
    RealD nF = norm2(D), rF = norm2(FM);
    bool okF = (nF <= tolF * rF);
    if (!okF) fail++;

    RealD SM = CfgM.logDetJacobianLevel(UM, smr);
    RealD SR = CfgR.logDetJacobianLevel(UR, smr);
    RealD dS = fabs(SM - SR), aS = fabs(SM);
    bool okS = (dS <= 1.0e-9 * aS);
    if (!okS) fail++;

    std::cout << GridLogMessage << "XCHECK smr=" << smr
              << " |dU|^2/|U|^2=" << nC / rC
              << " |dF|^2/|F|^2=" << ((rF > 0) ? nF / rF : 0.0)
              << " lndet M=" << SM << " R=" << SR << " rel=" << ((aS > 0) ? dS / aS : 0.0)
              << ((okC && okF && okS) ? "  ok" : "  FAIL") << std::endl;
  }

  CfgM.logDetJacobianForce(FM);
  CfgR.logDetJacobianForce(FR);
  D = FM - FR;
  RealD nF = norm2(D), rF = norm2(FM);
  bool okF = (nF <= tolF * rF);
  if (!okF) fail++;
  if (!okF && fabs(nF / rF - 0.25) < 1.0e-6)
    std::cout << GridLogMessage << "HINT: |dF|^2/|F|^2 = 1/4 exactly => the force-normalisation "
              << "fix (-0.5) is applied on ONE side only (Rect fixed, Masked not, or vice versa); "
              << "see FIXME EXTRA-FACTOR-2 in GaugeConfigurationRect.h" << std::endl;
  RealD SM = CfgM.logDetJacobian();
  RealD SR = CfgR.logDetJacobian();
  bool okS = (fabs(SM - SR) <= 1.0e-9 * fabs(SM));
  if (!okS) fail++;
  std::cout << GridLogMessage << "XCHECK TOTAL |dF|^2/|F|^2=" << nF / rF
            << " lndet M=" << SM << " R=" << SR
            << ((okF && okS) ? "  ok" : "  FAIL") << std::endl;

  if (fail) std::cout << GridLogMessage << "FAIL (" << fail << ")" << std::endl;
  else      std::cout << GridLogMessage << "PASS" << std::endl;

  Grid_finalize();
  return fail ? 1 : 0;
}
