/*
  Driver 5b: FD check of a SINGLE level (no chain rule): is
  logDetJacobianForceLevel(U,F,0) the gradient of logDetJacobianLevel(U,0)?

  Isolates the level routine from the AnalyticSmearedForce chain. Run for
  plq and rect schedules; eps sweep prints dS/dSpred (-> 1 if consistent,
  -> 1/2 if the level force is 2x).

  NOTE dSpred convention: the level force F is the same UdSdU-like object
  the total force assembles from (before the final Ta in the top level),
  so dSpred = -sum tr(P F) eps 2 HMC_MOMENTUM_DENOMINATOR as in ForceTest.
*/
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/GaugeConfigurationRect.h>

using namespace Grid;

typedef PeriodicGimplR Gimpl;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  std::cout << std::setprecision(14);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd, vComplex::Nsimd()),
                                                        GridDefaultMpi());
  std::vector<int> seeds({1, 2, 3, 5});
  GridSerialRNG sRNG;          sRNG.SeedFixedIntegers(seeds);
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds);

  LatticeGaugeField U0(UGrid), P(UGrid), U(UGrid), F(UGrid);
  LatticeColourMatrix Pmu(UGrid);
  SU<Nc>::HotConfiguration(RNG4, U0);
  Gimpl::generate_momenta(P, sRNG, RNG4);

  RealD rho = 0.1;
  Smear_Stout<Gimpl> SmearerS(rho);
  Rect_Stout<Gimpl>  SmearerR(0.0, rho, 0.0);
  const int Nsmear = 2 * Nd;

  for (int mask_type : {1, 2}) {
    std::vector<Smear_Stout<Gimpl> *> Stouts = {(mask_type == 1) ? &SmearerS : (Smear_Stout<Gimpl> *)&SmearerR};
    std::vector<int> mts = {mask_type};
    SmearedConfigurationRect<Gimpl> Cfg(UGrid, Nsmear, Stouts, mts);

    // S = -lndet to match the JacobianAction sign; F likewise
    F = Zero();
    Cfg.logDetJacobianForceLevel(U0, F, 0);
    F = Ta(F);

    for (RealD eps : {0.01, 0.005, 0.0025}) {
      RealD S1 = -Cfg.logDetJacobianLevel(U0, 0);
      U = U0;
      LatticeGaugeField Pc(UGrid); Pc = P;
      Gimpl::update_field(Pc, U, eps);
      RealD S2 = -Cfg.logDetJacobianLevel(U, 0);

      LatticeComplex dSl(UGrid); dSl = Zero();
      for (int mu = 0; mu < Nd; mu++) {
        auto Fmu = PeekIndex<LorentzIndex>(F, mu);
        Pmu = PeekIndex<LorentzIndex>(P, mu);
        dSl = dSl + trace(Pmu * Fmu) * eps * 2.0 * HMC_MOMENTUM_DENOMINATOR;
      }
      ComplexD dSpred = sum(dSl);
      std::cout << GridLogMessage << "FDLVL0 mask_type=" << mask_type
                << " eps=" << eps
                << " dS=" << S2 - S1
                << " dSpred=" << dSpred.real()
                << " dS/dSpred=" << (S2 - S1) / dSpred.real() << std::endl;
    }

    // TOTAL (all levels + chain), small-eps sweep
    LatticeGaugeField Ftot(UGrid);
    Cfg.set_Field(U0);
    Cfg.logDetJacobianForce(Ftot);
    Ftot = Ta(Ftot);
    for (RealD eps : {0.01, 0.005, 0.0025, 0.00125}) {
      Cfg.set_Field(U0);
      RealD S1 = -Cfg.logDetJacobian();
      U = U0;
      LatticeGaugeField Pc(UGrid); Pc = P;
      Gimpl::update_field(Pc, U, eps);
      Cfg.set_Field(U);
      RealD S2 = -Cfg.logDetJacobian();

      LatticeComplex dSl(UGrid); dSl = Zero();
      for (int mu = 0; mu < Nd; mu++) {
        auto Fmu = PeekIndex<LorentzIndex>(Ftot, mu);
        Pmu = PeekIndex<LorentzIndex>(P, mu);
        dSl = dSl + trace(Pmu * Fmu) * eps * 2.0 * HMC_MOMENTUM_DENOMINATOR;
      }
      ComplexD dSpred = sum(dSl);
      std::cout << GridLogMessage << "FDTOT mask_type=" << mask_type
                << " eps=" << eps
                << " dS=" << S2 - S1
                << " dSpred=" << dSpred.real()
                << " dS/dSpred=" << (S2 - S1) / dSpred.real() << std::endl;
    }
  }

  Grid_finalize();
  return 0;
}
