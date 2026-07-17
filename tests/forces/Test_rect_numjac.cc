/*
  Driver 5: ABSOLUTE check of logDetJacobianLevel against a brute-force
  numerical Jacobian of the actual level-0 smearing map (see README.md).

  The masking makes the active-active block of dU'/dU site-diagonal (the
  staple of an active link contains only inactive same-direction links), so

     lndet(true) = sum_{active y} log det[ J8(y) ],
     J8(y)_{ba}  = d omega^b(y) / d omega^a(y),

  where the input perturbation is U_mu(y) -> exp(i eps t^a) U_mu(y) and the
  output coordinate is omega^b = -2i tr( t^b dU' U'^dag ), central
  differences. The map is the production one: set_Field -> get_smeared_conf(0).

  Run for mask_types={1} (plq, CONTROL: validates the extraction against the
  FD-consistent plq implementation) and {2} (rect, the question):
  compare lndet_num with logDetJacobianLevel(U,0). Ratio ~1 => action right
  (force carries the x2). Ratio ~2 => action halved.
*/
#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#include <Grid/qcd/smearing/GaugeConfigurationRect.h>

using namespace Grid;

typedef PeriodicGimplR Gimpl;

// level-0 mask (mu=0, cb=0), same formulas as the class constructor
static bool active_site(const Coordinate &x, int mask_type)
{
  int s = 0;
  if (mask_type == 1) {
    for (int d = 0; d < Nd; d++) s += x[d];
  } else {
    s = x[0];
    for (int d = 1; d < Nd; d++) s += x[d] / 2;
  }
  return (s % 2) == 0;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  std::cout << std::setprecision(14);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd, vComplex::Nsimd()),
                                                        GridDefaultMpi());
  Coordinate latt = GridDefaultLatt();

  std::vector<int> seeds({1, 2, 3, 5});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds);

  LatticeGaugeField U0(UGrid);
  SU<Nc>::HotConfiguration(RNG4, U0);

  RealD rho = 0.1;
  Smear_Stout<Gimpl> SmearerS(rho);
  Rect_Stout<Gimpl>  SmearerR(0.0, rho, 0.0);
  const int Nsmear = 2 * Nd;
  const int mu = 0;          // level 0
  const RealD eps = 1.0e-4;  // central-difference step
  const int Ngen = 8;

  ColourMatrix ta[Ngen];
  for (int a = 0; a < Ngen; a++) SU<Nc>::generator(a, ta[a]);
  Complex ci(0.0, 1.0);

  int fail = 0;
  for (int mask_type : {1, 2}) {
    std::vector<Smear_Stout<Gimpl> *> Stouts = {(mask_type == 1) ? &SmearerS : (Smear_Stout<Gimpl> *)&SmearerR};
    std::vector<int> mts = {mask_type};
    SmearedConfigurationRect<Gimpl> Cfg(UGrid, Nsmear, Stouts, mts);

    // baseline smeared level-0 config
    LatticeGaugeField Uwork(UGrid);
    Uwork = U0;
    Cfg.set_Field(Uwork);
    LatticeGaugeField Up0(UGrid);
    Up0 = Cfg.get_smeared_conf(0);

    RealD lndet_impl = Cfg.logDetJacobianLevel(U0, 0);

    // enumerate sites
    RealD lndet_num = 0.0;
    int nactive = 0;
    Coordinate x(Nd);
    for (x[0] = 0; x[0] < latt[0]; x[0]++)
    for (x[1] = 0; x[1] < latt[1]; x[1]++)
    for (x[2] = 0; x[2] < latt[2]; x[2]++)
    for (x[3] = 0; x[3] < latt[3]; x[3]++) {
      if (!active_site(x, mask_type)) continue;
      nactive++;

      typedef typename LatticeGaugeField::scalar_object LorentzCM;
      LorentzCM Usite, UPsite, Us2, UPs2;
      peekSite(Usite, U0, x);
      peekSite(UPsite, Up0, x);
      ColourMatrix Uy  = peekIndex<LorentzIndex>(Usite, mu);
      ColourMatrix UPy = peekIndex<LorentzIndex>(UPsite, mu);

      Eigen::MatrixXd J(Ngen, Ngen);
      for (int a = 0; a < Ngen; a++) {
        ColourMatrix Upl[2];
        for (int s = 0; s < 2; s++) {
          RealD sgn = (s == 0) ? +1.0 : -1.0;
          // exp(i sgn eps t^a) to O(eps^4)
          ColourMatrix E = ci * (sgn * eps) * ta[a];
          ColourMatrix P;
          P = ComplexD(1.0);
          P = P + E + 0.5 * E * E + (1.0 / 6.0) * E * E * E;

          Us2 = Usite;
          pokeIndex<LorentzIndex>(Us2, (ColourMatrix)(P * Uy), mu);
          Uwork = U0;
          pokeSite(Us2, Uwork, x);
          Cfg.set_Field(Uwork);
          peekSite(UPs2, Cfg.get_smeared_conf(0), x);
          Upl[s] = peekIndex<LorentzIndex>(UPs2, mu);
        }
        ColourMatrix dUp = (0.5 / eps) * (Upl[0] - Upl[1]);
        ColourMatrix Pm = dUp * adj(UPy);
        for (int b = 0; b < Ngen; b++) {
          ComplexD tr = TensorRemove(trace(ta[b] * Pm));
          J(b, a) = real(ComplexD(0.0, -2.0) * tr);
        }
      }
      lndet_num += std::log(J.determinant());
    }

    bool ok = (fabs(lndet_num / lndet_impl - 1.0) < 1.0e-6);
    if (!ok) fail++;
    std::cout << GridLogMessage << "NUMJAC mask_type=" << mask_type
              << " active=" << nactive
              << " lndet_num=" << lndet_num
              << " lndet_impl=" << lndet_impl
              << " ratio num/impl=" << lndet_num / lndet_impl
              << (ok ? "  ok" : "  FAIL") << std::endl;
  }

  if (fail) std::cout << GridLogMessage << "FAIL (" << fail << " kernels)" << std::endl;
  else      std::cout << GridLogMessage << "PASS" << std::endl;

  Grid_finalize();
  return fail ? 1 : 0;
}
