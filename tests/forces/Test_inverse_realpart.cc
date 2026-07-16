/*
  Unit probe: Inverse vs Inverse_RealPart on an explicitly REAL, well
  conditioned 8x8 adjoint-like matrix field. Decides whether the observed
  force-level discrepancy comes from imaginary content in MpAd or from the
  Inverse_RealPart routine itself on this build.
*/
#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd, vComplex::Nsimd()),
                                                        GridDefaultMpi());
  GridParallelRNG RNG(UGrid);
  RNG.SeedFixedIntegers({1, 2, 3, 5});

  const int N = 8;
  typedef Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > AdjMatField;

  AdjMatField M(UGrid), Mre(UGrid), One(UGrid);
  gaussian(RNG, M);
  Mre = 0.5 * (M + conjugate(M));      // entrywise real
  M   = ComplexD(1.0, 0.0);
  M   = M + 0.1 * Mre;                 // well conditioned, exactly real
  One = ComplexD(1.0, 0.0);

  AdjMatField Minv_c(UGrid), Minv_r(UGrid), R(UGrid);

  Minv_c = Inverse(M);
  Minv_r = Inverse_RealPart(M);

  R = M * Minv_c - One;
  RealD res_c = norm2(R);
  R = M * Minv_r - One;
  RealD res_r = norm2(R);
  R = Minv_c - Minv_r;
  RealD diff = norm2(R);
  RealD nrm = norm2(Minv_c);

  std::cout << GridLogMessage << "|M Minv_complex - 1|^2  = " << res_c << std::endl;
  std::cout << GridLogMessage << "|M Minv_realpart - 1|^2 = " << res_r << std::endl;
  std::cout << GridLogMessage << "|Minv_c - Minv_r|^2 = " << diff
            << "  rel = " << diff / nrm << std::endl;

  if (res_r < 1e-20 && diff / nrm < 1e-24)
    std::cout << GridLogMessage << "PASS: Inverse_RealPart correct on real input" << std::endl;
  else
    std::cout << GridLogMessage << "FAIL: Inverse_RealPart deviates on real input" << std::endl;

  Grid_finalize();
  return 0;
}
