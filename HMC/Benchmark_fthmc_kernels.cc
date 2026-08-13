/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./HMC/Benchmark_fthmc_kernels.cc

    Copyright (C) 2026

Author: Shuhei Yamamoto <shuhei.yamamoto.2011@gmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */

// FTHMC kernel micro-benchmark: optimised vs non-optimised rect-FTHMC
// paths, ONE binary (Lattice 2026 "Code Optimization" backup slide).
//
// SmearedConfigurationRect keeps the pre-optimisation reference
// implementations as `int old` overloads next to the optimised defaults,
// so a single binary times both on identical state, interleaved
// opt,old,opt,old per sample:
//   lndet    : logDetJacobian()        vs logDetJacobianForce-era logDetJacobian(0)
//   jacforce : logDetJacobianForce(f)  vs logDetJacobianForce(0,f)
// Components with no old counterpart are timed single-path, for context:
//   forward   : F, set_Field / fill_smearedSet (plain smearers; untouched)
//   inverse   : F^-1, fixed-point inversion implemented here (see below);
//               iteration counts reported separately
//   xformforce: gauge-force pullback through F (smeared_force /
//               AnalyticSmearedForce; no old overload exists)
//   rawderiv  : raw gauge deriv on the smeared field — NOT FTHMC cost.
//
// REQUIREMENT: DEBUG must stay #undef in GaugeConfigurationRect.h (its
// default). With DEBUG defined the default routines self-check against
// the old ones, so "opt" timings would include old-path work.  Enforced
// below at compile time.
//
// Kernels are named by strings over {p,r}: p = plaquette (masked stout,
// mask type 1), r = short-side rectangle Rs (mask type 2).  Two-letter
// names are sequential compositions; the LEFTMOST letter is applied FIRST
// to the thin field (mask_types[0] in SmearedConfigurationRect).  Each
// letter is one full masked step = 2*Nd = 8 sub-levels.
//
// Usage (defaults in brackets):
//   Benchmark_fthmc_kernels --grid 24.24.24.40 --mpi ...
//     [--kernels p,r,pp,pr,rp,rr]   comma list; prefer ONE kernel per run
//                                   (each SmearedConfigurationRect instance
//                                   holds/leaks its full smeared tower)
//     [--rho 0.12]                  smearing parameter, both kernels
//     [--nwarm 5] [--nmeas 20]      warm-up / timed calls per component
//     [--beta 6.0]                  Wilson gauge action for the force pullback
//     [--fptol 1e-13] [--fpmaxit 50] fixed-point stopping (relative L2 of
//                                   successive iterates) and iteration cap
//     [--config file.scidac]        fixed gauge configuration (SciDAC);
//     [--cformat scidac|nersc]      reader selection [scidac]
//   Without --config: fixed-seed SU(3) hot start.
//
// NOTE: Aurora-GPU-written SciDAC files carry wrong stored checksums; the
// reader tolerates this only with commit c5fb9a59 (IldgIO.h: downgrade the
// checksum assert to a warning, keep the GRID_FIELD_NORM hard check).
// Cherry-pick it if --config points at an Aurora-written file.
//
// Machine-readable output: lines beginning "FTHMCBENCH" (boss rank only),
//   FTHMCBENCH kernel=<k> comp=<c> impl=<opt|old|single> n=<n> median_ms=...
//              q16_ms=... q84_ms=... min_ms=... max_ms=... [fp_iters_...]
// plus a per-kernel summary table.  The companion PBS script
// (systems/Aurora2/fthmc_kernel_bench.pbs) aggregates over repeats.

#include <Grid/Grid.h>
#include <Grid/qcd/smearing/GaugeConfigurationMasked.h>
#ifdef DEBUG
#error "DEBUG is defined after including GaugeConfigurationMasked.h: its debug blocks would pollute the timings. Keep '//#define DEBUG' commented there and rebuild."
#endif
#include <Grid/qcd/smearing/GaugeConfigurationRect.h>
#ifdef DEBUG
#error "GaugeConfigurationRect.h has DEBUG defined: the default routines would self-check against the old ones and the opt timings would be meaningless. Restore '#undef DEBUG' there and rebuild."
#endif

using namespace std;
using namespace Grid;

typedef PeriodicGimplR Gimpl;
typedef Gimpl::GaugeField     GaugeField;
typedef Gimpl::GaugeLinkField GaugeLinkField;

//////////////////////////////////////////////////////////////////////
// Timing: device + communicator sync around every sample
//////////////////////////////////////////////////////////////////////
static double now_sync(GridBase *grid)
{
  accelerator_barrier();
  grid->Barrier();
  return usecond();
}

template <class Work>
static std::vector<double> timeComponent(GridBase *grid, int nwarm, int nmeas, Work &&work)
{
  for (int i = 0; i < nwarm; i++) work();
  std::vector<double> ms(nmeas);
  for (int i = 0; i < nmeas; i++) {
    double t0 = now_sync(grid);
    work();
    double t1 = now_sync(grid);
    ms[i] = (t1 - t0) / 1.0e3;
  }
  return ms;
}

// opt,old,opt,old per sample: both paths see the same thermal/clock state
template <class WOpt, class WOld>
static void timeInterleaved(GridBase *grid, int nwarm, int nmeas,
                            WOpt &&wopt, WOld &&wold,
                            std::vector<double> &opt_ms, std::vector<double> &old_ms)
{
  for (int i = 0; i < nwarm + nmeas; i++) {
    double t0 = now_sync(grid);
    wopt();
    double t1 = now_sync(grid);
    wold();
    double t2 = now_sync(grid);
    if (i >= nwarm) {
      opt_ms.push_back((t1 - t0) / 1.0e3);
      old_ms.push_back((t2 - t1) / 1.0e3);
    }
  }
}

struct Stats { double med, q16, q84, mn, mx; int n; };

static Stats stats(std::vector<double> v)
{
  std::sort(v.begin(), v.end());
  int n = v.size();
  auto pct = [&](double p) { int i = (int)(p * (n - 1) + 0.5); return v[i]; };
  Stats s;
  s.n   = n;
  s.mn  = v.front();
  s.mx  = v.back();
  s.q16 = pct(0.16);
  s.q84 = pct(0.84);
  s.med = (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  return s;
}

static void report(const std::string &kernel, const std::string &comp,
                   const std::string &impl, const Stats &s,
                   const std::string &extra = "")
{
  std::cout << GridLogMessage << "FTHMCBENCH"
            << " kernel=" << kernel
            << " comp=" << comp
            << " impl=" << impl
            << " n=" << s.n
            << std::fixed << std::setprecision(3)
            << " median_ms=" << s.med
            << " q16_ms=" << s.q16
            << " q84_ms=" << s.q84
            << " min_ms=" << s.mn
            << " max_ms=" << s.mx
            << (extra.empty() ? "" : " " + extra)
            << std::defaultfloat << std::endl;
}

//////////////////////////////////////////////////////////////////////
// Per-level masks, transcribed from the SmearedConfigurationRect
// constructor (they are private there).  Level l of a tower:
//   block i = l/(2*Nd) selects mask_types[i]; j = l%(2*Nd);
//   mu = (j/2)%Nd; cb = j%2.  Only direction mu is masked.
// The driver cross-checks the transcription: F^-1(F(U)) is compared
// against U before any timing (see fp_relerr).
//////////////////////////////////////////////////////////////////////
static void levelMask(LatticeComplex &mask, int &mu,
                      GridCartesian *UGrid, GridRedBlackCartesian *UrbGrid,
                      int mask_type, int j)
{
  LatticeComplex ones(UGrid);  ones = ComplexD(1.0, 0.0);
  LatticeComplex zeros(UGrid); zeros = Zero();
  mu = (j / 2) % Nd;
  int cb = j % 2;
  mask = Zero();
  switch (mask_type) {
  case 1: {
    LatticeComplex tmpcb(UrbGrid);
    pickCheckerboard(cb, tmpcb, ones);
    setCheckerboard(mask, tmpcb);
    break;
  }
  case 2: {
    LatticeInteger coor_nu(UGrid), coor_sum(UGrid);
    coor_sum = Zero();
    for (int nu = 0; nu < Nd; nu++) {
      LatticeCoordinate(coor_nu, nu);
      if (nu != mu) coor_nu = div(coor_nu, 2);
      coor_sum = coor_sum + coor_nu;
    }
    mask = where(mod(coor_sum, 2) == (Integer)cb, ones, zeros);
    break;
  }
  default:
    assert(0 && "mask_type must be 1 (plq) or 2 (Rs)");
  }
}

//////////////////////////////////////////////////////////////////////
// Inverse of one masked stout step.
//
// The step is V = (1-m) U + m exp(iQ(U)) U on direction mu, where the
// staple C entering Q is built exclusively from frozen (unmasked) links
// -- that is what the masks are for -- but Q = Ta(C U^dag) still depends
// on the updated link itself.  Solve by fixed-point iteration
//     U_0 = V,   U_{k+1} = exp(-iQ(U_k)) V   on masked links,
// using exp(-iQ(U_k)) = U_k W_k^dag with W_k = Stout(U_k), i.e. the
// class's own smearer supplies Q; no staple code is duplicated here.
// Returns the iteration count; stops when the relative L2 distance of
// successive iterates drops below tol.
//////////////////////////////////////////////////////////////////////
static int invertLevel(GaugeField &U, const GaugeField &V,
                       Smear_Stout<Gimpl> *stout,
                       const LatticeComplex &mask, const LatticeComplex &cmask,
                       int mu, double tol, int maxit)
{
  GridBase *grid = V.Grid();
  GaugeField     W(grid);
  GaugeLinkField Umu(grid), Wmu(grid), Vmu(grid), Unew(grid);

  U   = V;
  Vmu = PeekIndex<LorentzIndex>(V, mu);
  RealD vnorm = norm2(Vmu);

  for (int k = 0; k < maxit; k++) {
    stout->smear(W, U);
    Umu  = PeekIndex<LorentzIndex>(U, mu);
    Wmu  = PeekIndex<LorentzIndex>(W, mu);
    Unew = Umu * adj(Wmu) * Vmu;          // exp(-iQ(U_k)) V on every site
    Unew = Unew * mask + Vmu * cmask;     // frozen links stay = V
    RealD diff = norm2(Unew - Umu);
    PokeIndex<LorentzIndex>(U, Unew, mu);
    if (diff <= tol * tol * vnorm) return k + 1;
  }
  std::cout << GridLogMessage << "FTHMCBENCH WARNING: fixed point not converged in "
            << maxit << " iterations (mu=" << mu << ")" << std::endl;
  return maxit;
}

// Full inverse map: peel the levels off in reverse order.
static void invertMap(GaugeField &U, const GaugeField &Vtop,
                      const std::vector<Smear_Stout<Gimpl> *> &stouts,
                      const std::vector<LatticeComplex> &masks,
                      const std::vector<LatticeComplex> &cmasks,
                      const std::vector<int> &mus,
                      double tol, int maxit, std::vector<int> &iters)
{
  int Nsmear = masks.size();
  GaugeField V(Vtop.Grid());
  V = Vtop;
  for (int l = Nsmear - 1; l >= 0; l--) {
    iters[l] = invertLevel(U, V, stouts[l / (2 * Nd)], masks[l], cmasks[l],
                           mus[l], tol, maxit);
    if (l > 0) V = U;
  }
}

//////////////////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  // -------- options
  std::string kernels = "p,r,pp,pr,rp,rr";
  std::string config  = "";
  std::string cformat = "scidac";
  double rho   = 0.12;
  double beta  = 6.0;
  double fptol = 1.0e-13;
  int fpmaxit  = 50;
  int nwarm    = 5;
  int nmeas    = 20;

  if (GridCmdOptionExists(argv, argv + argc, "--kernels")) kernels = GridCmdOptionPayload(argv, argv + argc, "--kernels");
  if (GridCmdOptionExists(argv, argv + argc, "--config"))  config  = GridCmdOptionPayload(argv, argv + argc, "--config");
  if (GridCmdOptionExists(argv, argv + argc, "--cformat")) cformat = GridCmdOptionPayload(argv, argv + argc, "--cformat");
  if (GridCmdOptionExists(argv, argv + argc, "--rho"))     rho     = std::stod(GridCmdOptionPayload(argv, argv + argc, "--rho"));
  if (GridCmdOptionExists(argv, argv + argc, "--beta"))    beta    = std::stod(GridCmdOptionPayload(argv, argv + argc, "--beta"));
  if (GridCmdOptionExists(argv, argv + argc, "--fptol"))   fptol   = std::stod(GridCmdOptionPayload(argv, argv + argc, "--fptol"));
  if (GridCmdOptionExists(argv, argv + argc, "--fpmaxit")) fpmaxit = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--fpmaxit"));
  if (GridCmdOptionExists(argv, argv + argc, "--nwarm"))   nwarm   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--nwarm"));
  if (GridCmdOptionExists(argv, argv + argc, "--nmeas"))   nmeas   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--nmeas"));

  GridCartesian *UGrid =
    SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                   GridDefaultSimd(Nd, vComplex::Nsimd()),
                                   GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::cout << GridLogMessage << "FTHMCBENCH setup"
            << " grid=" << GridDefaultLatt()
            << " kernels=" << kernels
            << " rho=" << rho
            << " beta=" << beta
            << " nwarm=" << nwarm << " nmeas=" << nmeas
            << " fptol=" << fptol << " fpmaxit=" << fpmaxit
            << " config=" << (config.empty() ? "<hot start>" : config)
            << std::endl;

  // -------- fixed gauge configuration
  LatticeGaugeField U(UGrid);
  if (config.empty()) {
    std::vector<int> seeds({1, 2, 3, 4, 5, 6, 7, 8});
    GridParallelRNG RNG4(UGrid);
    RNG4.SeedFixedIntegers(seeds);
    SU<Nc>::HotConfiguration(RNG4, U);
    std::cout << GridLogMessage << "FTHMCBENCH hot start, fixed seeds 1..8" << std::endl;
  } else if (cformat == "nersc") {
    FieldMetaData header;
    NerscIO::readConfiguration(U, header, config);
  } else {
    emptyUserRecord record;
    ScidacReader RD;
    RD.open(config);
    RD.readScidacFieldRecord(U, record);
    RD.close();
  }
  std::cout << GridLogMessage << "FTHMCBENCH input plaquette "
            << WilsonLoops<Gimpl>::avgPlaquette(U) << std::endl;

  WilsonGaugeActionR GaugeAction(beta);
  GaugeAction.is_smeared = true;

  // -------- kernel loop
  std::stringstream ss(kernels);
  std::string kname;
  while (std::getline(ss, kname, ',')) {
    for (char c : kname) assert((c == 'p' || c == 'r') && "kernel letters must be p or r");
    int nsteps = kname.size();
    int Nsmear = nsteps * 2 * Nd;

    std::cout << GridLogMessage << "==============================================" << std::endl;
    std::cout << GridLogMessage << "FTHMCBENCH kernel " << kname
              << " (" << nsteps << " step(s), " << Nsmear << " sub-levels)" << std::endl;
    std::cout << GridLogMessage << "==============================================" << std::endl;

    std::vector<Smear_Stout<Gimpl> *> stouts;
    std::vector<int> mask_types;
    for (char c : kname) {
      if (c == 'p') { stouts.push_back(new Smear_Stout<Gimpl>(rho));          mask_types.push_back(1); }
      else          { stouts.push_back(new Rect_Stout<Gimpl>(0.0, rho, 0.0)); mask_types.push_back(2); }
    }

    // masks for the driver-side inverse
    std::vector<LatticeComplex> masks, cmasks;
    std::vector<int> mus(Nsmear);
    LatticeComplex ones(UGrid); ones = ComplexD(1.0, 0.0);
    for (int l = 0; l < Nsmear; l++) {
      masks.emplace_back(UGrid);
      cmasks.emplace_back(UGrid);
      levelMask(masks[l], mus[l], UGrid, UrbGrid, mask_types[l / (2 * Nd)], l % (2 * Nd));
      cmasks[l] = ones - masks[l];
    }

    Stats sF, sI, sLo, sLl, sJo, sJl, sR, sX;   // o=opt, l=old
    int itot = 0, imax = 0;

    {
      SmearedConfigurationRect<Gimpl> smart(UGrid, Nsmear, stouts, mask_types);

      // ---- forward map F (single-path: fill_smearedSet uses the plain smearers)
      sF = stats(timeComponent(UGrid, nwarm, nmeas, [&] { smart.set_Field(U); }));
      report(kname, "forward", "single", sF);

      LatticeGaugeField V(UGrid);
      V = smart.get_U(true);

      // ---- inverse map F^-1 (single-path, driver-side):
      //      correctness + iteration counts first (untimed)
      LatticeGaugeField Uinv(UGrid);
      std::vector<int> iters(Nsmear);
      invertMap(Uinv, V, stouts, masks, cmasks, mus, fptol, fpmaxit, iters);
      for (int l = 0; l < Nsmear; l++) { itot += iters[l]; imax = std::max(imax, iters[l]); }
      RealD fprelerr = std::sqrt(norm2(Uinv - U) / norm2(U));
      std::cout << GridLogMessage << "FTHMCBENCH inverse check: |F^-1(F(U))-U|/|U| = "
                << std::scientific << fprelerr << std::defaultfloat
                << "  fixed-point iterations per level:";
      for (int l = 0; l < Nsmear; l++) std::cout << " " << iters[l];
      std::cout << std::endl;

      sI = stats(timeComponent(UGrid, nwarm, nmeas, [&] {
        invertMap(Uinv, V, stouts, masks, cmasks, mus, fptol, fpmaxit, iters);
      }));
      {
        std::stringstream ex;
        ex << "fp_iters_total=" << itot << " fp_iters_max=" << imax
           << " fp_relerr=" << std::scientific << std::setprecision(2) << fprelerr;
        report(kname, "inverse", "single", sI, ex.str());
      }

      // ---- ln det J: optimised default vs old reference, interleaved.
      //      Correctness first (untimed).
      RealD lndet_opt = smart.logDetJacobian();
      RealD lndet_old = smart.logDetJacobian(0);
      std::cout << GridLogMessage << "FTHMCBENCH lndet check: opt " << lndet_opt
                << " old " << lndet_old << " diff " << lndet_opt - lndet_old << std::endl;
      {
        std::vector<double> oms, lms;
        RealD ld;
        timeInterleaved(UGrid, nwarm, nmeas,
                        [&] { ld = smart.logDetJacobian(); },
                        [&] { ld = smart.logDetJacobian(0); }, oms, lms);
        sLo = stats(oms); sLl = stats(lms);
        report(kname, "lndet", "opt", sLo);
        report(kname, "lndet", "old", sLl);
      }

      // ---- Jacobian (log-det) force: optimised default vs old reference,
      //      interleaved.  Correctness first (untimed).
      LatticeGaugeField dJ(UGrid), dJold(UGrid);
      smart.logDetJacobianForce(dJ);
      smart.logDetJacobianForce(0, dJold);
      std::cout << GridLogMessage << "FTHMCBENCH jacforce check: |opt-old|/|old| = "
                << std::scientific << std::sqrt(norm2(dJ - dJold) / norm2(dJold))
                << std::defaultfloat << std::endl;
      {
        std::vector<double> oms, lms;
        timeInterleaved(UGrid, nwarm, nmeas,
                        [&] { smart.logDetJacobianForce(dJ); },
                        [&] { smart.logDetJacobianForce(0, dJold); }, oms, lms);
        sJo = stats(oms); sJl = stats(lms);
        report(kname, "jacforce", "opt", sJo);
        report(kname, "jacforce", "old", sJl);
      }

      // ---- transformed gauge force (single-path: no old overload of
      //      AnalyticSmearedForce); raw deriv on the smeared field timed
      //      separately as context (not FTHMC cost)
      LatticeGaugeField Sigma0(UGrid), Sigma(UGrid);
      sR = stats(timeComponent(UGrid, nwarm, nmeas, [&] { GaugeAction.deriv(V, Sigma0); }));
      report(kname, "rawderiv", "single", sR);

      std::vector<double> xms;
      for (int i = 0; i < nwarm + nmeas; i++) {
        Sigma = Sigma0;                    // smeared_force overwrites its argument
        double t0 = now_sync(UGrid);
        smart.smeared_force(Sigma);
        double t1 = now_sync(UGrid);
        if (i >= nwarm) xms.push_back((t1 - t0) / 1.0e3);
      }
      sX = stats(xms);
      report(kname, "xformforce", "single", sX);
    } // smart destructed before the next kernel's instance

    for (auto s : stouts) delete s;

    // ---- per-kernel summary table (median ms per call)
    char line[160];
    std::cout << GridLogMessage << "FTHMCBENCH-TABLE kernel " << kname
              << " : median ms/call over " << nmeas << " calls" << std::endl;
    snprintf(line, sizeof(line), "%-11s %12s %12s %9s", "component", "non-opt[ms]", "opt[ms]", "speedup");
    std::cout << GridLogMessage << "FTHMCBENCH-TABLE " << line << std::endl;
    auto row1 = [&](const char *c, const Stats &s, const char *note) {
      snprintf(line, sizeof(line), "%-11s %12s %12.3f %9s  %s", c, "-", s.med, "-", note);
      std::cout << GridLogMessage << "FTHMCBENCH-TABLE " << line << std::endl;
    };
    auto row2 = [&](const char *c, const Stats &so, const Stats &sl) {
      snprintf(line, sizeof(line), "%-11s %12.3f %12.3f %8.2fx", c, sl.med, so.med,
               (so.med > 0) ? sl.med / so.med : 0.0);
      std::cout << GridLogMessage << "FTHMCBENCH-TABLE " << line << std::endl;
    };
    row1("forward", sF, "");
    { std::stringstream n; n << "(fp_iters_total=" << itot << ")";
      row1("inverse", sI, n.str().c_str()); }
    row2("lndet",    sLo, sLl);
    row2("jacforce", sJo, sJl);
    row1("xformforce", sX, "");
    row1("rawderiv", sR, "(gauge deriv on smeared field; not FTHMC cost)");
  }

  std::cout << GridLogMessage << "FTHMCBENCH done" << std::endl;
  Grid_finalize();
  return 0;
}
