/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/hmc/Test_WilsonFlow.cc

Copyright (C) 2017

Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Shuhei Yamamoto <syamamoto@bnl.gov>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
#include <string>
#include "ACC.hpp"


int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);

  typedef typename PeriodicGimplR::ComplexField ComplexField;

  typedef Grid::XmlReader Serialiser;
  Serialiser              Reader("input_ACF.xml", false, "root");
  WFParameters            WFPar(Reader);
  ConfParameters          CPar(Reader);
  ACFParameters           APar(Reader);

  const std::string file_path = WFPar.path;
  const int tau = WFPar.tau;
  const int W   = APar.MScut;  // Madras-Sokal cutoff: G[t>W] unused in MS error estimate

  const int total_configs = CPar.EndConfiguration - CPar.StartConfiguration + 1;
  const int T             = total_configs / APar.MDtime_div_fac;  // configs per bin
  const int n_bin         = APar.MDtime_div_fac;                  // number of MD-time bins
  const int arr_size      = T * n_bin;                            // total stored ACF entries

  if (APar.MDtime_div_fac == 1) assert(APar.R > 0 || W >= 100);

  int debug = GridCmdOptionExists(argv, argv+argc, "--debug");

  ComplexField one(&Grid); one = ComplexField::scalar_type(1.0, 0.0);
  ComplexField A0(&Grid), A1(&Grid), avg(&Grid);

  // F[i]:   raw observable field for configuration i
  // G[it]:  NVS autocovariance — subtracts ensemble mean <<a>> (a uniform scalar per bin)
  // G2[it]: NVS autocovariance — subtracts local mean field avg(x) (spatially varying)
  std::vector<ComplexField> F(total_configs, &Grid);
  std::vector<ComplexField> G(arr_size, &Grid), G2(arr_size, &Grid);

  std::cout << std::setprecision(15);

  ////////////////////////
  // Retrieve the fields
  ////////////////////////
  for (int i = 0; i < total_configs; i++) {
    int conf = CPar.StartConfiguration + i;
    std::string fname = file_path + WFPar.data_name + "_" + std::to_string(tau)
                      + "_" + CPar.conf_prefix + "." + std::to_string(conf);
    readFile(F[i], fname);
    if (debug) {
      RealD out = real(sum(F[i]));
      std::cout << GridLogMessage << "NVS " << WFPar.data_name
                << " (conf, tau, val): " << conf << " " << tau << " "
                << out / RealD(Grid.gSites()) << std::endl;
    }
  }

  /////////////////////////////////////////
  // Compute NVS autocovariance function
  /////////////////////////////////////////
  // NVS: use the ensemble mean <<a>> as the subtracted mean, rather than the per-config
  // spatial mean \bar a(x).  Two estimators are computed:
  //   G[it](x)  = avg over i of  (A(x,i) - <<a>>)(A(x,i+t) - <<a>>)   using scalar <<a>>
  //   G2[it](x) = avg over i of  (A(x,i) - avg(x))(A(x,i+t) - avg(x)) using field  avg(x)
  for (int i_bin = 0; i_bin < n_bin; i_bin++) {

    // Estimate the ensemble mean within this bin
    avg = Zero();
    for (int t = 0; t < T; t++)
      avg = avg + (1.0 / RealD(T)) * F[i_bin * T + t];

    // Scalar ensemble mean: spatial average of avg(x)
    ComplexD avg_scalar = TensorRemove(sum(avg)) / RealD(Grid.gSites());

    for (int t = 0; t < T; t++) {
      int it = i_bin * T + t;
      G[it] = Zero();  G2[it] = Zero();

      // n_src: number of source times used (full time average or single source)
      int n_src = APar.isFullTimeAvg ? T - t : 1;
      for (int i = 0; i < n_src; i++) {
        A0 = F[i_bin * T + i];
        A1 = F[i_bin * T + i + t];

        G[it]  = G[it]  + (A0 - avg_scalar * one) * (A1 - avg_scalar * one);
        G2[it] = G2[it] + (A0 - avg)               * (A1 - avg);
      }

      G[it]  = (1.0 / RealD(n_src)) * G[it];
      G2[it] = (1.0 / RealD(n_src)) * G2[it];
    }
  }

  ///////////////////////////
  // Data for Error Estimate
  ///////////////////////////
  for (const int& bs : APar.space_block_sizes) {

    // Build coarsened lattice with block size bs in all directions
    Coordinate clatt_size(Nd);
    for (int i = 0; i < Nd; i++) clatt_size[i] = Grid.FullDimensions()[i] / bs;
    GridCartesian Coarse(clatt_size, simd_layout, mpi_layout);

    ////////////////////// Block averaging //////////////////////
    std::vector<ComplexField> G_B(arr_size, &Coarse), G2_B(arr_size, &Coarse);
    for (int i = 0; i < arr_size; i++) {
      blockSum(G_B[i],  G[i]);   G_B[i]  = (1.0 / RealD(std::pow(bs, Nd))) * G_B[i];
      blockSum(G2_B[i], G2[i]);  G2_B[i] = (1.0 / RealD(std::pow(bs, Nd))) * G2_B[i];
    }

    ////////////////////// Sparse sampling //////////////////////
    // Retain only sites where every coordinate is a multiple of bs
    std::vector<ComplexField> G_s(arr_size, &Coarse), G2_s(arr_size, &Coarse);
    {
      LatticeInteger coor(&Grid);
      ComplexField filter(&Grid), zero(&Grid); filter = one; zero = Zero();
      for (int d = 0; d < Nd; d++) {
        LatticeCoordinate(coor, d);
        filter = where(mod(coor, bs) == Integer(0), filter, zero);
      }
      ComplexField tmp(&Grid);
      for (int i = 0; i < arr_size; i++) {
        tmp = filter * G[i];   blockSum(G_s[i],  tmp);
        tmp = filter * G2[i];  blockSum(G2_s[i], tmp);
      }
    }

    /***********   Error estimation  *************************************
      No binning  => Madras-Sokal approximation (valid when t << T)
      With binning => sample variance over MD-time bins
    *********************************************************************/
    std::string tag  = WFPar.data_name + " ACC";
    std::string tag2 = WFPar.data_name + " ACC2";

    if (total_configs == T) {
      // Madras-Sokal approximation
      assert(APar.R > 0 || T >= W);
      MS_approx(&Coarse, G_B,  T, W, bs, tau, "Blocked " + tag,  "NVS");
      MS_approx(&Coarse, G_s,  T, W, bs, tau, "Sparsed " + tag,  "NVS");
      MS_approx(&Coarse, G2_B, T, W, bs, tau, "Blocked " + tag2, "NVS");
      MS_approx(&Coarse, G2_s, T, W, bs, tau, "Sparsed " + tag2, "NVS");
    } else {
      // Binning over MD-time bins
      binning( &Coarse, G_B,  T, n_bin, bs, tau, "Blocked " + tag,  "NVS");
      binning( &Coarse, G_s,  T, n_bin, bs, tau, "Sparsed " + tag,  "NVS");
      binning2(&Coarse, G_B,  T, n_bin, bs, tau, "Blocked " + tag,  "NVS");
      binning2(&Coarse, G_s,  T, n_bin, bs, tau, "Sparsed " + tag,  "NVS");

      binning( &Coarse, G2_B, T, n_bin, bs, tau, "Blocked " + tag2, "NVS");
      binning( &Coarse, G2_s, T, n_bin, bs, tau, "Sparsed " + tag2, "NVS");
      binning2(&Coarse, G2_B, T, n_bin, bs, tau, "Blocked " + tag2, "NVS");
      binning2(&Coarse, G2_s, T, n_bin, bs, tau, "Sparsed " + tag2, "NVS");
    }
  }

  Grid_finalize();
}  // main


/*
Input file example


JSON

{
    "WilsonFlow":{
	"steps": 200,
	"step_size": 0.01,
	"meas_interval": 50,
  "maxTau": 2.0
    },
    "Configurations":{
	"conf_prefix": "ckpoint_lat",
	"rng_prefix": "ckpoint_rng",
	"StartConfiguration": 3000,
	"EndConfiguration": 3000,
	"Skip": 5
    }
}


*/
