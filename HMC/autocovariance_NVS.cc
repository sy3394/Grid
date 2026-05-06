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

  // F[i]:        raw observable field for configuration i
  // G_cent[it]:  NVS autocovariance — CENTERED form using scalar bin-mean <<a>>
  //                G_cent[it](x) = avg_i (A(x,i) - <<a>>)(A(x,i+t) - <<a>>)
  // G_conn[it]:  NVS autocovariance — variant using the spatially-varying field
  //                avg(x) as the subtracted mean
  //                G_conn[it](x) = avg_i (A(x,i) - avg(x))(A(x,i+t) - avg(x))
  // (Naming aligned with the VS driver: cent vs conn distinguishes which
  //  "mean" is subtracted.  See autocova_usage.md.)
  std::vector<ComplexField> F(total_configs, &Grid);
  std::vector<ComplexField> G_cent(arr_size, &Grid), G_conn(arr_size, &Grid);

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
  //   G_cent[it](x) = avg over i of  (A(x,i) - <<a>>)(A(x,i+t) - <<a>>)   scalar mean
  //   G_conn[it](x) = avg over i of  (A(x,i) - avg(x))(A(x,i+t) - avg(x)) field mean
  for (int i_bin = 0; i_bin < n_bin; i_bin++) {

    // Estimate the ensemble mean within this bin
    avg = Zero();
    for (int t = 0; t < T; t++)
      avg = avg + (1.0 / RealD(T)) * F[i_bin * T + t];

    // Scalar ensemble mean: spatial average of avg(x)
    ComplexD avg_scalar = TensorRemove(sum(avg)) / RealD(Grid.gSites());

    for (int t = 0; t < T; t++) {
      int it = i_bin * T + t;
      G_cent[it] = Zero();  G_conn[it] = Zero();

      // n_src: number of source times used (full time average or single source)
      int n_src = APar.isFullTimeAvg ? T - t : 1;
      for (int i = 0; i < n_src; i++) {
        A0 = F[i_bin * T + i];
        A1 = F[i_bin * T + i + t];

        G_cent[it] = G_cent[it] + (A0 - avg_scalar * one) * (A1 - avg_scalar * one);
        G_conn[it] = G_conn[it] + (A0 - avg)               * (A1 - avg);
      }

      G_cent[it] = (1.0 / RealD(n_src)) * G_cent[it];
      G_conn[it] = (1.0 / RealD(n_src)) * G_conn[it];
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

    ////////////////////// Block averaging (spatial-lattice coarsening) //////////////////////
    std::vector<ComplexField> G_cent_B(arr_size, &Coarse), G_conn_B(arr_size, &Coarse);
    for (int i = 0; i < arr_size; i++) {
      blockSum(G_cent_B[i], G_cent[i]);  G_cent_B[i] = (1.0 / RealD(std::pow(bs, Nd))) * G_cent_B[i];
      blockSum(G_conn_B[i], G_conn[i]);  G_conn_B[i] = (1.0 / RealD(std::pow(bs, Nd))) * G_conn_B[i];
    }

    ////////////////////// Sparse sampling (spatial-lattice coarsening) //////////////////////
    // Retain only sites where every coordinate is a multiple of bs.
    // (NOTE: "binning" is reserved for MD-chain partition — see ACC.hpp header.)
    std::vector<ComplexField> G_cent_s(arr_size, &Coarse), G_conn_s(arr_size, &Coarse);
    {
      LatticeInteger coor(&Grid);
      ComplexField filter(&Grid), zero(&Grid); filter = one; zero = Zero();
      for (int d = 0; d < Nd; d++) {
        LatticeCoordinate(coor, d);
        filter = where(mod(coor, bs) == Integer(0), filter, zero);
      }
      ComplexField tmp(&Grid);
      for (int i = 0; i < arr_size; i++) {
        tmp = filter * G_cent[i];  blockSum(G_cent_s[i], tmp);
        tmp = filter * G_conn[i];  blockSum(G_conn_s[i], tmp);
      }
    }

    /***********   Error estimation  *************************************
      No MD-time binning  =>  Madras-Sokal approximation (placeholder only;
                              see autocova_usage.md re: structural circularity)
      With MD-time binning =>  sample variance over MD-time bins
                               (binning_avg_rho preferred; binning_avg_cov as cross-check)
    *********************************************************************/
    std::string tag_cent = WFPar.data_name + " G_cent";
    std::string tag_conn = WFPar.data_name + " G_conn";

    if (total_configs == T) {
      // Madras-Sokal approximation — currently retained as a placeholder;
      // do not treat as a trustworthy error bar.
      assert(APar.R > 0 || T >= W);
      MS_approx(&Coarse, G_cent_B, T, W, bs, tau, "Blocked " + tag_cent, "NVS");
      MS_approx(&Coarse, G_cent_s, T, W, bs, tau, "Sparsed " + tag_cent, "NVS");
      MS_approx(&Coarse, G_conn_B, T, W, bs, tau, "Blocked " + tag_conn, "NVS");
      MS_approx(&Coarse, G_conn_s, T, W, bs, tau, "Sparsed " + tag_conn, "NVS");
    } else {
      // MD-time binning.  binning_avg_rho is the preferred reported estimator.
      // For binning, spatial blocked vs sparse is essentially a non-choice —
      // the bin-to-bin variance is what gives the error.
      binning_avg_cov(&Coarse, G_cent_B, T, n_bin, bs, tau, "Blocked " + tag_cent, "NVS");
      binning_avg_cov(&Coarse, G_cent_s, T, n_bin, bs, tau, "Sparsed " + tag_cent, "NVS");
      binning_avg_rho(&Coarse, G_cent_B, T, n_bin, bs, tau, "Blocked " + tag_cent, "NVS");
      binning_avg_rho(&Coarse, G_cent_s, T, n_bin, bs, tau, "Sparsed " + tag_cent, "NVS");

      binning_avg_cov(&Coarse, G_conn_B, T, n_bin, bs, tau, "Blocked " + tag_conn, "NVS");
      binning_avg_cov(&Coarse, G_conn_s, T, n_bin, bs, tau, "Sparsed " + tag_conn, "NVS");
      binning_avg_rho(&Coarse, G_conn_B, T, n_bin, bs, tau, "Blocked " + tag_conn, "NVS");
      binning_avg_rho(&Coarse, G_conn_s, T, n_bin, bs, tau, "Sparsed " + tag_conn, "NVS");
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
