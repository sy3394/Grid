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

namespace Grid{
  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
	    int, tau,
	    std::string, data_name,
	    std::string, path);
       

    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
           std::string, conf_prefix,
	   int, StartConfiguration,
	   int, EndConfiguration);
  
    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };

  struct ACFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ACFParameters,
	   int, MDtime_div_fac,
	   int, MScut,	// Cutoff time for Madras-Sokal approx. to error in ACC
	   int, space_block_size);

    template <class ReaderClass >
    ACFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Autocorrelations", *this);
    }

  };
}

template <class T> void readFile(T& out, std::string const fname){
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
}

template <class L, typename A> void MS_approx(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int W, int block_size, int tau, std::string const e_name){
  using namespace Grid;
  // applicable only if blocking the src fields
  // Otherwise, for block averaged G_B(t), its error is hard to estimate as we need to consider correlation of G_x(t) at diff. sites wihting a block

  RealD avgG[in.size()], ACC;
  for (int t=0; t<T; t++) {
    avgG[t] = TensorRemove(sum(in[t])).real();
    std::cout << GridLogMessage << "NVS MFACCC " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << t << " "
              << avgG[t]/avgG[0] << std::endl;
  }
  
  L var(Coarse), tmp(Coarse), Gt(Coarse), G0(Coarse);
  for (int t=0; 2*t<T-W; t++){ // \sigma_{\rho}(t) uses \rho(t') up to 2t+W
    ACC = avgG[t]/avgG[0];
    var = in[0] - in[t];
    var = var*var;
    for (int k=1; k<W+t; k++) {
      tmp = (1.0/avgG[t]/avgG[0]+2.0/avgG[0]/avgG[0])*in[k]*in[k] + in[std::abs(k-t)]*( (1.0/avgG[t]/avgG[t])*in[k+t] - (4.0/avgG[0]/avgG[t])*in[k]);
      var = var + (2.0/RealD(T))*tmp;
    }
    std::cout << GridLogMessage << "NVS Variance " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << t << " "
	      << real(sum(var))*ACC*ACC << std::endl;
  }
}

template <class L, typename A> void binning(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int n_bin, int block_size, int tau, std::string const e_name){
  using namespace Grid;

  L avg(Coarse), avg0(Coarse), var(Coarse), var0(Coarse), cov(Coarse), tmp(Coarse);
  for (int t=0; t<T; t++){
    
    avg = Zero(); var = Zero(); cov = Zero();
    
    // Find avg field over binns
    for (int b=0; b<n_bin; b++)
      avg = avg + (1/RealD(n_bin))*in[b*T + t];
    if ( t== 0 ) avg0 = avg;
    RealD ACC = TensorRemove(sum(avg)).real()/TensorRemove(sum(avg0)).real();
    std::cout << GridLogMessage << "LVS MFACC " + e_name + " (Binning): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << ACC << std::endl;

    // Find variance of Gt/G0
    for	(int b=0; b<n_bin; b++){
      tmp = in[b*T + t] - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
      cov = cov + (1/RealD((n_bin-1)*n_bin))*tmp*(in[b*T] - avg0);
    }
    RealD G0 = real(sum(avg0)), Gt = real(sum(avg)); // the factor of 1/V is cancel out by 1/V^2 from Cov[G0,Gt] in the expression of var of Gt/G0
    if ( t == 0 ) var0 = var;
    std::cout << GridLogMessage << "LVS Variance " + e_name + " (Binning): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0 - 2.0*real(sum(cov))/G0/Gt)*ACC*ACC 
	      << " " << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0)*ACC*ACC << " " << - 2.0*real(sum(cov))/G0/Gt*ACC*ACC << std::endl;
  }
}

template <class L, typename A> void binning2(Grid::GridBase *Coarse, std::vector<L,A>  const& in, std::vector<L,A>  const& in2, int T, int n_bin, int block_size, int tau, std::string const e_name){
  using namespace Grid;
  
  RealD avg, var, tmp;
  std::vector<RealD> in_sum(in.size()), in2_sum(in.size());

  for (int i=0; i<in.size(); i++) {
    in_sum[i] = TensorRemove(sum(in[i])).real();
    in2_sum[i] = TensorRemove(sum(in2[i])).real();
  }
  
  // Find avg field over binns
  for (int t=0; t<T; t++){
    avg = 0.0; var = 0.0;
    for (int b=0; b<n_bin; b++)
      avg = avg + (1.0/RealD(n_bin))*in_sum[b*T + t]/sqrt(in2_sum[b*T]*in2_sum[b*T+t]);
    std::cout << GridLogMessage << "LVS MFACC " + e_name + " (Binning2): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << avg << std::endl;

    // Find variance of Gt/G0 
    for (int b=0; b<n_bin; b++){
      tmp = in_sum[b*T + t]/sqrt(in2_sum[b*T]*in2_sum[b*T+t]) - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
    }
    std::cout << GridLogMessage << "LVS Variance " + e_name + " (Binning2): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << var << " " << (1-avg*avg)/sqrt(RealD(Coarse->gSites())-3) << std::endl;
  }
}


int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);

  typedef typename PeriodicGimplR::ComplexField ComplexField;
  
  typedef Grid::XmlReader       Serialiser;
  Serialiser Reader("input_ACF.xml", false, "root");
  WFParameters WFPar(Reader);
  ConfParameters CPar(Reader);
  ACFParameters APar(Reader);

  std::string fname;
  std::string file_path  = WFPar.path;
  int W = APar.MScut, tau = WFPar.tau;
  ComplexD coeff0, coeff1;
  int total_configs = CPar.EndConfiguration - CPar.StartConfiguration + 1;
  int T = total_configs/APar.MDtime_div_fac; // #processed configs per bin
  int arr_size = T*APar.MDtime_div_fac;  
  //if ( APar.MDtime_div_fac == 1 ) assert(W >= 100);

  ComplexField A0(&Grid), A1(&Grid), one(&Grid); one=ComplexField::scalar_type(1.0,0.0);
  std::vector<ComplexField> G(total_configs,&Grid), G2(total_configs,&Grid), A(total_configs,&Grid);
  
  std::cout << std::setprecision(15);
  for (int conf_s=CPar.StartConfiguration, i_bin=0, it=0; conf_s+T-1<=CPar.EndConfiguration; conf_s+=T, i_bin++){
    for (int t=0; t<T; t++, it++){
	
      G[it] = Zero(); G2[it] = Zero();
	
      for (int i=0; i<1; i++){ // if we store G_x(t) in a vector, we can loop over i first and then over t
	  
	int conf_0 = conf_s + i;
	int conf_1 = conf_s + i + t;
	  
	fname = file_path + WFPar.data_name + "_" + std::to_string(tau) + "_" + CPar.conf_prefix + ".";
	readFile(A0, fname + std::to_string(conf_0));
	readFile(A1, fname + std::to_string(conf_1));
	if (i==0){ //debug
	  RealD out = real(sum(A0));
	  std::cout << GridLogMessage << "LVS " + WFPar.data_name + " (conf, tau, val):   " << " " << conf_0 << " " << tau << " "
		    <<  out/Real(Grid.gSites()) << std::endl;
	}
	
	coeff0 = TensorRemove(sum(A0))/RealD(Grid.gSites());
	coeff1 = TensorRemove(sum(A1))/RealD(Grid.gSites());
	G2[it] = G2[it] + A0*A1 - (coeff0*coeff1)*one;
	  
	A0 = A0 - coeff0*one;
	A1 = A1 - coeff1*one;
	G[it] = G[it] + A0*A1;
	A[it] = A1*A1;

	if (i==0){ //debug
	  Coordinate scoor; for (int mu=0; mu < Nd; mu++) scoor[mu] = 0;
	  RealD a0 = real(peekSite(A1,scoor));
          std::cout << GridLogMessage << "LVS (at origin) " + WFPar.data_name + " (conf, tau, val):   " << " " << conf_0 << " " << tau << " "
                    <<  a0 << std::endl;
        }
      }
      //removed!!!
    }// END(t): loop within a bin for MD time
  }// END(cong_s): loop over bins
	
  // Error Estimate
  int bs = APar.space_block_size;

  ///// Reduce to sub-lattice first
  Coordinate clatt_size(Nd);
  for(int i=0;i<Nd;i++) clatt_size[i] = Grid.FullDimensions()[i]/bs;
  GridCartesian Coarse(clatt_size,simd_layout,mpi_layout);
      
  // blocking
  std::vector<ComplexField> G_B(arr_size,&Coarse), G2_B(arr_size,&Coarse), A_B(arr_size,&Coarse);
  for (int i=0; i<arr_size; i++){
    blockSum(G_B[i],G[i]); blockSum(G2_B[i],G2[i]); blockSum(A_B[i],A[i]);
    G_B[i] = (1/RealD(std::pow(bs,Nd)))*G_B[i]; G2_B[i] = (1/RealD(std::pow(bs,Nd)))*G2_B[i]; A_B[i] = (1/RealD(std::pow(bs,Nd)))*A_B[i];
  }
      
  // sparse sampling
  std::vector<ComplexField> G_s(arr_size,&Coarse), G2_s(arr_size,&Coarse), A_s(arr_size,&Coarse);
  for (int i=0; i<arr_size; i++){
	
    LatticeInteger coor(&Grid);
    ComplexField tmp(&Grid), filter(&Grid), zero(&Grid); filter = one; zero = Zero();
    for (int d=0; d<Nd; d++) {
      LatticeCoordinate(coor,d);
      filter = where(mod(coor,bs)==Integer(0),filter,zero);
    }
    tmp = filter*G[i]; blockSum(G_s[i],tmp);
    tmp = filter*G2[i]; blockSum(G2_s[i],tmp);
    tmp = filter*A[i]; blockSum(A_s[i],tmp);
  }

  /*
    No Binning => Use Madras-Sokal approximation for error estimation
    Otherwise  => error estimate via sample variance by binning over MD time
  */
  int n_bin = total_configs/T;
  if ( total_configs == T ){
    // Madras-Sokal Approximation
    //   Valid: when t << T
    assert( T > W );
    MS_approx(&Coarse, G_B, T, W, bs, tau, "Blocked " + WFPar.data_name + " ACC");
    MS_approx(&Coarse, G2_B, T, W, bs, tau, "Blocked " + WFPar.data_name + " ACC2");
    
    MS_approx(&Coarse, G_s, T, W, bs, tau, "Sparsed " + WFPar.data_name + " ACC");
    MS_approx(&Coarse, G2_s, T, W, bs, tau, "Sparsed " + WFPar.data_name + " ACC2");
  }
  else {
    // Binning
    //binning(&Coarse, G_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC");
    //binning(&Coarse, G2_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC2");
    
    //binning(&Coarse, G_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC");
    //binning(&Coarse, G2_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC");

    binning2(&Coarse, G_B, A_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC");
    binning2(&Coarse, G2_B, A_B, T, n_bin, bs, tau, "Blocked " + WFPar.data_name + " ACC2");

    binning2(&Coarse, G_s, A_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC");
    binning2(&Coarse, G2_s, A_s, T, n_bin, bs, tau, "Sparsed " + WFPar.data_name + " ACC2");
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
