/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Original file: ./tests/Test_dwf_G5R5.cc

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>
Author: Shuhei Yamamoto <syamamoto@bnl.gov>
From Duo and Bob's Chirality study

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


template <class T> void writeFile(T& in, std::string const fname){
#ifdef HAVE_LIME
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0);
  WR.close();
#endif
  // What is the appropriate way to throw error?
}

namespace Grid {

  struct LanczosParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParameters,
				    RealD, mass , 
				    Integer, Nstop,
				    Integer, Nk,
				    Integer, Np,
				    RealD, ChebyLow,
				    RealD, ChebyHigh,
				    Integer, ChebyOrder,
				    Integer, StartTrajectory,
				    Integer, Trajectories, /* @brief Number of configs processed in this run */
				    std::string, fpath,
				    std::string, fname,
				    std::string, outpath)
    
    LanczosParameters() {
      // Default values
      mass              = 0;
      Nk                = 20;
      Nstop             = Nk;
      Np                = 80;
      StartTrajectory   = 0;
      Trajectories      = 1;
    }

  };

  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
	    bool, is_flow,
            int, steps,
            double, step_size,
            int, meas_interval,
            double, maxTau, // for the adaptive algorithm
            std::string, path);


    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct SFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(SFParameters,
                                    double, massIncr,
                                    double, maxMass
                                    );


    template <class ReaderClass >
    SFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "SpectralFlow", *this);
    }

  };
  
}

int main(int argc, char** argv) {

  using namespace Grid;

  typedef WilsonFermionD FermionOp;
  typedef typename WilsonFermionD::FermionField FermionField;
  

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  
  XmlReader  Reader("SFWParams.xml");
  LanczosParameters LanParams;   
  read(Reader,"LanczosParameters",LanParams);
  WFParameters WFParams(Reader);
  SFParameters SFParams(Reader);
  
  double mass  = LanParams.mass;
  int    Nk    = LanParams.Nk;
  int    Nstop = LanParams.Nstop;
  int    Np    = LanParams.Np;

  GridCartesian *         UGrid   = &Grid;
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = UGrid;
  GridRedBlackCartesian * FrbGrid = UrbGrid;

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG4(UGrid);     RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5(FGrid);     RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG5rb(FrbGrid); RNG5.SeedFixedIntegers(seeds5);

  std::vector<Complex> boundary = {1,1,1,-1};
  FermionOp::ImplParams Params(boundary);
  
  LatticeGaugeField Umu(UGrid), Uflow(UGrid);
  WilsonFlow<PeriodicGimplR> WF(WFParams.step_size,WFParams.steps,WFParams.meas_interval);
  int tmp = std::round(WFParams.maxTau);
  std::string tau = std::to_string(tmp);

  /*************  Finds eigenvectors of D_H = \gamma_5 D_W for each conf and variable mass  ******************/
  std::vector<std::vector<FermionField>> conv_evecs_all;
  for(int i_conf=LanParams.StartTrajectory; i_conf<LanParams.StartTrajectory + LanParams.Trajectories; i_conf++) {

    FieldMetaData header;
    std::string file(LanParams.fpath + "/" + LanParams.fname + "." + std::to_string(i_conf));
    NerscIO::readConfiguration(Umu,header,file);

    if( WFParams.is_flow )
      WF.smear(Uflow, Umu);
    else
      Uflow = Umu;

    // TODO: add the following in the measurement for WF if to be repeated for diff flow times
    std::cout << GridLogMessage << "Start: " << file << std::endl;
    
    int   Nm      = Nk + Np;
    int   MaxIt   = 10000;
    RealD resid   = 1.0e-5;
    
    FermionField src(FGrid);
    gaussian(RNG5, src);

    while( mass > - SFParams.maxMass ){
      FermionOp                                                WilsonOperator(Umu,*FGrid,*FrbGrid,mass,Params);
      MdagMLinearOperator<FermionOp,FermionField>              HermOp(WilsonOperator);
      Gamma5HermitianLinearOperator<FermionOp, LatticeFermion> HermOp2(WilsonOperator);

      Chebyshev<FermionField>      Cheby   (LanParams.ChebyLow,LanParams.ChebyHigh,LanParams.ChebyOrder);
      FunctionHermOp<FermionField> OpCheby (Cheby,HermOp);
      PlainHermOp<FermionField>    Op2     (HermOp2);
    
      ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op2, Nstop, Nk, Nm, resid, MaxIt);
      
      /***********************************************************************/
      /*                    compute eigenvectors of D_H = G_5 D_W            */
      /***********************************************************************/
      int Nconv;
      std::vector<RealD> eval(Nm);
      std::vector<FermionField> evec(Nm, FGrid);
      IRL.calc(eval, evec, src, Nconv);
      
      std::cout << GridLogMessage << "mass eval: " << mass <<" : " << eval        << std::endl;
      std::cout << GridLogMessage << " #evecs "   << evec.size() << std::endl;
      std::cout << GridLogMessage << " Nconv  "   << Nconv       << std::endl;
      std::cout << GridLogMessage << " Nm     "   << Nm          << std::endl;
      if ( Nconv > evec.size() ) Nconv = evec.size();
    
      /*
      std::string sp_file = LanParams.outpath + "/" + std::to_string(i_conf) + "/sp_sum_tau_"+tau+"."+std::to_string(i_conf);
      LatticeComplexD sp_sum(FGrid); sp_sum = Zero();
      for(int i = 0; i < Nconv; i++) {
	RealD sign = (eMe[i]>=0)? 1.0 : -1.0;
	RealD abs_lambda = sqrt(eMe[i]*eMe[i] - mass*mass);
	sp_sum = sp_sum - localInnerProduct(finalevec[i],G5evec[i]) + 0.5*sign*abs_lambda*localInnerProduct(finalevec[i],finalevec[i]);
      }
      LatticeComplexD sp_sum4D(UGrid), tmp_F(UGrid); sp_sum4D = Zero();
      for(int i=0; i<Ls;i++){
	ExtractSlice(tmp_F,sp_sum,i,0);
	sp_sum4D = sp_sum4D + tmp_F;
      }
      writeFile(sp_sum4D,sp_file);
      */
      
      /***********************************************************************/
      /*                   calculate chirality matrix                        */
      /***********************************************************************/
      Gamma g5(Gamma::Algebra::Gamma5) ;
      ComplexD dot;
      std::vector<FermionField> G5evec(Nconv, FGrid);

      for (int i = 0; i < Nconv ; i++)
	G5evec[i] = g5*evec[i];

      std::vector<std::vector<ComplexD>> chiral_matrix(Nconv);
      std::vector<std::vector<RealD>>    chiral_matrix_real(Nconv);

      // Compute spectral reconstruction of topological charge density
      for(int i = 0; i < Nconv; i++){
	chiral_matrix_real[i].resize(Nconv);
	chiral_matrix[i].resize(Nconv);

	auto evdensity = localInnerProduct(evec[i],evec[i] );
	writeFile(evdensity,
		  LanParams.outpath + "/" + std::to_string(i_conf) + "/evec_density_Wilson_" +
		  std::to_string(i)+"_tau_"+tau+"_m_"+std::to_string(mass)+"."+std::to_string(i_conf));
	
	for(int j = 0; j < Nconv; j++){
	  chiral_matrix[i][j] = innerProduct(evec[i], G5evec[j]);
	  chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);

	  std::cout << GridLogMessage << " chiral_matrix_cplx "<<i<<" "<<j<<" "<< real(chiral_matrix[i][j]) <<" "<< imag(chiral_matrix[i][j]) << std::endl;
	  std::cout << GridLogMessage << " chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
	  
	  if ( chiral_matrix_real[i][j] > 0.8 ) {
	    auto g5density = localInnerProduct(evec[i], G5evec[j]);
	    writeFile(g5density,
		      LanParams.outpath + "/" + std::to_string(i_conf) + "/chiral_density_Wilson_" +
		      std::to_string(i)+"_"+std::to_string(j)+"_tau_"+tau+"_m_"+std::to_string(mass)+"."+std::to_string(i_conf));
	  }
	}
      }
      for(int i = 0; i < Nconv; i++){
	if(chiral_matrix[i][i].real() < 0.){
	  chiral_matrix_real[i][i] = -1. * chiral_matrix_real[i][i];
	}
      }
    
      // Save Chiral matrix for the config as a text file
      if( UGrid->IsBoss()){
	FILE *fp = fopen((LanParams.outpath + "/" + std::to_string(i_conf) +
			  "/chiral_matrix_real_Wilson_"+"tau_"+tau+"_m_"+std::to_string(mass)+"_"+std::to_string(i_conf)).c_str(),"w");
	assert(fp!=NULL);
	for(int i = 0; i < Nconv; i++){
	  for(int j = 0; j < Nconv; j++){
	    fprintf(fp,"%lf ",chiral_matrix_real[i][j]);
	}
	  fprintf(fp,"\n");
	}
	fclose(fp);
      }
      src  = evec[0]+evec[1]+evec[2];
      mass += -SFParams.massIncr;
    }
  }
  
  Grid_finalize();
}
