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

using namespace Grid;

typedef DomainWallFermionD FermionOp;
typedef typename DomainWallFermionD::FermionField FermionField;

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
				    RealD, M5 , 
				    Integer, Ls,
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
      mass              = 0.01;
      M5                = 1.8;
      Ls                = 16;
      Nk                = 20;
      Nstop             = Nk;
      Np                = 80;
      StartTrajectory   = 0;
      Trajectories      = 1;
    }

    void print_parameters() const {
      std::cout << GridLogMessage << "[Config parameters] Trajectories            : " << Trajectories << "\n";
      std::cout << GridLogMessage << "[Config parameters] Start trajectory        : " << StartTrajectory << "\n";
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
  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  
  XmlReader  Reader("SFParams.xml");
  LanczosParameters LanParams;   
  read(Reader,"LanczosParameters",LanParams);
  WFParameters WFParams(Reader);
  SFParameters SFParams(Reader);

  double mass  = LanParams.mass;
  int    Ls    = LanParams.Ls;
  RealD  M5    = LanParams.M5;
  int    Nk    = LanParams.Nk;
  int    Nstop = LanParams.Nstop;
  int    Np    = LanParams.Np;
  
  GridCartesian *         UGrid   = &Grid;
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);     RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);     RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5rb(FrbGrid); RNG5.SeedFixedIntegers(seeds5);
  
  LatticeGaugeField Umu(UGrid), Uflow(UGrid);
  WilsonFlow<PeriodicGimplR> WF(WFParams.step_size,WFParams.steps,WFParams.meas_interval);
  int tmp = std::round(WFParams.maxTau);
  std::string tau = std::to_string(tmp);

  /*************  Finds eigenvectors of D_H = \gamma_5 R_5 D_dwf for each conf and variable m_f  ******************/
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
      FermionOp                                                  Ddwf(Uflow,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
      MdagMLinearOperator<FermionOp,FermionField>                HermOp(Ddwf);
      Gamma5R5HermitianLinearOperator<FermionOp, LatticeFermion> G5R5Herm(Ddwf);

      Chebyshev<FermionField>      Cheby   (LanParams.ChebyLow+abs(mass),LanParams.ChebyHigh+abs(mass),LanParams.ChebyOrder);
      FunctionHermOp<FermionField> OpCheby (Cheby,HermOp);
      PlainHermOp<FermionField>    Op      (HermOp);
      PlainHermOp<FermionField>    Op2     (G5R5Herm);
    
      ImplicitlyRestartedLanczos<FermionField> IRL(OpCheby, Op, Nstop, Nk, Nm, resid, MaxIt);
      
      /***********************************************************************/
      /*                    compute eigenvectors of D_H^2                    */
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
    
      /***********************************************************************/
      /*                       orthogonalization                             */
      /***********************************************************************/
      std::cout << GridLogMessage << "Start orthogonalization "     << std::endl;
      // calculat the matrix
      std::cout << GridLogMessage << "calculate the matrix element" << std::endl;
      std::vector<LatticeFermion> G5R5Mevec(Nconv, FGrid);
      std::vector<LatticeFermion> finalevec(Nconv, FGrid);
      std::vector<RealD> eMe(Nconv), eMMe(Nconv);
      for(int i = 0; i < Nconv; i++){
	std::cout << GridLogMessage << "calculate the matrix element["<<i<<"]" << std::endl;
	G5R5Herm.HermOpAndNorm(evec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
      }
      std::cout << GridLogMessage << "Re<evec, G5R5M(evec)>: "    << std::endl;
      std::cout << GridLogMessage << eMe                          << std::endl;
      std::cout << GridLogMessage << "<G5R5M(evec), G5R5M(evec)>" << std::endl;
      std::cout << GridLogMessage << eMMe                         << std::endl;
      std::vector<std::vector<ComplexD>> VevecG5R5Mevec(Nconv);
      Eigen::MatrixXcd evecG5R5Mevec = Eigen::MatrixXcd::Zero(Nconv, Nconv);
      for(int i = 0; i < Nconv; i++){
	VevecG5R5Mevec[i].resize(Nconv); //can be just RealD tmp; or static conversion?
	for(int j = 0; j < Nconv; j++){
	  VevecG5R5Mevec[i][j] = innerProduct(evec[i], G5R5Mevec[j]);
	  evecG5R5Mevec(i, j) = VevecG5R5Mevec[i][j];
	}
      }
      // calculate eigenvector
      std::cout << GridLogMessage << "Eigen solver" << std::endl;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(evecG5R5Mevec);
      std::vector<RealD> eigeneval(Nconv);
      std::vector<std::vector<ComplexD>> eigenevec(Nconv);
      for(int i = 0; i < Nconv; i++){
	eigeneval[i] = eigensolver.eigenvalues()[i];
	eigenevec[i].resize(Nconv);
	for(int j = 0; j < Nconv; j++){
	  eigenevec[i][j] = eigensolver.eigenvectors()(i, j);
	}
      }
      //rotation
      std::cout << GridLogMessage << "Do rotation" << std::endl;
      for(int i = 0; i < Nconv; i++){
	finalevec[i] = finalevec[i] - finalevec[i];
	for(int j = 0; j < Nconv; j++){
	finalevec[i] = eigenevec[j][i]*evec[j] + finalevec[i];
      }
      }
      //normalize again;
      for(int i = 0; i < Nconv; i++){
	RealD tmp_RealD = norm2(finalevec[i]);
	tmp_RealD = 1./pow(tmp_RealD, 0.5);
	finalevec[i] = finalevec[i]*tmp_RealD;
      }
      
      //check
      for(int i = 0; i < Nconv; i++){
	G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
      }
      
      /***********************************************************************/
      /*                    sort the eigenvectors                            */
      /***********************************************************************/
      std::vector<LatticeFermion> finalevec_copy(Nconv, FGrid);
      for(int i = 0; i < Nconv; i++){
	finalevec_copy[i] = finalevec[i];
      }
      std::vector<RealD> eMe_copy(eMe);
      for(int i = 0; i < Nconv; i++){
	eMe[i] = fabs(eMe[i]);
	eMe_copy[i] = eMe[i];
      }

      sort(eMe_copy.begin(), eMe_copy.end());
      for(int i = 0; i < Nconv; i++){
	for(int j = 0; j < Nconv; j++){
	  if(eMe[j] == eMe_copy[i]){
	    finalevec[i] = finalevec_copy[j];
	  }
	}
      }
      for(int i = 0; i < Nconv; i++){
	G5R5Herm.HermOpAndNorm(finalevec[i], G5R5Mevec[i], eMe[i], eMMe[i]);
      }
      std::cout << GridLogMessage << "Sorted Re<evec, G5R5M(evec)>: "               << std::endl;
      std::cout << GridLogMessage << "Sorted G5R5M Evals: m= "<< mass << " " << eMe << std::endl;
      std::cout << GridLogMessage << "Sorted <G5R5M(evec), G5R5M(evec)>"            << std::endl;
      std::cout << GridLogMessage << eMMe                                           << std::endl;
      //conv_evecs_all.push_back(finalevec);    

    
      /***********************************************************************/
      /*                   calculate chirality matrix                        */
      /***********************************************************************/
      std::vector<LatticeFermion>        G5evec(Nconv, FGrid);
      std::vector<std::vector<ComplexD>> chiral_matrix(Nconv);
      std::vector<std::vector<RealD>>    chiral_matrix_real(Nconv);
      for(int i = 0; i < Nconv; i++){
	G5evec[i] = Zero();//Did it take into account R5?
      for(int j = 0; j < Ls/2; j++){
	axpby_ssp(G5evec[i], -1., finalevec[i], 0., G5evec[i], j, j);
      }
      for(int j = Ls/2; j < Ls; j++){
	axpby_ssp(G5evec[i], 1., finalevec[i], 0., G5evec[i], j, j);
      }
      }
      // Compute spectral reconstruction of topological charge density (q_naive / sp_sum).
      // NOTE: this archived block has a sign-convention bug — G5evec[i] is built (above)
      // with eps_code(s) = -Gamma_5(s) (i.e. -1 for s<Ls/2, +1 for s>=Ls/2), so the
      // per-mode "-localInnerProduct(u, G5evec)" actually evaluates to +chi^B(x),
      // opposite to what the variable name suggests. The corrected reconstruction
      // (with cleanly-named accumulators) lives in
      // visualisation/FieldDensityEigen.cxx (q_naive accumulator); see
      // sp_sum_vs_qB_bulk_analysis.tex §1.3 for the derivation.
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
      for(int i = 0; i < Nconv; i++){
	chiral_matrix_real[i].resize(Nconv);
	chiral_matrix[i].resize(Nconv);

	auto evdensity = localInnerProduct(finalevec[i],finalevec[i] );
	writeFile(evdensity,
		  LanParams.outpath + "/" + std::to_string(i_conf) + "/evec_density" +
		  "_"+std::to_string(i)+"_tau_"+tau+"_m_"+std::to_string(mass)+"."+std::to_string(i_conf));
	
	for(int j = 0; j < Nconv; j++){
	  chiral_matrix[i][j] = innerProduct(finalevec[i], G5evec[j]);
	  chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);

	  std::cout << GridLogMessage << " chiral_matrix_cplx "<<i<<" "<<j<<" "<< real(chiral_matrix[i][j]) <<" "<< imag(chiral_matrix[i][j]) << std::endl;
	  std::cout << GridLogMessage << " chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
	  
	  if ( chiral_matrix_real[i][j] > 0.8 ) {
	    auto g5density = localInnerProduct(finalevec[i], G5evec[j]);
	    writeFile(g5density,
		      LanParams.outpath + "/" + std::to_string(i_conf) + "/chiral_density_" +
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
			  "/chiral_matrix_real_"+"tau_"+tau+"_m_"+std::to_string(mass)+"_"+std::to_string(i_conf)).c_str(),"w");
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
