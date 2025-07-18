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
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  LanczosParameters LanParams;   
  XmlReader  LanReader("LanParams.xml");
  read(LanReader,"LanczosParameters",LanParams);
#if 0 //DEBUG
  {
    std::cout << GridLogMessage << LanParams <<std::endl;
    XmlWriter HMCwr("LanParams.xml.out");
    write(HMCwr,"LanczosParameters",LanParams);
  }
#endif
  WFParameters WFParams(LanReader);

  int   mass  = LanParams.mass;
  int   Ls    = LanParams.Ls;
  RealD M5    = LanParams.M5;
  int   Nk    = LanParams.Nk;
  int   Nstop = LanParams.Nstop;
  int   Np    = LanParams.Np;
  
  GridCartesian *         UGrid   = SpaceTimeGrid::makeFourDimGrid(
								   GridDefaultLatt(),
								   GridDefaultSimd(Nd, vComplex::Nsimd()),
								   GridDefaultMpi()
								   );
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

  /*************  Finds eigenvectors of D_H = \gamma_5 R_5 D_dwf  ******************/
  std::vector<std::vector<FermionField>> conv_evecs_all;
  for(int i_conf=LanParams.StartTrajectory; i_conf<LanParams.StartTrajectory + LanParams.Trajectories; i_conf++) {

    FieldMetaData header;
    std::string file(LanParams.fpath + "/" + LanParams.fname + "." + std::to_string(i_conf));
    NerscIO::readConfiguration(Umu,header,file);

    /**********    Compute E density and TC density during WF     ****************/
    std::string file_pre  = WFParams.path + "/";
    std::string file_post = LanParams.fname + "." + std::to_string(i_conf);
    WF.addMeasurement(WFParams.meas_interval, [&file_pre,&file_post,&i_conf](int step, RealD t, const typename PeriodicGimplR::GaugeField &U){

      int tmp = std::round(t);
      std::string tau = std::to_string(tmp);
  
      typedef typename PeriodicGimplR::GaugeLinkField GaugeMat;
      typedef typename PeriodicGimplR::ComplexField ComplexField;
      
      assert(Nd == 4);
      
      GaugeMat F(U.Grid());
      ComplexField R(U.Grid());
      R = Zero();
      
      for(int mu=0;mu<3;mu++){
	for(int nu=mu+1;nu<4;nu++){
	  WilsonLoops<PeriodicGimplR>::FieldStrength(F, U, mu, nu);
	  R = R + trace(F*F);
	}
      }
      R = (-1.0) * R;
      
      //// Taken from qcd/utils/WilsonLoops.h
      
      // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
      GaugeMat Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
      WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, U, Ydir, Zdir);
      WilsonLoops<PeriodicGimplR>::FieldStrength(By, U, Zdir, Xdir);
      WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, U, Xdir, Ydir);
      
      // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
      GaugeMat Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
      WilsonLoops<PeriodicGimplR>::FieldStrength(Ex, U, Tdir, Xdir);
      WilsonLoops<PeriodicGimplR>::FieldStrength(Ey, U, Tdir, Ydir);
      WilsonLoops<PeriodicGimplR>::FieldStrength(Ez, U, Tdir, Zdir);
      
      double coeff = 8.0/(32.0*M_PI*M_PI);
      ComplexField qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
      
      std::string efile = file_pre + "E_dnsty_" + tau + "_" + file_post;
      writeFile(R,efile);
      std::string tfile = file_pre + "Top_dnsty_" + tau + "_" + file_post;
      writeFile(qfield,tfile);
      
      RealD WFlow_TC5Li   = WilsonLoops<PeriodicGimplR>::TopologicalCharge5Li(U);
      RealD E = real(sum(R))/ RealD(U.Grid()->gSites());
      RealD T = real( sum(qfield) );
      Coordinate scoor; for (int mu=0; mu < Nd; mu++) scoor[mu] = 0;
      RealD E0 = real(peekSite(R,scoor));
      RealD T0 = real(peekSite(qfield,scoor));
      std::cout << GridLogMessage << "[WilsonFlow] Saved energy density (clover) & topo. charge density: "  << i_conf << " " << step << "  " << tau << "  "
		<< "(E_avg,T_sum) " << E << " " << T << " (E, T at origin) " << E0 << " " << T0 << " 5Li " << WFlow_TC5Li << std::endl;    
    });

    if( WFParams.is_flow )
      WF.smear(Uflow, Umu);
    else
      Uflow = Umu;

    // TODO: add the following in the measurement for WF if to be repeated for diff flow times
    std::cout << GridLogMessage << "Start: " << file << std::endl;
    
    int   Nm      = Nk + Np;
    int   MaxIt   = 10000;
    RealD resid   = 1.0e-5;

    FermionOp                                                  Ddwf(Uflow,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
    MdagMLinearOperator<FermionOp,FermionField>                HermOp(Ddwf);
    Gamma5R5HermitianLinearOperator<FermionOp, LatticeFermion> G5R5Herm(Ddwf);

    Chebyshev<FermionField>      Cheby   (LanParams.ChebyLow,LanParams.ChebyHigh,LanParams.ChebyOrder);
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
    FermionField src(FGrid);
    gaussian(RNG5, src);
    IRL.calc(eval, evec, src, Nconv);

    std::cout << GridLogMessage << mass <<" : " << eval        << std::endl;
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
    int pair_flag = 1;
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
    std::cout << GridLogMessage << "Sorted Re<evec, G5R5M(evec)>: "    << std::endl;
    std::cout << GridLogMessage << eMe                                 << std::endl;
    std::cout << GridLogMessage << "Sorted <G5R5M(evec), G5R5M(evec)>" << std::endl;
    std::cout << GridLogMessage << eMMe                                << std::endl;
    conv_evecs_all.push_back(finalevec);    

    
    /***********************************************************************/
    /*                   calculate chirality matrix                        */
    /***********************************************************************/
    std::vector<LatticeFermion>        G5evec(Nconv, FGrid);
    std::vector<std::vector<ComplexD>> chiral_matrix(Nconv);
    std::vector<std::vector<RealD>>    chiral_matrix_real(Nconv);
    for(int i = 0; i < Nconv; i++){
      G5evec[i] = Zero();
      for(int j = 0; j < Ls/2; j++){
	axpby_ssp(G5evec[i], -1., finalevec[i], 0., G5evec[i], j, j);
      }
      for(int j = Ls/2; j < Ls; j++){
	axpby_ssp(G5evec[i], 1., finalevec[i], 0., G5evec[i], j, j);
      }
    }
    // Compute spectral reconstruction of topological charge density
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

    for(int i = 0; i < Nconv; i++){
      chiral_matrix_real[i].resize(Nconv);
      chiral_matrix[i].resize(Nconv);
      
      std::string evfile(LanParams.outpath + "/" + std::to_string(i_conf) + "/evec_density");
      evfile = evfile+"_"+std::to_string(i)+"_tau_"+tau+"."+std::to_string(i_conf);
      auto evdensity = localInnerProduct(finalevec[i],finalevec[i] );
      writeFile(evdensity,evfile);
      
      for(int j = 0; j < Nconv; j++){
	chiral_matrix[i][j] = innerProduct(finalevec[i], G5evec[j]);
	std::cout << GridLogMessage << " chiral_matrix_real signed "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
	chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);
	std::cout << GridLogMessage << " chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
	if ( chiral_matrix_real[i][j] > 0.8 ) {
	  auto g5density = localInnerProduct(finalevec[i], G5evec[j]);
	  std::string chfile(LanParams.outpath + "/" + std::to_string(i_conf) + "/chiral_density_");
	  chfile = chfile +std::to_string(i)+"_"+std::to_string(j)+"_tau_"+tau+"."+std::to_string(i_conf);
	  writeFile(g5density,chfile);
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
      FILE *fp = fopen((LanParams.outpath + "/" + std::to_string(i_conf) + "/chiral_matrix_real_"+"tau_"+tau+"_"+std::to_string(i_conf)).c_str(),"w");
      assert(fp!=NULL);
      for(int i = 0; i < Nconv; i++){
	for(int j = 0; j < Nconv; j++){
	  fprintf(fp,"%lf ",chiral_matrix_real[i][j]);
	}
	fprintf(fp,"\n");
      }
      fclose(fp);
    }
  }
  // Compute tensor of <evecs_i(ii), evecs_j(jj)> where evecs_i is the conv'ed evecs for i^th config
  // row major
  FILE *fp;
  if( UGrid->IsBoss()) fp = fopen((LanParams.outpath + "/evec_tensor_tau_"+tau).c_str(),"w");
  for(int i=0; i<conv_evecs_all.size()-1; i++)
    for(int j=i+1; j<conv_evecs_all.size(); j++)
      for(int ii=0; ii<conv_evecs_all[i].size(); ii++){
	if( UGrid->IsBoss()) for(int jj=0; jj<ii; jj++) fprintf(fp,"%lf ", 0.0);
	for(int jj=ii; jj<conv_evecs_all[j].size(); jj++){
	  RealD tmp = abs(innerProduct(conv_evecs_all[i][ii],conv_evecs_all[j][jj]));
	  if( UGrid->IsBoss()) fprintf(fp,"%lf ", tmp);
	  }
	  if( UGrid->IsBoss()) fprintf(fp,"\n");
	}
  if( UGrid->IsBoss()) fclose(fp);
  
  Grid_finalize();
}
