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

// User record embedded in each evec_density Scidac file.
// Stores the H_DWF eigenvalue and mode index for self-documentation.
// tau (Wilson-flow time) is already encoded in the file name; not repeated here.
namespace Grid {
  struct H_DWF_EvalRecord : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(H_DWF_EvalRecord,
      double, eval,   // eigenvalue of H_DWF = gamma5*R5*D_DWF(mass), i.e. eMe[i]
      int,    n       // mode index sorted by |mu_n| ascending (0 = lowest)
    );
  };
}

template <class T, class RecordT>
void writeFile(T& in, std::string const fname, RecordT& record){
#ifdef HAVE_LIME
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0);
  WR.close();
#endif
  // What is the appropriate way to throw error?
}

template <class T> void writeFile(T& in, std::string const fname){
  Grid::emptyUserRecord record;
  writeFile(in, fname, record);
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
				    bool, take_meas,
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
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size, simd_layout, mpi_layout);
  
  XmlReader  Reader("LanParams.xml");
  LanczosParameters LanParams;   
  read(Reader,"LanczosParameters",LanParams);
#if 0 //DEBUG
  {
    std::cout << GridLogMessage << LanParams <<std::endl;
    XmlWriter HMCwr("LanParams.xml.out");
    write(HMCwr,"LanczosParameters",LanParams);
  }
#endif
  WFParameters WFParams(Reader);

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

  /*************  Finds eigenvectors of D_H = \gamma_5 R_5 D_dwf  ******************/
  std::vector<std::vector<FermionField>> conv_evecs_all;
  for(int i_conf=LanParams.StartTrajectory; i_conf<LanParams.StartTrajectory + LanParams.Trajectories; i_conf++) {

    FieldMetaData header;
    std::string file(LanParams.fpath + "/" + LanParams.fname + "." + std::to_string(i_conf));
    NerscIO::readConfiguration(Umu,header,file);

    /**********    Compute E density and TC density during WF     ****************/
    if(WFParams.take_meas){
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
    }
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
    std::cout << GridLogMessage << "Sorted G5R5M Evals: " << eMe       << std::endl;
    std::cout << GridLogMessage << "Sorted <G5R5M(evec), G5R5M(evec)>" << std::endl;
    std::cout << GridLogMessage << eMMe                                << std::endl;

    // Write eigenvalue text file: one value per line, sorted by |mu_n| ascending.
    // Read back in FieldDensityEigen via --evals for q_Bp (m_gap/mu_n) weighting.
    if( UGrid->IsBoss() ){
      std::string eval_file = LanParams.outpath + "/" + std::to_string(i_conf) +
                              "/eigenvalues_tau_" + tau + "." + std::to_string(i_conf);
      FILE *fp_eval = fopen(eval_file.c_str(), "w");
      assert(fp_eval != NULL);
      for(int i = 0; i < Nconv; i++)
        fprintf(fp_eval, "%.17g\n", eMe[i]);
      fclose(fp_eval);
      std::cout << GridLogMessage << "Wrote eigenvalues to: " << eval_file << std::endl;
    }

    conv_evecs_all.push_back(finalevec);    

    
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
    /***********************************************************************/
    /*   Four topological charge density estimators (q_A, q_B, q_B', q_C) */
    /*                                                                      */
    /*  q_B : sign(mu_n) weight, all 5D slices   — Bulk (Formula B)        */
    /*  q_Bp: m_gap/mu_n weight, all 5D slices   — Bulk improved (B')      */
    /*  q_C : sign(mu_n) weight, boundary slices — Boundary (Formula C)    */
    /*  q_A : sign(mu_n) weight, midpoint slices — Midpoint (Formula A)    */
    /*                                                                      */
    /*  m_gap = min_n |mu_n|  ~  m_f + m_res  (spectral gap of H_DWF).    */
    /*  For q_Bp: bulk modes are suppressed by m_gap/|mu_n| << 1,          */
    /*  making q_Bp the recommended estimator for pointwise comparison      */
    /*  with gradient-flowed gauge q(x).  q_A/B/C are recommended for      */
    /*  global Q and topological susceptibility chi_t.                      */
    /***********************************************************************/

    // --- Spectral-gap estimate: m_gap ~ m_f + m_res ---
    // The smallest |mu_n| among converged eigenvectors approximates m_gap.
    RealD m_gap = std::fabs(eMe[0]);
    for(int i = 1; i < Nconv; i++)
      if(std::fabs(eMe[i]) < m_gap) m_gap = std::fabs(eMe[i]);
    std::cout << GridLogMessage << "m_gap estimate (min|mu_n|) = " << m_gap << std::endl;

    // --- 4D accumulators ---
    LatticeComplexD q_B_4D(UGrid), q_Bp_4D(UGrid), q_C_4D(UGrid), q_A_4D(UGrid);
    q_B_4D = Zero(); q_Bp_4D = Zero(); q_C_4D = Zero(); q_A_4D = Zero();

    for(int i = 0; i < Nconv; i++){
      RealD mu_n    = eMe[i];
      RealD sign_mu = (mu_n >= 0.0) ? 1.0 : -1.0;
      // m_gap/mu_n weight for Formula B'.  Guard against mu_n==0 (should not occur).
      RealD w_Bp    = (mu_n != 0.0) ? (m_gap / mu_n) : 0.0;

      // 5D scalar density: rho_n(x,s) = |psi_n(x,s)|^2  (already stored as G5evec)
      // G5evec[i] was built with eps_code(s) applied:
      //   s <  Ls/2: G5evec[i] = -finalevec[i]   (eps_code = -1)
      //   s >= Ls/2: G5evec[i] = +finalevec[i]   (eps_code = +1)
      // So localInnerProduct(finalevec[i], G5evec[i]) = eps_code(s)*rho_n(x,s).
      // Summing over s gives chi_n^B(x) = sum_s eps_code(s) * rho_n(x,s).
      LatticeComplexD chi_B(UGrid); chi_B = Zero();
      {
        LatticeComplexD g5rho5D = localInnerProduct(finalevec[i], G5evec[i]); // 5D field
        LatticeComplexD sl(UGrid);
        for(int s = 0; s < Ls; s++){
          ExtractSlice(sl, g5rho5D, s, 0);
          chi_B = chi_B + sl;   // accumulate sum_s eps_code(s)*rho_n(x,s)
        }
      }
      // Formula B:  q_B(x)  += sign(mu_n) * chi_B(x)
      q_B_4D  = q_B_4D  + sign_mu * chi_B;
      // Formula B': q_B'(x) += (m_gap/mu_n) * chi_B(x)   [proper bulk suppression]
      q_Bp_4D = q_Bp_4D + w_Bp    * chi_B;

      // 5D scalar density rho_n(x,s) = |psi_n(x,s)|^2 (no eps_code factor)
      LatticeComplexD rho5D = localInnerProduct(finalevec[i], finalevec[i]);

      // Formula C: q_C(x) += -sign(mu_n) * [rho_n(x,Ls-1) - rho_n(x,0)]
      {
        LatticeComplexD bdy_s0(UGrid), bdy_sLs(UGrid);
        ExtractSlice(bdy_s0,  rho5D, 0,    0);  // left wall  s=0
        ExtractSlice(bdy_sLs, rho5D, Ls-1, 0);  // right wall s=Ls-1
        q_C_4D = q_C_4D - sign_mu * (bdy_sLs - bdy_s0);
      }

      // Formula A: q_A(x) += -sign(mu_n) * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
      if(Ls >= 2){
        LatticeComplexD mid_lo(UGrid), mid_hi(UGrid);
        ExtractSlice(mid_lo, rho5D, Ls/2-1, 0);  // below midpoint
        ExtractSlice(mid_hi, rho5D, Ls/2,   0);  // above midpoint
        q_A_4D = q_A_4D - sign_mu * 0.5 * (mid_hi - mid_lo);
      }
    }

    // --- Report global charges ---
    std::cout << GridLogMessage
              << "TCD estimators:"
              << "  Q_B="  << real(TensorRemove(sum(q_B_4D)))
              << "  Q_B'=" << real(TensorRemove(sum(q_Bp_4D)))
              << "  Q_C="  << real(TensorRemove(sum(q_C_4D)))
              << "  Q_A="  << real(TensorRemove(sum(q_A_4D))) << std::endl;

    // --- Write 4D fields ---
    std::string obase = LanParams.outpath + "/" + std::to_string(i_conf) + "/";
    writeFile(q_B_4D,  obase + "topo_q_B_tau_"  + tau + "." + std::to_string(i_conf));
    writeFile(q_Bp_4D, obase + "topo_q_Bp_tau_" + tau + "." + std::to_string(i_conf));
    writeFile(q_C_4D,  obase + "topo_q_C_tau_"  + tau + "." + std::to_string(i_conf));
    writeFile(q_A_4D,  obase + "topo_q_A_tau_"  + tau + "." + std::to_string(i_conf));
    /******************* end four-estimator block ****************************/

    for(int i = 0; i < Nconv; i++){
      chiral_matrix_real[i].resize(Nconv);
      chiral_matrix[i].resize(Nconv);

      auto evdensity = localInnerProduct(finalevec[i],finalevec[i] );
      Grid::H_DWF_EvalRecord eval_rec;
      eval_rec.eval = eMe[i];
      eval_rec.n    = i;
      writeFile(evdensity,
		LanParams.outpath + "/" + std::to_string(i_conf) + "/evec_density" +
		"_"+std::to_string(i)+"_tau_"+tau+"."+std::to_string(i_conf),
		eval_rec);

      auto diag_g5density = localInnerProduct(finalevec[i],G5evec[i] );
      writeFile(diag_g5density,
                LanParams.outpath + "/" + std::to_string(i_conf) + "/g5_density" +
                "_"+std::to_string(i)+"_tau_"+tau+"."+std::to_string(i_conf));
      std::cout << GridLogMessage << "evec G5 evec: " << i << " " << TensorRemove(sum(diag_g5density)) << std::endl;
      
      for(int j = 0; j < Nconv; j++){
	chiral_matrix[i][j]      = innerProduct(finalevec[i], G5evec[j]);
	chiral_matrix_real[i][j] = abs(chiral_matrix[i][j]);

	std::cout << GridLogMessage << " chiral_matrix_cplx "<<i<<" "<<j<<" "<< real(chiral_matrix[i][j]) <<" "<< imag(chiral_matrix[i][j]) << std::endl;
	std::cout << GridLogMessage << " chiral_matrix_real "<<i<<" "<<j<<" "<< chiral_matrix_real[i][j] << std::endl;
	if ( chiral_matrix_real[i][j] > 0.8 ) {
	  auto g5density = localInnerProduct(finalevec[i], G5evec[j]);
	  writeFile(g5density,
		    LanParams.outpath + "/" + std::to_string(i_conf) + "/chiral_density_" +
		    std::to_string(i)+"_"+std::to_string(j)+"_tau_"+tau+"."+std::to_string(i_conf));
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
  int Nconf = LanParams.Trajectories;
  RealD tensor[(Nconf-1)*Nconf/2][Nk][Nk]; for(int i=0; i<(Nconf-1)*Nconf/2; i++)for(int ii=0; ii<Nk; ii++)for(int jj=0; jj<Nk; jj++) tensor[i][ii][jj] = 0.0;
  int counter = 0;
  for(int i=0; i<conv_evecs_all.size()-1; i++){
    for(int j=i+1; j<conv_evecs_all.size(); j++)
      for(int ii=0; ii<conv_evecs_all[i].size(); ii++)
        for(int jj=0; jj<conv_evecs_all[j].size(); jj++){
          tensor[counter+j][ii][jj] = abs(innerProduct(conv_evecs_all[i][ii],conv_evecs_all[j][jj]));
	}
    counter += Nconf - 1 - i;
  }
  
  FILE *fp;
  if( UGrid->IsBoss()) fp = fopen((LanParams.outpath + "/evec_tensor_tau_"+tau).c_str(),"w");
  for(int i=0; i<conv_evecs_all.size()-1; i++)
    for(int j=i+1; j<conv_evecs_all.size(); j++)
      for(int ii=0; ii<Nk; ii++){
	for(int jj=0; jj<Nk; jj++)
	  if( UGrid->IsBoss()) fprintf(fp,"%lf ", tensor[i*(Nconf-1)-i*(i-1)/2+j][ii][jj]);
	if( UGrid->IsBoss()) fprintf(fp,"\n");
      }
  if( UGrid->IsBoss()) fclose(fp);
  
  Grid_finalize();
}
