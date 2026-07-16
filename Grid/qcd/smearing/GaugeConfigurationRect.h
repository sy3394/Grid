
/*!
  @file GaugeConfigurationRect.h
  @brief Declares the GaugeConfigurationRect class

  Conventions:
    - Luscher's normalization and inner product
    - Explicitly, T^a = -i*t^a where T^a:Luscher, t^a:Grid (c.f. SUnAdjoint.h)
*/
#pragma once

NAMESPACE_BEGIN(Grid);

/*!
  @brief Smeared configuration masked container for rectagle flow
  Modified for a multi-subset smearing (aka Luscher Flowed HMC)
*/
#undef DEBUG
//#define DEBUG

template <class Gimpl>
class SmearedConfigurationRect : public SmearedConfigurationMasked<Gimpl>
// TODO:
//can be a child of SmearedConfiguration but wanted to keep JacobianAction template form
// we can possibly add a virtual method logDetJacobian in the SmearedConfiguration class to make it work
// To use this class as it is, you need to include both this file and GaugeConfigurationMasked.h
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

private:
  
  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
  typedef typename SU3::vAlgebraMatrix vAlgebraMatrix;

  // These live in base class
  //  const unsigned int smearingLevels;
  //  Smear_Stout<Gimpl> *StoutSmearing;
  //  std::vector<GaugeField> SmearedSet;

  // Conventions:
  //   - flow kernel: 1 =:= Wilson,    2 =:= short-side rectangle
  //   -  mask types: 1 =:= red-black, 2 =:= red-black in 2x2 blocks of the dirs perp. to \mu,
  //                                       alternating site-by-site along \mu:
  //                                       mask = ( x_\mu + \sum_{\nu!=\mu} floor(x_\nu/2) ) mod 2

  int Nsmr_one_step = 2*Nd; // = #filterings(even colour, odd colour) x #dirs of smearing
  std::vector<int> mask_types;
  std::vector<Smear_Stout<Gimpl> *> Stouts;
  std::vector<LatticeLorentzComplex> masks; // should we turn this to poiners?????????

  // Optimised rect staple: ghost exchange of depth 2 (paths reach two hops in nu)
  // + per-mu stencils, set up once in the constructor.
  PaddedCell GhostRect;
  std::vector<GeneralLocalStencil> gStencils_rectsmear;
  // Force-level PlaqL/PlaqR stencils on the depth-2 padded grid.
  // Entry order == kernel read order in logDetJacobianForceLevel; keep in step.
  std::vector<GeneralLocalStencil> gStencils_plqforce;  //  6 per (mu,nu) pair, 3-4 entries
  std::vector<GeneralLocalStencil> gStencils_rectforce; // 10 per (mu,nu) pair, 6 entries
  std::vector<GeneralLocalStencil> gStencils_plqsmear;  //  1 per mu, 6 entries per nu!=mu (plq staple)

  void ApplyMask(GaugeField &U,int smr)
  {
    LatticeComplex tmp(U.Grid());
    GaugeLinkField Umu(U.Grid());
    for(int mu=0;mu<Nd;mu++){
      Umu=PeekIndex<LorentzIndex>(U,mu);
      tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
      Umu=Umu*tmp;
      PokeIndex<LorentzIndex>(U, Umu, mu);
    }
  }

  // flow_kernel = 1 : plq, 2 : Rs
  void BaseSmear(GaugeLinkField& Cmu, const GaugeField& U,int mu,RealD rho, int flow_kernel) {
    GridBase *grid = U.Grid();
    WilsonLoops<Gimpl> WL;



    Cmu = Zero();
    switch (flow_kernel) { 
    case 1:
      {
	WL.Staple(Cmu, U, mu);  //nb staple conventions of IroIro and Grid differ by a dagger
	break;
      }
    case 2:
      {
	// TODO: prepare optimized version
	WL.RectStapleUnoptimisedRs(Cmu, U, mu);
	break;
      }
    }

    Cmu = adj(rho * Cmu);
#ifdef DEBUG
    std::cout << GridLogMessage << "BaseSmear: " << mu<<" "<<rho<<" "<<flow_kernel<<" "<<norm2(Cmu) << std::endl;//DEBUG
#endif
  }

  // Optimised Rs staple: same result as BaseSmear with flow_kernel=2, but computed
  // on the padded grid with the pre-built stencil. gU = GhostRect.ExchangePeriodic(U).
  void BaseSmear_ghost_rect(GaugeLinkField& Cmu, const GaugeField& gU, int mu, RealD rho) {
    GRID_TRACE("BaseSmear_ghost_rect");
    assert((int)gStencils_rectsmear.size() == Nd);
    GridBase *ggrid = gU.Grid();
    GaugeLinkField gC(ggrid);
    Rect_Stout<Gimpl>::RectStaplePaddedRs(gC, gU, gStencils_rectsmear[mu], mu, rho);
    Cmu = GhostRect.Extract(gC);
  }

  // Optimised plq staple on the padded grid: Cmu = rho * adj(staple), the
  // flw_knl==1 result of BaseSmear. Kernel as BaseSmear_ghost in
  // GaugeConfigurationMasked.h, but on the full grid (no checkerboard).
  // gU = GhostRect.ExchangePeriodic(U).
  void BaseSmear_ghost_plq(GaugeLinkField& Cmu, const GaugeField& gU, int mu, RealD rho) {
    GRID_TRACE("BaseSmear_ghost_plq");
    GridBase *ggrid = gU.Grid();
    GaugeLinkField gtmp(ggrid);
    {
      autoView( gtmp_v , gtmp, AcceleratorWrite);
      autoView( gU_v , gU, AcceleratorRead);
      autoView( gStencil_v, gStencils_plqsmear[mu], AcceleratorRead);
      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
	  typedef decltype(coalescedRead(gtmp_v[0])) LinkMat;

	  LinkMat tmp = Zero();
	  for(int nu=0;nu<Nd;nu++){
	    int inc = 6*(nu - (mu<=nu));
	    if (nu != mu) {
	      GeneralStencilEntry const* e = gStencil_v.GetEntry(0+inc,ss);
	      auto U_nu_x = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(nu)();
	      e = gStencil_v.GetEntry(1+inc,ss);
	      auto U_mu_xpnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(mu)();
	      e = gStencil_v.GetEntry(2+inc,ss);
	      auto Udag_nu_xpmu = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();

	      tmp()() = tmp()() + U_nu_x * U_mu_xpnu * Udag_nu_xpmu;

	      e = gStencil_v.GetEntry(3+inc,ss);
	      auto Udag_nu_xmnu = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
	      e = gStencil_v.GetEntry(4+inc,ss);
	      auto U_mu_xmnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(mu)();
	      e = gStencil_v.GetEntry(5+inc,ss);
	      auto U_nu_xpmu_mnu = coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd)(nu)();

	      tmp()() = tmp()() + Udag_nu_xmnu * U_mu_xmnu * U_nu_xpmu_mnu;
	    }
	  }
	  coalescedWrite(gtmp_v[ss],rho*tmp);
	});
    }
    Cmu = GhostRect.Extract(gtmp);
  }

  void BaseSmearDerivativeP(GaugeField& SigmaTerm,
			    const GaugeField& iLambda,
			    const GaugeField& U,
			    int mmu, RealD rho)
  {
    // Reference
    // Morningstar, Peardon, Phys.Rev.D69,054501(2004)
    // Equation 75
    // Computing Sigma_mu, derivative of S[fat links] with respect to the thin links
    // Output SigmaTerm

    GridBase *grid = U.Grid();

    WilsonLoops<Gimpl> WL;
    GaugeLinkField staple(grid), u_tmp(grid);
    GaugeLinkField iLambda_mu(grid), iLambda_nu(grid);
    GaugeLinkField U_mu(grid), U_nu(grid);
    GaugeLinkField sh_field(grid), temp_Sigma(grid);
    Real rho_munu, rho_numu;

    rho_munu = rho;
    rho_numu = rho;
    for(int mu = 0; mu < Nd; ++mu){
      U_mu       = peekLorentz(      U, mu);
      iLambda_mu = peekLorentz(iLambda, mu);

      for(int nu = 0; nu < Nd; ++nu){
	if(nu==mu) continue;

	U_nu       = peekLorentz(      U, nu);

	// Nd(nd-1) = 12 staples normally.
	// We must compute 6 of these
	// in FTHMC case
	if ( (mu==mmu)||(nu==mmu) )
	  WL.StapleUpper(staple, U, mu, nu);
	
	if(nu==mmu) {
	  iLambda_nu = peekLorentz(iLambda, nu);

	  temp_Sigma = -rho_numu*staple*iLambda_nu;  //ok
	  //-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);

	  sh_field = Cshift(iLambda_nu, mu, 1);// general also for Gparity?

	  temp_Sigma = rho_numu*sh_field*staple; //ok
	  //r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	}

	if ( mu == mmu ) { 
	  sh_field = Cshift(iLambda_mu, nu, 1);

	  temp_Sigma = -rho_munu*staple*U_nu*sh_field*adj(U_nu); //ok
	  //-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	}

	//	staple = Zero();
	sh_field = Cshift(U_nu, mu, 1);

	temp_Sigma = Zero();

	if ( mu == mmu )
	  temp_Sigma = -rho_munu*adj(sh_field)*adj(U_mu)*iLambda_mu*U_nu;

	if ( nu == mmu ) {
	  temp_Sigma += rho_numu*adj(sh_field)*adj(U_mu)*iLambda_nu*U_nu;

	  u_tmp = adj(U_nu)*iLambda_nu;
	  sh_field = Cshift(u_tmp, mu, 1);
	  temp_Sigma += -rho_numu*sh_field*adj(U_mu)*U_nu;
	}
	
	sh_field = Cshift(temp_Sigma, nu, -1);
	Gimpl::AddLink(SigmaTerm, sh_field, mu);

      }
    }
  }

  void BaseSmearDerivativeRs(GaugeField& SigmaTerm,
                            const GaugeField& iLambda,
                            const GaugeField& U,
                            int mmu, RealD rho)
  {
    // mmu: dir in which U is updated
    
    GridBase *grid = U.Grid();

    GaugeLinkField staple(grid), u_tmp(grid);
    GaugeLinkField iLambda_mu(grid), iLambda_nu(grid);
    GaugeLinkField U_mu(grid), U_nu(grid);
    GaugeLinkField sh_field(grid), temp_Sigma(grid);
    Real rho_munu, rho_numu;

    rho_munu = rho;
    rho_numu = rho;
    for(int mu = 0; mu < Nd; ++mu){
      U_mu       = peekLorentz(      U, mu);
      iLambda_mu = peekLorentz(iLambda, mu);

      for(int nu = 0; nu < Nd; ++nu){
        if(nu==mu) continue;

        U_nu = peekLorentz(U, nu);
	
	if ( (mu==mmu)||(nu==mmu) )
	  Rect_Stout<Gimpl>::RectStapleUnoptimisedRsUpper(staple, U, mu, nu);
	
	if(nu==mmu) {
	  iLambda_nu = peekLorentz(iLambda, nu);

	  // 1st: -
          temp_Sigma = -rho_numu*staple*iLambda_nu;
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);

	  // 2nd: -
	  sh_field = adj(U_mu)*Cshift(temp_Sigma*U_mu, mu, -1);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);

	  // 3rd: +
	  sh_field = Cshift(iLambda_nu, mu, 1);
	  temp_Sigma = rho_numu*sh_field*adj(U_mu)*Cshift(staple*U_mu,mu,-1);
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);

	  // 4th: +
	  sh_field = Cshift(U_mu*temp_Sigma, mu,1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);

	  // 5th: +
	  temp_Sigma = rho_numu*adj(U_mu)*Cshift(adj(staple*U_nu)*adj(U_mu)*iLambda_nu*U_nu,nu,-1);
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);

	  // 6th: +
	  sh_field = adj(U_mu)*Cshift(temp_Sigma*U_mu,mu,-1);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);

	  // 7th: -
	  temp_Sigma = -rho_numu*Cshift(Cshift(adj(U_nu)*iLambda_nu*adj(adj(U_nu)*Cshift(adj(U_mu)*Cshift(staple*U_mu,mu,-1)*U_mu,mu,-1)),mu,1),nu,-1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm,temp_Sigma, mu);

	  // 8th: -
	  sh_field = Cshift(U_mu*temp_Sigma,mu,1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);
	}

	if ( mu == mmu ) {
	  // 9th: -
	  sh_field = Cshift(U_nu,mu,1);
	  temp_Sigma = -rho_munu*sh_field*Cshift(sh_field*Cshift(adj(U_mu)*iLambda_mu,nu,1)*adj(U_nu),nu,1)*adj(U_nu);
	  Gimpl::AddLink(SigmaTerm,temp_Sigma, mu);

	  // 10th: -
	  sh_field = adj(sh_field);
	  temp_Sigma = -rho_munu*Cshift(sh_field*Cshift(sh_field*adj(U_mu)*iLambda_mu*U_nu,nu,-1)*U_nu,nu,-1);
	  Gimpl::AddLink(SigmaTerm,temp_Sigma, mu);
	}
      }
    }
  }

  void BaseSmearDerivative(GaugeField& SigmaTerm,
                            const GaugeField& iLambda,
                            const GaugeField& U,
			   int mmu, RealD rho, int flow_kernel)
  {
    switch(flow_kernel){
    case 1:
      BaseSmearDerivativeP(SigmaTerm,iLambda,U,mmu,rho);
      break;
    case 2:
      BaseSmearDerivativeRs(SigmaTerm,iLambda,U,mmu,rho);
      break;
    }
  }

  // Adjoint vector to GaugeField force
  //tmp note: the output deviates from Luscher's convntion by -1 as it is
  //            to account for deviation by -1 from force calculation
  void InsertForce(GaugeField &Fdet,AdjVectorField &Fdet_nu,int nu)
  {
    Complex ci(0,1);
    GaugeLinkField Fdet_pol(Fdet.Grid());
    Fdet_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generator(e, te);
      auto tmp=peekColour(Fdet_nu,e);
      Fdet_pol=Fdet_pol + ci*tmp*te; // Fdet_pol + ci*tmp*te; // <- to be changed to this 
    }
    pokeLorentz(Fdet, Fdet_pol, nu);
  }

  // tmp comment: no extra factor
  // Old implementation, kept for consistency checks against the fused
  // default; to be deleted once confirmed. The int old argument only
  // selects this overload.
  void ComputeNxy(int old, const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GRID_TRACE("ComputeNxy_old");
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, tb);
      tb = 2.0 * ci * tb; // - ci * tb; in Lucher's convention but multiplied the missing factor from below
      Nx = Ta( adj(PlaqL)*tb * PlaqR );
      SU3::LieAlgebraProject(NxAd,Nx,b);
    }
  }

  // tmp comment: orig. extra factor of (-2)*(-2)/(-2) = -2 <- multiplied the result by 2 but still deviation from Luscher by -1
  // Old implementation, kept for consistency checks against the fused
  // default; to be deleted once confirmed. The int old argument only
  // selects this overload.
  void Compute_MpInvJx_dNxxdSy(int old, const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy_old");
    GaugeLinkField UtaU(PlaqL.Grid());
    GaugeLinkField D(PlaqL.Grid());
    AdjMatrixField Dbc(PlaqL.Grid());
    AdjMatrixField Dbc_opt(PlaqL.Grid());
    LatticeComplex tmp(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    RealD t=0;
    RealD tp=0;
    RealD tta=0;
    RealD tpk=0;
    t-=usecond();
    for(int a=0;a<Ngen;a++) {
      tta-=usecond();
      SU3::generator(a, ta);
      ta = ci * ta; //2.0 * ci * ta;
      UtaU= adj(PlaqL)*ta*PlaqR; // 6ms
      tta+=usecond();
      ////////////////////////////////////////////
      // Could add this entire C-loop to a projection routine
      // for performance. Could also pick checkerboard on UtaU
      // and set checkerboard on result for 2x perf
      ////////////////////////////////////////////
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	tc = (2.0)* ci * tc; //2.0*ci*tc;//changd so that the entire method is consistent with Luscher's convention
	tp-=usecond(); 
	D = Ta( tc *UtaU); // 2ms
#if 1
	SU3::LieAlgebraProject(Dbc_opt,D,c); // 5.5ms
#else // extra factor of -1/2 from Lucher's convention
	for(int b=0;b<Ngen;b++){
	  SU3::generator(b, tb);
	  tmp =-trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // Adjoint rep
	}
#endif
	tp+=usecond();
      }
      //      Dump(Dbc_opt,"Dbc_opt");
      //      Dump(Dbc,"Dbc");
      tpk-=usecond();
      tmp = trace(MpInvJx * Dbc_opt);
      PokeIndex<ColourIndex>(Fdet2,tmp,a);
      tpk+=usecond();
    }
    t+=usecond();
    std::cout << GridLogPerformance << " Compute_MpInvJx_dNxxdSy_old " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
	      << " ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms"<<std::endl;
  }

  // Default (fused) implementation, port of the plaquette-kernel
  // optimisation in GaugeConfigurationMasked.h. Same result as the old
  // overload above: note this class's convention ta = i t^a (the
  // masked/plaquette class uses 2i t^a); the c-loop convention T'^c = 2i t^c
  // is internal to the one-argument LieAlgebraProject.
  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, const AdjMatrixField &MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy");
    GridBase *grid = PlaqL.Grid();
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    RealD t=-usecond();

    autoView(Fdet2_v,Fdet2,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    autoView(MpInvJx_v,MpInvJx,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for2d(ss,grid->oSites(),a,Ngen,nsimd,{
        typedef decltype(coalescedRead(MpInvJx_v[0])) adj_mat;

        adj_mat Dbc;
        ColourMatrix ta;

	SU3::generator(a, ta);
	ta = ci * ta;
	auto UtaU = adj(PlaqL_v(ss))*ta*PlaqR_v(ss);
	SU3::LieAlgebraProject(Dbc,UtaU);

        coalescedWrite(Fdet2_v[ss]()()(a),traceProduct(MpInvJx_v(ss),Dbc)()()());
      });
    t+=usecond();
    std::cout << GridLogPerformance << " Compute_MpInvJx_dNxxdSy " << t/1e3 <<" ms"<<std::endl;
  }

  // Default (fused) implementation of ComputeNxy (same conventions as the
  // old overload above: tb = 2i t^b is internal to the two-argument
  // LieAlgebraProject).
  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GRID_TRACE("ComputeNxy");
    GridBase *grid = PlaqL.Grid();
    RealD t=-usecond();

    autoView(NxAd_v,NxAd,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
        typedef decltype(coalescedRead(NxAd_v[0]))  adj_mat;
        adj_mat NxAd_site;
	SU3::LieAlgebraProject(NxAd_site,PlaqL_v(ss),PlaqR_v(ss));
        coalescedWrite(NxAd_v[ss],NxAd_site);
      });
    t+=usecond();
    std::cout << GridLogPerformance << " ComputeNxy " << t/1e3 <<" ms"<<std::endl;
  }

  // The site-local real-part inverse used inside the fused force kernel
  // lives in Lattice_trace.h (Inverse_RealPartSite), next to LUdcmp/solve
  // and Inverse_RealPart, which is now implemented on top of it.

  void linkTracer(const std::vector<GaugeLinkField> &Umu, const GaugeLinkField &Umskd, const std::vector<int> dirs0, int ind, Real rho, GaugeLinkField &rect){
    // dir in dirs is 1+mu where mu=0,..3 to put sign on dir
    // ind: index of dirs for which masked gauge link field should be used.  The indexing starts from 1, not 0
    //    actually ind is always equal to size of dirs0 => can become just a flag
    GaugeLinkField tmp(Umu[0].Grid()), U(Umu[0].Grid());

    std::vector<int> dirs = dirs0; 
    std::reverse(dirs.begin(), dirs.end());
    
    for(int i=0; i<dirs.size(); i++){
      int mu  = dirs[i];
      int sgn = mu/std::abs(mu);
      mu = std::abs(mu) - 1; // \in {0,..,Nd-1}

      if(ind>0 && dirs.size()-i-1==ind) U = Umskd;
      else                U = Umu[mu];
      if(i == 0){
	if(sgn>0)
	  tmp = Gimpl::CovShiftIdentityForward(U,mu);
	else
	  tmp = Gimpl::CovShiftIdentityBackward(U,mu);
      }else{
	if(sgn>0)
	  tmp = Gimpl::CovShiftForward(U,mu,tmp);
	else
	  tmp = Gimpl::CovShiftBackward(U,mu,tmp);
      }
    }
    rect = rho*tmp;
  }
    
public:

  // Old implementation, kept for consistency checks against the optimised
  // default below; to be deleted once confirmed. The int old argument only
  // selects this overload.
  void logDetJacobianForceLevel(int old, const GaugeField &U, GaugeField &force ,int smr)
  {
    GridBase* grid = U.Grid();
    ColourMatrix   tb;
    ColourMatrix   tc;
    ColourMatrix   ta;
    GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(grid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(grid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Nxx(grid);  // Nxx fundamental space
    GaugeLinkField Utmp(grid);
    GaugeLinkField PlaqL(grid);
    GaugeLinkField PlaqR(grid);
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    ColourMatrix Ident;
    LatticeComplex  cplx(grid);
    
    AdjVectorField  dJdXe_nMpInv(grid); 
    AdjVectorField  dJdXe_nMpInv_y(grid); 
    AdjMatrixField  MpAd(grid);    // Mprime luchang's notes
    AdjMatrixField  MpAdInv(grid); // Mprime inverse
    AdjMatrixField  NxxAd(grid);    // Nxx in adjoint space
    AdjMatrixField  JxAd(grid);     
    AdjMatrixField  ZxAd(grid);
    AdjMatrixField  mZxAd(grid);
    AdjMatrixField  X(grid);
    Complex ci(0,1);

    RealD t0 = usecond();
    Ident = ComplexD(1.0);
    for(int d=0;d<Nd;d++){
      Umu[d] = peekLorentz(U, d);
    }
    int mu= (smr/2) %Nd;

    ////////////////////////////////////////////////////////////////////////////////
    // Mask the gauge field
    ////////////////////////////////////////////////////////////////////////////////
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);
    std::cout << GridLogMessage << "norm Utmp: "<<mu<<" "<<norm2(Utmp)<<std::endl;
    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    int smr_ind = smr/Nsmr_one_step;
    int flw_knl = mask_types[smr_ind];
    double rho;
    
    switch(flw_knl){
    case 1:
      rho=this->Stouts[smr_ind]->SmearRho[1];
      break;
    case 2:
      rho=((Rect_Stout<Gimpl> *) this->Stouts[smr_ind])->SmearRhoRs[1];
      break;
    }
    int idx=0;
    for(int mu=0;mu<4;mu++){
      for(int nu=0;nu<4;nu++){
	double rho1;
	switch(flw_knl){
	case 1:
	  rho1=this->Stouts[smr_ind]->SmearRho[idx];
	  break;
	case 2:
	  rho1=((Rect_Stout<Gimpl> *) this->Stouts[smr_ind])->SmearRhoRs[idx];
	  break;
	}

	if ( mu!=nu) assert(rho1==rho);
	else         assert(rho1==0.0);
	idx++;
      }}
    std::cout << GridLogMessage <<"JacobActionREct: "<<rho<<" "<<smr_ind<<" "<<flw_knl<<" "<<mu<<std::endl;
    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only and do it here
    RealD time;
    time=-usecond();
    BaseSmear(Cmu, U,mu,rho, flw_knl);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Umu[mu]));
    time+=usecond();
    std::cout << GridLogMessage << "Z took "<<time<< " us"<<std::endl;

    time=-usecond();
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's
      SU3::generator(b, tb);         // <- traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      cplx = 2.0*trace(ci*tb*Zx);    // Luscher's norm conv. (c.f. top comment) // 2.0*trace(ci*tb*Zx); orig
      ZxAd = ZxAd - cplx * TRb;      // orig: ZxAd + cplx * TRb; after negating TRb, i.e., TRb=-TRb;
    }
    time+=usecond();
    std::cout << GridLogMessage << "ZxAd took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
#if 1
    time=-usecond();
    X=1.0; 
    JxAd = X;
    mZxAd = (-1.0)*ZxAd; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZxAd;
      kpfac = kpfac /(k+1);
      JxAd = JxAd + X * kpfac;
    }
    time+=usecond();
    std::cout << GridLogMessage << "Jx took "<<time<< " us"<<std::endl;
#endif
    //////////////////////////////////////
    // dJ(x)/dxe
    //////////////////////////////////////
    time=-usecond();
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    std::vector<AdjMatrix> TRb_s; TRb_s.resize(8);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);

#if 1
    // tmp comment: an extra factor of -2 here (just compare with theoretical form) <- the factor of -1 removed
    // Norm: the remaning factor of 2 is accouted for in the end of this function
    for(int b=0;b<8;b++){
      SU3Adjoint::generator(b, TRb_s[b]);
      dJdX[b] = TRb_s[b];
    }
    aunit = ComplexD(1.0);
    // Could put into an accelerator_for
    X  = (-1.0)*ZxAd; 
    t2 = X;
    for (int j = 12; j > 1; --j) {
      t3  = t2*(1.0 / (j + 1))  + aunit;
      t2  = X * t3;
      for(int b=0;b<8;b++){
	dJdX[b]= TRb_s[b] * t3 + X * dJdX[b]*(1.0 / (j + 1));
      }
    }
    //  The above computation comes with an additinal factor of 2 (just compared with the theoretical form or below)
    for(int b=0;b<8;b++){
      dJdX[b] = -0.5*dJdX[b];
    }
#endif

#if 0  //testing combined Horner's rule
    std::vector<AdjMatrixField>  XB_dJdX;    XB_dJdX.resize(8,grid);
    //std::vector<AdjMatrix> TRb_s; TRb_s.resize(8);
    AdjMatrixField XB_tbXn(grid);
    AdjMatrixField XB_sumXtbX(grid);
    AdjMatrixField XB_t2(grid);
    AdjMatrixField XB_dt2(grid);
    AdjMatrixField XB_t3(grid);
    AdjMatrixField XB_dt3(grid);
    AdjMatrixField XB_JxAd(grid);
    //AdjMatrixField aunit(grid);
    
    // 1. Pre-loop Negation
    X = (-1.0) *ZxAd;
    
    aunit = ComplexD(1.0);
    XB_t3 = aunit; 
    XB_t2 = X;
    for(int b=0; b<8; b++) {
      // Note: Since Y = -X, the derivative w/r/t X includes a -1 factor
      // As TRb_s[b] is - adj_basis in Luscher's convention, -1 is already included
      XB_dJdX[b] = Zero(); //TRb_s[b]; 
    }
    
    // 2. The Horner Loop 
    for (int j = 12; j >= 1; --j) {
      for(int b=0; b<8; b++) {
        // dJdX[b] correctly accumulates the negative contributions
        XB_dJdX[b] = (TRb_s[b] * XB_t3 + X * XB_dJdX[b]) *(1.0 / (j + 1));
      }

      XB_t3 = XB_t2 * (1.0 / (j + 1)) + aunit;
      XB_t2 = X * XB_t3;
      
    }
    XB_JxAd = XB_t3;
    std::cout << GridLogMessage << "DEBUG: Horner's method JxAd"<<norm2(XB_JxAd-JxAd)<< " djdx ";
    for(int i =0;i<dJdX.size();i++) std::cout<<i<<" "<<norm2(XB_dJdX[i] + dJdX[i])<<" "<<norm2(XB_dJdX[i])<<" "<<norm2(dJdX[i])<<" ";//note: sign deviation of dJdX from Luscher's convention
    std::cout <<std::endl;
    for(int i =0;i<dJdX.size();i++){
      Dump(dJdX[i], "DEBUG HO (dJdX");
      Dump(XB_dJdX[i], "DEBUG HO (XB_dJdX");
    }

#if 0

    for(int b=0;b<8;b++){
      dJdX[b] = -XB_dJdX[b];
    }
#endif
#endif
    
    time+=usecond();
    std::cout << GridLogMessage << "dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(old,PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "ComputeNxy took "<<time<< " us"<<std::endl;
    
    ////////////////////////////
    // Mab
    ////////////////////////////
    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxxAd;

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    time=-usecond();
    MpAdInv = Inverse(MpAd); //Inverse_RealPart(MpAd); //Inverse(MpAd);
    time+=usecond();
    std::cout << GridLogMessage << "MpAdInv took "<<time<< " us"<<std::endl;
    
    RealD t3a = usecond();
    /////////////////////////////////////////////////////////////////
    // Nxx Mp^-1
    /////////////////////////////////////////////////////////////////
    AdjVectorField  FdetV(grid);
    AdjVectorField  Fdet1_nu(grid);
    AdjVectorField  Fdet2_nu(grid);
    AdjVectorField  Fdet2_mu(grid);
    AdjVectorField  Fdet1_mu(grid);

    AdjMatrixField nMpInv(grid);
    nMpInv= NxxAd *MpAdInv;

    AdjMatrixField MpInvJx(grid);
    AdjMatrixField MpInvJx_nu(grid);
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor

    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx,FdetV);
    Fdet2_mu=FdetV;
    Fdet1_mu=Zero();
    
    for(int e =0 ; e<8 ; e++){
      LatticeComplexD tr(grid);
      //      ColourMatrix te;
      //      SU3::generator(e, te);
      tr = trace(dJdX[e] * nMpInv);
      pokeColour(dJdXe_nMpInv,tr,e);
    }
    ///////////////////////////////
    // Mask it off
    ///////////////////////////////
    auto tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
    dJdXe_nMpInv = dJdXe_nMpInv*tmp;
    
    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1)
    //       Nxy_mu one site forward  in each nu direction (3)
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.

    AdjMatrixField Nxy(grid);

    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);
    GaugeLinkField Fdet_pol(grid); // one polarisation

    std::vector<int> dirs;
    
    RealD t4 = usecond();
    for(int nu=0;nu<Nd;nu++){
      ////////////////////  Conventions  ///////////////////
      // || or == <- smeared link
      // : or ..  <- link w/r/t/ which derivative is taken
      //
      // nu dir
      // ^
      // |   
      // mu dir
      // ->
      // (x,mu): updating link
      // (y,nu): current link for which force is computed
      /////////////////////////////////////////////////////
      if (nu!=mu) {
	switch(flw_knl){
	case 1:
	  {

	    ///////////////// +ve nu /////////////////
	    //     __
	    //    :  |
	    //    x==    // nu polarisation -- clockwise
	    // x = y
	    
	    time=-usecond();
	    PlaqL=Ident;
	    PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
						Gimpl::CovShiftForward(Umu[mu], mu,
								       Gimpl::CovShiftBackward(Umu[nu], nu,
											       Gimpl::CovShiftIdentityBackward(Utmp, mu))));
	    time+=usecond();
	    std::cout << GridLogMessage << "PlaqLR took "<<time<< " us"<<std::endl;
	    
#if 0 //DEBUG
	    GaugeLinkField PlaqR2(grid);
	    dirs = {(nu+1),mu+1,-(nu+1),-(mu+1)};
	    linkTracer(Umu, Utmp, dirs, 3, -rho, PlaqR2);
	    std::cout << GridLogMessage << "DEBUG: PlaqR linkT: "<<norm2(PlaqR2-PlaqR)<<std::endl;
#endif
	
	    time=-usecond();
	    dJdXe_nMpInv_y =   dJdXe_nMpInv;
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	    time+=usecond();
	    std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;
	    
	    time=-usecond();
	    PlaqR=(-1.0)*PlaqR;
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx,FdetV);
	    Fdet2_nu = FdetV;
	    time+=usecond();
	    std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	    
	    //     __
	    //    |  :
	    //    x==y    // nu polarisation -- anticlockwise
	    
	    PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
					       Gimpl::CovShiftBackward(Umu[mu], mu,
								       Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));
#if 0 //DEBUG
	    //GaugeLinkField PlaqR2(grid);
	    dirs = {(nu+1),-(mu+1),-(nu+1)};
	    linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR2);
	    std::cout << GridLogMessage << "DEBUG: PlaqR linkT: "<<norm2(PlaqR2-PlaqR)<<std::endl;
#endif
	    PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    

	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	    
	    ///////////////// -ve nu /////////////////
	    // x==
	    // :  |
	    // y__|          // nu polarisation -- clockwise
	    
	    PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftForward(Umu[nu], nu,
								       Gimpl::CovShiftIdentityBackward(Utmp, mu)));
	    
	    PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	    
	    // x==
	    // |  :
	    // |__y         // nu polarisation
	    
	    PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
						Gimpl::CovShiftIdentityBackward(Utmp, mu));
	    
	    PlaqR=Gimpl::CovShiftBackward(Umu[mu], mu,
					  Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);
	    
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    /////////////////////////////////////////////////////////////////////
	    // Set up the determinant force contribution in 3x3 algebra basis
	    /////////////////////////////////////////////////////////////////////
	    InsertForce(Fdet1,Fdet1_nu,nu);
	    InsertForce(Fdet2,Fdet2_nu,nu);
	    
	    //////////////////////////////////////////////////
	    // Parallel direction terms
	    //////////////////////////////////////////////////
	    
	    //    y..
	    //    |  |
	    //    x==   // mu polarisation
	    PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftBackward(Umu[nu], nu,
									Gimpl::CovShiftIdentityBackward(Utmp, mu)));
	    
	    PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);
	    
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	    
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;

	    // x==
	    // |  |
	    // y..          // mu polarisation
	    
	    PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftForward(Umu[nu], nu,
								       Gimpl::CovShiftIdentityBackward(Utmp, mu)));

	    PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	    
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,1);
	    
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;

	    break;
	  }
	case 2:
	  {
	    ///////////////// +ve nu /////////////////
	    //     ->
	    //    |  |
	    //    :  |
	    //    x==    
	    // x = y : Computes contr. from this type to the force for U_nu(y) and U_nu(y+nu) <- similar effect happens in all calc. below
	  
	    time=-usecond();
	    PlaqL=Ident;
	    
	    dirs = {nu+1,nu+1,mu+1,-(nu+1),-(nu+1),-(mu+1) };
	    linkTracer(Umu, Utmp, dirs, 5, -rho, PlaqR);
	    time+=usecond();
	    std::cout << GridLogMessage << "PlaqLR took "<<time<< " us Rect"<<norm2(PlaqR)<<" mu "<<mu<<" nu "<<nu<<std::endl;
	    
	    time=-usecond();
	    dJdXe_nMpInv_y =   dJdXe_nMpInv;
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	    time+=usecond();
	    std::cout << GridLogMessage << "ComputeNxy (occurs 10x) took "<<time<< " us"<<std::endl;
	    
	    time=-usecond();
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx,FdetV);
	    Fdet2_nu = FdetV;
	    time+=usecond();
	    std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 10x) took "<<time<< " us"<<std::endl;

	    //     <-
	    //    |  |
	    //    |  :
	    //    x==y    
	    // x = y - mu 

	    dirs = {nu+1,nu+1,-(mu+1),-(nu+1),-(nu+1)};
	    //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	    PlaqR = rho * Gimpl::CovShiftForward(Umu[nu],nu,
						 Gimpl::CovShiftForward(Umu[nu],nu,
									Gimpl::CovShiftBackward(Umu[mu],mu,
												Gimpl::CovShiftBackward(Umu[nu],nu,
															Gimpl::CovShiftIdentityBackward(Umu[nu],nu)))));
	    PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu); // Note: adj(PlaqL) is used in ComputeNx & Compute_MpInvJx_dNxxdSy
	  
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	
	
	    //    :->
	    //    y  |
	    //    |  |
	    //    x==    
	    // x = y - nu
	    
	    PlaqL = Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	    dirs = {nu+1,(mu+1),-(nu+1),-(nu+1),-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 4,-rho, PlaqR);
	    PlaqR = (-rho) * Gimpl::CovShiftForward(Umu[nu],nu,
						    Gimpl::CovShiftForward(Umu[mu],mu,
									   Gimpl::CovShiftBackward(Umu[nu],nu,
												   Gimpl::CovShiftBackward(Umu[nu],nu,
															   Gimpl::CovShiftIdentityBackward(Utmp,mu)))));
												 
	  
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	  
	    MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	
	    //     <-:
	    //    |  y
	    //    |  |
	    //    x==    
	    // x = y - mu - nu

	    dirs = {-(nu+1),-(mu+1)};
	    //linkTracer(Umu,Utmp, dirs, 1, 1.0,PlaqL);
	    PlaqL = Gimpl::CovShiftBackward(Umu[nu],nu,
					   Gimpl::CovShiftIdentityBackward(Utmp,mu));
	    dirs = {nu+1,-(mu+1),-(nu+1),-(nu+1)};
	    //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	    PlaqR = rho * Gimpl::CovShiftForward(Umu[nu],nu,
						 Gimpl::CovShiftBackward(Umu[mu],mu,
									 Gimpl::CovShiftBackward(Umu[nu],nu,
												 Gimpl::CovShiftIdentityBackward(Umu[nu],nu))));
					       


	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,-1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,-1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	    
	
	    ///////////////// -ve nu /////////////////
	    
	    //    x==
	    //    :  |
	    //    y
	    //    |  |
	    //     <-
	    // x = y + nu: Computes contr. from this type to the force for U_nu(y) and U_nu(y+nu) <- similar effect happens in all calc. below
	    
	    dirs = {-(nu+1),(mu+1),(nu+1),(nu+1),-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 4, rho, PlaqL);
	    PlaqL = rho* Gimpl::CovShiftBackward(Umu[nu],nu,
						 Gimpl::CovShiftForward(Umu[mu],mu,
									Gimpl::CovShiftForward(Umu[nu],nu,
											       Gimpl::CovShiftForward(Umu[nu],nu,
														      Gimpl::CovShiftIdentityBackward(Utmp,mu)))));
														      
	    PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	    
	
	    //    x==
	    //    |  :
	    //       y
	    //    |  |
	    //     ->
	    // x = y - mu + nu

	    dirs = {-(nu+1),-(mu+1),(nu+1),(nu+1)};
	    //linkTracer(Umu, Utmp, dirs, -1, -rho, PlaqL);
	    PlaqL = (-rho)*Gimpl::CovShiftBackward(Umu[nu],nu,
						   Gimpl::CovShiftBackward(Umu[mu],mu,
									   Gimpl::CovShiftForward(Umu[nu],nu,
												  Gimpl::CovShiftIdentityForward(Umu[nu], nu))));
	    dirs = {(nu+1),-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 1, 1.0,PlaqR);
	    PlaqR = Gimpl::CovShiftForward(Umu[nu],nu,
					   Gimpl::CovShiftIdentityBackward(Utmp,mu));
	    
	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,1);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,1);
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;
	    

	    //    x==
	    //    |  |
	    //    :  |
	    //    y->
	    // x = y + 2nu
	    
	    dirs = {nu+1,nu+1};
	    //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	    PlaqR = rho*Gimpl::CovShiftForward(Umu[nu],nu,
					       Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	    dirs = {mu+1,-(nu+1),(nu+1),(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 3, 1.0, PlaqL);
	    PlaqL = Gimpl::CovShiftForward(Umu[mu],mu,
					   Gimpl::CovShiftForward(Umu[nu],nu,
								  Gimpl::CovShiftForward(Umu[nu],nu,
											 Gimpl::CovShiftIdentityBackward(Utmp,mu))));
	    
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,2);
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,2);
	    Compute_MpInvJx_dNxxdSy(old,PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu + FdetV;
	  
	    //    x==
	    //    |  |
	    //    |  :
	    //     ->y
	    // x = y + 2nu - mu
	    
	    dirs = {nu+1,nu+1,-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 2, -rho, PlaqR);
	    PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[nu],nu,
						  Gimpl::CovShiftForward(Umu[nu],nu,
									 Gimpl::CovShiftIdentityBackward(Utmp,mu)));
	    dirs = {-(mu+1),nu+1,nu+1};
	    //linkTracer(Umu, Utmp, dirs, -1, 1.0, PlaqL);
	    PlaqL = Gimpl::CovShiftBackward(Umu[mu],mu,
					    Gimpl::CovShiftForward(Umu[nu],nu,
								   Gimpl::CovShiftIdentityForward(Umu[nu], nu)));
	    
	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,2);
	    ComputeNxy(old,PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,2);
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    /////////////////////////////////////////////////////////////////////
	    // Set up the determinant force contribution in 3x3 algebra basis
	    /////////////////////////////////////////////////////////////////////
	    InsertForce(Fdet1,Fdet1_nu,nu);
	    InsertForce(Fdet2,Fdet2_nu,nu);
	    
	  
	    //////////////////////////////////////////////////
	    // Parallel direction terms
	    //////////////////////////////////////////////////
	    
	    //    y..
	    //    |  |
	    //    |  |
	    //    x=<=
	    // x = y - 2nu : Computes contr. from this type to the force for U_nu(y) and U_nu(y+nu) <- similar effect happens in all calc. below
	  
	    PlaqL=Gimpl::CovShiftBackward(Umu[nu],nu,Gimpl::CovShiftIdentityBackward(Umu[nu], nu));
	    dirs = {(mu+1),-(nu+1),-(nu+1),-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 3, -rho,PlaqR);
	    PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[mu],mu,
						  Gimpl::CovShiftBackward(Umu[nu],nu,
									  Gimpl::CovShiftBackward(Umu[nu],nu,
												  Gimpl::CovShiftIdentityBackward(Utmp,mu))));
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-2);
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	  
	    MpInvJx_nu = Cshift(MpInvJx,nu,-2);
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;
	  
	
	    //    x<=
	    //    |  |
	    //    |  |
	    //    y..
	    // x = y + 2nu

	    PlaqL=Gimpl::CovShiftForward(Umu[nu],nu,Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	    dirs = {(mu+1),(nu+1),(nu+1),-(mu+1)};
	    //linkTracer(Umu, Utmp, dirs, 3, -rho,PlaqR);
	    PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[mu],mu,
						  Gimpl::CovShiftForward(Umu[nu],nu,
									 Gimpl::CovShiftForward(Umu[nu],nu,
											      Gimpl::CovShiftIdentityBackward(Utmp,mu))));

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,2);
	    ComputeNxy(old,PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	    
	    MpInvJx_nu = Cshift(MpInvJx,nu,2);
	    Compute_MpInvJx_dNxxdSy(old,PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;
	    break;
	  }
	default:
	  {
	    assert(1!=1 && " At present, the only valid choice of flow kernel is either 1 or 2");
	    break;
	  }
	}
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    // The overall extra factor of 2 (from dJdX and Compute_MpInvJx_dNxxdSy) is removed
    // each also came with an extra factor of -1 from Luscher's convention, which is kept for the moment
    //   This is because InertForce come with extra factor of -1, which can be adjasted
    // The overall minus sign in force is necessary: Trivializing Maps, the Wilson Flow and the HMC Algorithm (Lushcer) Eq. (6.1)
    force=-1.0*(Fdet1 + Fdet2); //-0.5*(Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t4-t3 dJdXe_nMpInv "<<t4-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }
  
  // Default (optimised) implementation. Same result and per-term structure
  // as the old overload above; differences are performance-only:
  //  - flw_knl==2 staple from the depth-2 padded cell (RectStaplePaddedRs)
  //  - ZxAd via make_adjoint_rep; JxAd Taylor series fused sitewise
  //  - dJdX Horner + trace(dJdX nMpInv) fused in one kernel (factor -0.5 as
  //    in the old overload; the masked/plaquette class uses -1.0)
  //  - fused ComputeNxy / Compute_MpInvJx_dNxxdSy
  //  - MpAd inversion via Inverse_RealPart (production choice, GPU-resident;
  //    adjoint rep is real). The old overload keeps the complex Inverse, so
  //    the consistency check covers this difference (tolerance, not bitwise).
  void logDetJacobianForceLevel(const GaugeField &U, GaugeField &force ,int smr)
  {
    GRID_TRACE("logDetJacobianForceLevel");
    GridBase* grid = U.Grid();
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(grid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(grid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Utmp(grid);
    GaugeLinkField PlaqL(grid);
    GaugeLinkField PlaqR(grid);
    const int Ngen = SU3Adjoint::Dimension;
    ColourMatrix Ident;

    AdjVectorField  dJdXe_nMpInv(grid);
    AdjVectorField  dJdXe_nMpInv_y(grid);
    AdjMatrixField  NxxAd(grid);   // Nxx in adjoint space
    AdjMatrixField  ZxAd(grid);
    // Jx, Mab, MpAdInv, nMpInv live only inside the fused kernel below
    Complex ci(0,1);

    RealD t0 = usecond();
    Ident = ComplexD(1.0);
    for(int d=0;d<Nd;d++){
      Umu[d] = peekLorentz(U, d);
    }
    int mu= (smr/2) %Nd;

    ////////////////////////////////////////////////////////////////////////////////
    // Mask the gauge field
    ////////////////////////////////////////////////////////////////////////////////
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu);

    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);

    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    int smr_ind = smr/Nsmr_one_step;
    int flw_knl = mask_types[smr_ind];
    double rho;

    switch(flw_knl){
    case 1:
      rho=this->Stouts[smr_ind]->SmearRho[1];
      break;
    case 2:
      rho=((Rect_Stout<Gimpl> *) this->Stouts[smr_ind])->SmearRhoRs[1];
      break;
    }

    //////////////////////////////////////////////////////////////////
    // Padded (ghost) gauge field, per-direction links and masked link:
    // used by the staple (flw_knl==2) and by every stencil-built
    // PlaqL/PlaqR term in the nu loop below.
    //////////////////////////////////////////////////////////////////
    GridBase *ggrid = GhostRect.grids[Nd-1];
    GaugeField gU(grid);
    std::vector<GaugeLinkField> gUmu(Nd,grid);
    GaugeLinkField gUtmp(grid);
    {
      GRID_TRACE("ExchangePeriodicRect");
      gU = GhostRect.ExchangePeriodic(U);
      for(int d=0; d<Nd; d++) gUmu[d] = peekLorentz(gU, d);
      gUtmp = GhostRect.ExchangePeriodic(Utmp);
    }

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    switch(flw_knl){
    case 1:
      BaseSmear_ghost_plq(Cmu, gU, mu, rho);
      break;
    case 2:
      BaseSmear_ghost_rect(Cmu, gU, mu, rho);
      break;
    }

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Umu[mu]));

    // Move Z to the adjoint rep: ZxAd = -sum_b 2 tr(i t^b Zx) TRb
    {GRID_TRACE("ZxAdOpt");
      SU3Adjoint::make_adjoint_rep(ZxAd, Zx);
    }

    /////////////////////////////////////////////////////////////////
    // NxxAd (needed before the fused J/Mab kernel below)
    /////////////////////////////////////////////////////////////////
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(PlaqL,PlaqR,NxxAd);

    RealD t3a = usecond();
    /////////////////////////////////////////////////////////////////
    // Nxx Mp^-1
    /////////////////////////////////////////////////////////////////
    AdjVectorField  FdetV(grid);
    AdjVectorField  Fdet1_nu(grid);
    AdjVectorField  Fdet2_nu(grid);
    AdjVectorField  Fdet2_mu(grid);
    AdjVectorField  Fdet1_mu(grid);

    AdjMatrixField MpInvJx(grid);
    AdjMatrixField MpInvJx_nu(grid);

    /////////////////////////////////////////////////////////////////
    // ONE fused kernel: a single Horner recursion yields BOTH dJdX_b and
    // J (= 1 + t2/2 after the loop — same truncation as the reference
    // Taylor series, cf. the XB scheme in the old overload); then
    // Mab = 1 - Jx Nxx; its real-part 8x8 LU inverse in-kernel
    // (Inverse_RealPartSite, Lattice_trace.h — the adjoint rep is real;
    // production choice as in GaugeConfigurationMasked.h, the old overload
    // keeps the complex Inverse); the dJdXe traces; and
    // MpInvJx = -MpAdInv Jx. Jx, Mab, MpAdInv, nMpInv are kernel-local —
    // no intermediate lattice fields. On CPU builds only the LU section
    // serialises SIMD lanes (data-dependent pivoting); the rest stays
    // vectorised. NB Horner-J equals Taylor-J only in exact arithmetic:
    // the FP summation order differs (last-bit), visible as ~1e-27 rel^2
    // in the old-vs-default force check.
    /////////////////////////////////////////////////////////////////
    {GRID_TRACE("J_Mab_Inv_dJdX_fusedOpt");
      autoView(dJdXe_nMpInv_v,dJdXe_nMpInv,AcceleratorWrite);
      autoView(MpInvJx_v,MpInvJx,AcceleratorWrite);
      autoView(ZxAd_v,ZxAd,AcceleratorRead);
      autoView(NxxAd_v,NxxAd,AcceleratorRead);
      const int nsimd = vAlgebraMatrix::Nsimd();
      accelerator_for(ss,grid->oSites(),nsimd,{
	  typedef decltype(coalescedRead(ZxAd_v[0]))         adj_mat;
	  typedef decltype(coalescedRead(dJdXe_nMpInv_v[0])) adj_vec;
	  adj_mat X, t3, t2, aunit, JxAd_site, MpAd_site, MpAdInv_site, nMpInv_site;
	  adj_vec dJdXe_nMpInv_site;
	  iVector<adj_mat,Ngen> iTas;
	  iVector<adj_mat,Ngen> dJdX_b;

	  // One Horner recursion yields BOTH dJdX_b and J (the XB scheme of
	  // the old overload): after the j-loop, J = 1 + t2/2 reproduces the
	  // reference Taylor truncation sum_{k=0..11} X^k/(k+1)! exactly.
	  for(int b=0;b<Ngen;b++){
	    SU3Adjoint::generator(b, iTas(b));
	    dJdX_b(b) = iTas(b);
	  }
	  aunit = ComplexD(1.0);
	  X  = (-1.0)*ZxAd_v(ss);
	  t2 = X;
	  for (int j = 12; j > 1; --j) {
	    t3  = t2*(1.0 / (j + 1))  + aunit;
	    t2  = X * t3;
	    for(int b=0;b<Ngen;b++){
	      dJdX_b(b)= iTas(b) * t3 + X * dJdX_b(b)*(1.0 / (j + 1));
	    }
	  }
	  JxAd_site = aunit + t2*0.5;

	  // Mab = 1 - Jx Nxx and its real-part inverse, in registers
	  MpAd_site = Complex(1.0,0.0);
	  MpAd_site = MpAd_site - JxAd_site * NxxAd_v(ss);
	  Inverse_RealPartSite(MpAdInv_site, MpAd_site);

	  nMpInv_site= NxxAd_v(ss) * MpAdInv_site;
	  // factor -0.5: this class's dJdX normalisation, c.f. reference above
	  for(int e=0;e<Ngen;e++){
	    dJdXe_nMpInv_site()()(e) = traceProduct((-0.5)*dJdX_b(e),nMpInv_site)()()();
	  }
	  coalescedWrite(dJdXe_nMpInv_v[ss],dJdXe_nMpInv_site);

	  // MpInvJx = -MpAdInv Jx (rho is on the plaq factor)
	  coalescedWrite(MpInvJx_v[ss],(-1.0)*(MpAdInv_site*JxAd_site));
	});
    }

    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
    Fdet2_mu=FdetV;
    Fdet1_mu=Zero();
    ///////////////////////////////
    // Mask it off
    ///////////////////////////////
    {
      auto tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
      dJdXe_nMpInv = dJdXe_nMpInv*tmp;
    }

    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1)
    //       Nxy_mu one site forward  in each nu direction (3)
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.

    AdjMatrixField Nxy(grid);

    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);

    ///////////////////////////////////////////////////////////////////
    // Padded workspaces + views for the stencil-built PlaqL/PlaqR terms
    ///////////////////////////////////////////////////////////////////
    GaugeLinkField gPlaqL(ggrid), gPlaqR(ggrid);
    autoView( gPlaqL_v , gPlaqL, AcceleratorWrite);
    autoView( gPlaqR_v , gPlaqR, AcceleratorWrite);
    autoView( gU_mu_v  , gUmu[mu], AcceleratorRead);
    autoView( gUtmp_v  , gUtmp,    AcceleratorRead);

    RealD t4 = usecond();
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	autoView( gU_nu_v , gUmu[nu], AcceleratorRead);
	const int p6  = (mu*(Nd-1) + (nu-(mu<=nu)))*6;  // plq  stencil base
	const int p10 = (mu*(Nd-1) + (nu-(mu<=nu)))*10; // rect stencil base
	switch(flw_knl){
	case 1:
	  {
	    ///////////////// +ve nu /////////////////
	    //     __
	    //    :  |
	    //    x==    // nu polarisation -- clockwise
	    // x = y

	    PlaqL=Ident;
	    // PlaqR = -rho U_nu(x) U_mu(x+nu) U_nu^d(x+mu) msk^d(x)
	    {
	      GRID_TRACE("PlaqP1");
	      autoView( gStencil_v, gStencils_plqforce[p6+0], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x       =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_mu_xpnu    =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_nu_xpmu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(3,ss);
		  auto Mdag_x       = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_nu_x * U_mu_xpnu * Udag_nu_xpmu * Mdag_x);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y =   dJdXe_nMpInv;
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;

	    PlaqR=(-1.0)*PlaqR;
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	    Fdet2_nu = FdetV;

	    //     __
	    //    |  :
	    //    x==y    // nu polarisation -- anticlockwise

	    // PlaqR = rho U_nu(x) U_mu^d(x-mu+nu) U_nu^d(x-mu) ; PlaqL = msk^d(x-mu)
	    {
	      GRID_TRACE("PlaqP2");
	      autoView( gStencil_v, gStencils_plqforce[p6+1], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x         =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto Udag_mu_xmmupnu= adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_nu_xmmu   = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (rho) * U_nu_x * Udag_mu_xmmupnu * Udag_nu_xmmu);
		  e = gStencil_v.GetEntry(3,ss);
		  auto Mdag_xmmu      = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], Mdag_xmmu);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	      PlaqL = GhostRect.Extract(gPlaqL);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    ///////////////// -ve nu /////////////////
	    // x==
	    // :  |
	    // y__|          // nu polarisation -- clockwise

	    // PlaqL = rho U_mu(x) U_nu(x+mu) msk^d(x+nu) ; PlaqR = U_nu(x)
	    {
	      GRID_TRACE("PlaqP3");
	      autoView( gStencil_v, gStencils_plqforce[p6+2], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_mu_x     =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpmu  =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Mdag_xpnu  = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], (rho) * U_mu_x * U_nu_xpmu * Mdag_xpnu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = Umu[nu];
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    // x==
	    // |  :
	    // |__y         // nu polarisation

	    // PlaqL = -rho U_nu(x) msk^d(x-mu+nu) ; PlaqR = U_mu^d(x-mu) U_nu(x-mu)
	    {
	      GRID_TRACE("PlaqP4");
	      autoView( gStencil_v, gStencils_plqforce[p6+3], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x       =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto Mdag_xmmupnu = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], (-rho) * U_nu_x * Mdag_xmmupnu);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_mu_xmmu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(3,ss);
		  auto U_nu_xmmu    =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqR_v[ss], Udag_mu_xmmu * U_nu_xmmu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);

	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    /////////////////////////////////////////////////////////////////////
	    // Set up the determinant force contribution in 3x3 algebra basis
	    /////////////////////////////////////////////////////////////////////
	    InsertForce(Fdet1,Fdet1_nu,nu);
	    InsertForce(Fdet2,Fdet2_nu,nu);

	    //////////////////////////////////////////////////
	    // Parallel direction terms
	    //////////////////////////////////////////////////

	    //    y..
	    //    |  |
	    //    x==   // mu polarisation
	    // PlaqL = -rho U_mu(x) U_nu^d(x+mu-nu) msk^d(x-nu) ; PlaqR = U_nu^d(x-nu)
	    {
	      GRID_TRACE("PlaqP5");
	      autoView( gStencil_v, gStencils_plqforce[p6+4], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_mu_x         =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto Udag_nu_xpmumnu= adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(2,ss);
		  auto Mdag_xmnu      = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], (-rho) * U_mu_x * Udag_nu_xpmumnu * Mdag_xmnu);
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xmnu   = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], Udag_nu_xmnu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;

	    // x==
	    // |  |
	    // y..          // mu polarisation

	    // PlaqL = -rho U_mu(x) U_nu(x+mu) msk^d(x+nu) ; PlaqR = U_nu(x)
	    {
	      GRID_TRACE("PlaqP6");
	      autoView( gStencil_v, gStencils_plqforce[p6+5], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_mu_x     =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpmu  =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Mdag_xpnu  = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], (-rho) * U_mu_x * U_nu_xpmu * Mdag_xpnu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = Umu[nu];
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,1);

	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;

	    break;
	  }
	case 2:
	  {
	    ///////////////// +ve nu /////////////////
	    //     ->
	    //    |  |
	    //    :  |
	    //    x==
	    // x = y

	    PlaqL=Ident;

	    // PlaqR = -rho U_nu(x) U_nu(x+nu) U_mu(x+2nu) U_nu^d(x+mu+nu) U_nu^d(x+mu) msk^d(x)
	    {
	      GRID_TRACE("PlaqR1");
	      autoView( gStencil_v, gStencils_rectforce[p10+0], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x         =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpnu      =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_mu_xp2nu     =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xpmupnu= adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(4,ss);
		  auto Udag_nu_xpmu   = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_x         = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_nu_x * U_nu_xpnu * U_mu_xp2nu * Udag_nu_xpmupnu * Udag_nu_xpmu * Mdag_x);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y =   dJdXe_nMpInv;
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;

	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx,FdetV);
	    Fdet2_nu = FdetV;

	    //     <-
	    //    |  |
	    //    |  :
	    //    x==y
	    // x = y - mu

	    // PlaqR = rho U_nu(x) U_nu(x+nu) U_mu^d(x-mu+2nu) U_nu^d(x-mu+nu) U_nu^d(x-mu)
	    // PlaqL = msk^d(x-mu)   [adj(PlaqL) is used in ComputeNxy & Compute_MpInvJx_dNxxdSy]
	    {
	      GRID_TRACE("PlaqR2");
	      autoView( gStencil_v, gStencils_rectforce[p10+1], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x           =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpnu        =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_mu_xmmup2nu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xmmupnu  = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(4,ss);
		  auto Udag_nu_xmmu     = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (rho) * U_nu_x * U_nu_xpnu * Udag_mu_xmmup2nu * Udag_nu_xmmupnu * Udag_nu_xmmu);
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xmmu        = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], Mdag_xmmu);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	      PlaqL = GhostRect.Extract(gPlaqL);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    //    :->
	    //    y  |
	    //    |  |
	    //    x==
	    // x = y - nu

	    // PlaqR = -rho U_nu(x) U_mu(x+nu) U_nu^d(x+mu) U_nu^d(x+mu-nu) msk^d(x-nu)
	    // PlaqL = U_nu^d(x-nu)
	    {
	      GRID_TRACE("PlaqR3");
	      autoView( gStencil_v, gStencils_rectforce[p10+2], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x          =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_mu_xpnu       =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_nu_xpmu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xpmumnu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(4,ss);
		  auto Mdag_xmnu       = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_nu_x * U_mu_xpnu * Udag_nu_xpmu * Udag_nu_xpmumnu * Mdag_xmnu);
		  e = gStencil_v.GetEntry(5,ss);
		  auto Udag_nu_xmnu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], Udag_nu_xmnu);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	      PlaqL = GhostRect.Extract(gPlaqL);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    //     <-:
	    //    |  y
	    //    |  |
	    //    x==
	    // x = y - mu - nu

	    // PlaqR = rho U_nu(x) U_mu^d(x-mu+nu) U_nu^d(x-mu) U_nu^d(x-mu-nu)
	    // PlaqL = U_nu^d(x-nu) msk^d(x-mu-nu)
	    {
	      GRID_TRACE("PlaqR4");
	      autoView( gStencil_v, gStencils_rectforce[p10+3], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x          =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto Udag_mu_xmmupnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(2,ss);
		  auto Udag_nu_xmmu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xmmumnu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (rho) * U_nu_x * Udag_mu_xmmupnu * Udag_nu_xmmu * Udag_nu_xmmumnu);
		  e = gStencil_v.GetEntry(4,ss);
		  auto Udag_nu_xmnu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xmmumnu    = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], Udag_nu_xmnu * Mdag_xmmumnu);
		});
	      PlaqR = GhostRect.Extract(gPlaqR);
	      PlaqL = GhostRect.Extract(gPlaqL);
	    }

	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,-1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,-1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    ///////////////// -ve nu /////////////////

	    //    x==
	    //    :  |
	    //    y
	    //    |  |
	    //     <-
	    // x = y + nu

	    // PlaqL = rho U_nu^d(x-nu) U_mu(x-nu) U_nu(x-nu+mu) U_nu(x+mu) msk^d(x+nu)
	    // PlaqR = U_nu(x)
	    {
	      GRID_TRACE("PlaqR5");
	      autoView( gStencil_v, gStencils_rectforce[p10+4], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto Udag_nu_xmnu   = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_mu_xmnu      =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_nu_xmnupmu   =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto U_nu_xpmu      =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(4,ss);
		  auto Mdag_xpnu      = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], (rho) * Udag_nu_xmnu * U_mu_xmnu * U_nu_xmnupmu * U_nu_xpmu * Mdag_xpnu);
		  e = gStencil_v.GetEntry(5,ss);
		  auto U_nu_x         =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqR_v[ss], U_nu_x);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,1);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    //    x==
	    //    |  :
	    //       y
	    //    |  |
	    //     ->
	    // x = y - mu + nu

	    // PlaqL = -rho U_nu^d(x-nu) U_mu^d(x-mu-nu) U_nu(x-mu-nu) U_nu(x-mu)
	    // PlaqR = U_nu(x) msk^d(x-mu+nu)
	    {
	      GRID_TRACE("PlaqR6");
	      autoView( gStencil_v, gStencils_rectforce[p10+5], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto Udag_nu_xmnu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(1,ss);
		  auto Udag_mu_xmmumnu = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_nu_xmmumnu    =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto U_nu_xmmu       =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqL_v[ss], (-rho) * Udag_nu_xmnu * Udag_mu_xmmumnu * U_nu_xmmumnu * U_nu_xmmu);
		  e = gStencil_v.GetEntry(4,ss);
		  auto U_nu_x          =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xmmupnu    = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], U_nu_x * Mdag_xmmupnu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,1);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,1);
	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    //    x==
	    //    |  |
	    //    :  |
	    //    y->
	    // x = y + 2nu

	    // PlaqL = U_mu(x) U_nu(x+mu) U_nu(x+mu+nu) msk^d(x+2nu)
	    // PlaqR = rho U_nu(x) U_nu(x+nu)
	    {
	      GRID_TRACE("PlaqR7");
	      autoView( gStencil_v, gStencils_rectforce[p10+6], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_mu_x        =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpmu     =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_nu_xpmupnu  =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto Mdag_xp2nu    = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], U_mu_x * U_nu_xpmu * U_nu_xpmupnu * Mdag_xp2nu);
		  e = gStencil_v.GetEntry(4,ss);
		  auto U_nu_x        =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(5,ss);
		  auto U_nu_xpnu     =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqR_v[ss], (rho) * U_nu_x * U_nu_xpnu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,2);
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,2);
	    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu + FdetV;

	    //    x==
	    //    |  |
	    //    |  :
	    //     ->y
	    // x = y + 2nu - mu

	    // PlaqL = U_mu^d(x-mu) U_nu(x-mu) U_nu(x-mu+nu)
	    // PlaqR = -rho U_nu(x) U_nu(x+nu) msk^d(x-mu+2nu)
	    {
	      GRID_TRACE("PlaqR8");
	      autoView( gStencil_v, gStencils_rectforce[p10+7], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto Udag_mu_xmmu   = adj(coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xmmu      =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_nu_xmmupnu   =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqL_v[ss], Udag_mu_xmmu * U_nu_xmmu * U_nu_xmmupnu);
		  e = gStencil_v.GetEntry(3,ss);
		  auto U_nu_x         =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(4,ss);
		  auto U_nu_xpnu      =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xmmup2nu  = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_nu_x * U_nu_xpnu * Mdag_xmmup2nu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,2);
	    ComputeNxy(PlaqL, PlaqR,Nxy);
	    Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,2);
	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_nu = Fdet2_nu+FdetV;

	    /////////////////////////////////////////////////////////////////////
	    // Set up the determinant force contribution in 3x3 algebra basis
	    /////////////////////////////////////////////////////////////////////
	    InsertForce(Fdet1,Fdet1_nu,nu);
	    InsertForce(Fdet2,Fdet2_nu,nu);

	    //////////////////////////////////////////////////
	    // Parallel direction terms
	    //////////////////////////////////////////////////

	    //    y..
	    //    |  |
	    //    |  |
	    //    x=<=
	    // x = y - 2nu

	    // PlaqL = U_nu^d(x-nu) U_nu^d(x-2nu)
	    // PlaqR = -rho U_mu(x) U_nu^d(x+mu-nu) U_nu^d(x+mu-2nu) msk^d(x-2nu)
	    {
	      GRID_TRACE("PlaqR9");
	      autoView( gStencil_v, gStencils_rectforce[p10+8], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto Udag_nu_xmnu     = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(1,ss);
		  auto Udag_nu_xm2nu    = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqL_v[ss], Udag_nu_xmnu * Udag_nu_xm2nu);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_mu_x           =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto Udag_nu_xpmumnu  = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(4,ss);
		  auto Udag_nu_xpmum2nu = adj(coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd));
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xm2nu       = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_mu_x * Udag_nu_xpmumnu * Udag_nu_xpmum2nu * Mdag_xm2nu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-2);
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,-2);
	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;

	    //    x<=
	    //    |  |
	    //    |  |
	    //    y..
	    // x = y + 2nu

	    // PlaqL = U_nu(x) U_nu(x+nu)
	    // PlaqR = -rho U_mu(x) U_nu(x+mu) U_nu(x+mu+nu) msk^d(x+2nu)
	    {
	      GRID_TRACE("PlaqR10");
	      autoView( gStencil_v, gStencils_rectforce[p10+9], AcceleratorRead);
	      accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
		  GeneralStencilEntry const* e = gStencil_v.GetEntry(0,ss);
		  auto U_nu_x        =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(1,ss);
		  auto U_nu_xpnu     =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  coalescedWrite(gPlaqL_v[ss], U_nu_x * U_nu_xpnu);
		  e = gStencil_v.GetEntry(2,ss);
		  auto U_mu_x        =     coalescedReadGeneralPermute(gU_mu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(3,ss);
		  auto U_nu_xpmu     =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(4,ss);
		  auto U_nu_xpmupnu  =     coalescedReadGeneralPermute(gU_nu_v[e->_offset], e->_permute, Nd);
		  e = gStencil_v.GetEntry(5,ss);
		  auto Mdag_xp2nu    = adj(coalescedReadGeneralPermute(gUtmp_v[e->_offset], e->_permute, Nd));
		  coalescedWrite(gPlaqR_v[ss], (-rho) * U_mu_x * U_nu_xpmu * U_nu_xpmupnu * Mdag_xp2nu);
		});
	      PlaqL = GhostRect.Extract(gPlaqL);
	      PlaqR = GhostRect.Extract(gPlaqR);
	    }

	    dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,2);
	    ComputeNxy(PlaqL,PlaqR,Nxy);
	    Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	    MpInvJx_nu = Cshift(MpInvJx,nu,2);
	    Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);
	    Fdet2_mu = Fdet2_mu+FdetV;
	    break;
	  }
	default:
	  {
	    assert(1!=1 && " At present, the only valid choice of flow kernel is either 1 or 2");
	    break;
	  }
	}
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    // Sign conventions as in the reference routine above
    force=-1.0*(Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogPerformance << " logDetJacobianForceLevelOpt took "<<t1-t0<<" us"
	      << " (prelim "<<t3a-t0<<" us, dJdXe "<<t4-t3a<<" us, nu loop "<<t5-t4<<" us)"<<std::endl;
  }

  // Old top-level, kept for consistency checks: identical chain rule to the
  // default logDetJacobianForce, but with the old level routine. To be
  // deleted together with the old level routines.
  void logDetJacobianForce(int old, GaugeField &force)
  {
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond();

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);

	// Add back U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = peekLorentz(this->get_smeared_conf(ismr - 1), mu) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

	// Get this levels determinant force
	force_det = Zero();
	logDetJacobianForceLevel(old,this->get_smeared_conf(ismr-1),force_det,ismr);

	// Sum the contributions
	force = force + force_det;
      }

      // remove U in UdSdU...
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = adj(peekLorentz(this->get_smeared_conf(0), mu)) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);

      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(old,*this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta

      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationRect: lnDetJacobianForce(old) took " << time << " ms" << std::endl;
    }  // if smearingLevels = 0 do nothing
  }

  // Default (optimised) implementation. Same result as the old overload
  // below: fused Ncb via the two-argument LieAlgebraProject (with PlaqL = 1
  // it reproduces Nb = 2 Ta(i T^b U C^dag) exactly), Zac via
  // make_adjoint_rep, and J Taylor + Mab + Determinant + log fused in one
  // sitewise kernel (cf. the plaquette version in GaugeConfigurationMasked.h,
  // which runs on the half grid; here full grid + mask before the sum —
  // compression of the masked sum is a later lever).
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
  {
    GRID_TRACE("logDetJacobianLevel");
    GridBase* grid = U.Grid();
    GaugeLinkField Umu(grid), Cmu(grid), PlaqL(grid);
    GaugeLinkField Z(grid);
    AdjMatrixField  Ncb(grid);
    AdjMatrixField  Zac(grid);
    LatticeComplex ln_det(grid);
    ColourMatrix Ident;

    int mu= (smr/2) %Nd; // both smearing types are of 2 colouring
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu);
    Ident = ComplexD(1.0);

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    int smr_ind = smr/Nsmr_one_step;
    int flw_knl = mask_types[smr_ind];
    double rho;
    switch(flw_knl){
    case 1:
      rho=this->Stouts[smr_ind]->SmearRho[1];
      break;
    case 2:
      rho=((Rect_Stout<Gimpl> *) this->Stouts[smr_ind])->SmearRhoRs[1];
      break;
    }

    {
      GRID_TRACE("ExchangePeriodicRect");
      GaugeField gU(grid);
      gU = GhostRect.ExchangePeriodic(U);
      switch(flw_knl){
      case 1:
	BaseSmear_ghost_plq(Cmu, gU, mu, rho);
	break;
      case 2:
	BaseSmear_ghost_rect(Cmu, gU, mu, rho);
	break;
      }
    }

    Umu = peekLorentz(U, mu);
    PlaqL = Ident;
    ComputeNxy(PlaqL, Umu*adj(Cmu), Ncb);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra; move to the adjoint rep
    Z  = Ta(Cmu * adj(Umu));
    SU3Adjoint::make_adjoint_rep(Zac, Z);

    //////////////////////////////////////////////////////////////////
    // J(x) = 1 + Sum_k (-Zac)^k/(k+1)!, Mab, det, log: one kernel
    //////////////////////////////////////////////////////////////////
    {GRID_TRACE("J_Mab_lnDet");
      autoView(ln_det_v,ln_det,AcceleratorWrite);
      autoView(Zac_v,Zac,AcceleratorRead);
      autoView(Ncb_v,Ncb,AcceleratorRead);
      accelerator_for(ss,grid->oSites(),grid->Nsimd(),{
	  typedef decltype(coalescedRead(Zac_v(0)))    adj_mat;
	  adj_mat X, Jac, Mab_ss;
	  RealD kpfac = 1;

	  X=1.0;
	  Jac = X;
	  for(int k=1;k<12;k++){
	    X=(-1.0)*X*Zac_v(ss);
	    kpfac = kpfac /((RealD) (k+1));
	    Jac = Jac + X * kpfac;
	  }

	  Mab_ss = Complex(1.0,0.0);
	  Mab_ss = Mab_ss - Jac * Ncb_v(ss);

	  auto detD = Determinant(Mab_ss);
	  coalescedWrite(ln_det_v[ss],log(detD));
	});
    }

    ////////////////////////////
    // Masked sum
    ////////////////////////////
    ln_det = ln_det * mask;
    Complex result = sum(ln_det);
    return result.real();
  }

  // Old implementation, kept for consistency checks against the optimised
  // default above; to be deleted once confirmed. The int old argument only
  // selects this overload.
  RealD logDetJacobianLevel(int old, const GaugeField &U,int smr)
  {
    GridBase* grid = U.Grid();
    GaugeField C(grid);
    GaugeLinkField Nb(grid);
    GaugeLinkField Z(grid);
    GaugeLinkField Umu(grid), Cmu(grid);
    ColourMatrix   Tb;
    ColourMatrix   Tc;
    typedef typename SU3Adjoint::AMatrix AdjMatrix;
    typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
    typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    LatticeComplex  cplx(grid); 
    AdjVectorField  AlgV(grid); 
    AdjMatrixField  Mab(grid);
    AdjMatrixField  Ncb(grid);
    AdjMatrixField  Jac(grid);
    AdjMatrixField  Zac(grid);
    AdjMatrixField  mZac(grid);
    AdjMatrixField  X(grid);

    int mu= (smr/2) %Nd; // both smearing types are of 2 colouring

    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu);

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    int smr_ind = smr/Nsmr_one_step;
    int flw_knl = mask_types[smr_ind];
    double rho;
    switch(flw_knl){
    case 1:
      rho=this->Stouts[smr_ind]->SmearRho[1];
      break;
    case 2:
      rho=((Rect_Stout<Gimpl> *) this->Stouts[smr_ind])->SmearRhoRs[1];
      break;
    }
    BaseSmear(Cmu, U,mu,rho,flw_knl);

    std::cout << GridLogMessage << "logDetJacobianLevel " << mu <<" "<< rho<<" "<<flw_knl<<" "<<smr_ind <<std::endl;
    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, Tb);
      // Qlat Tb = 2i Tb^Grid
      Nb = (2.0)*Ta( ci*Tb * Umu * adj(Cmu));
      // FIXME -- replace this with LieAlgebraProject
#if 0
      SU3::LieAlgebraProject(Ncb,tmp,b);
#else
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, Tc);
	auto tmp = -trace(ci*Tc*Nb); // Luchang's norm: (2Tc) (2Td) N^db = -2 delta cd N^db // - was important
	PokeIndex<ColourIndex>(Ncb,tmp,c,b); 
      }
#endif
    }      

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Z  = Ta(Cmu * adj(Umu));

    // Move Z to the Adjoint Rep == make_adjoint_representation
    Zac = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      // Is the mapping of these the same? Same structure constants
      // Might never have been checked.
      SU3::generator(b, Tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*Tb*Z); // my convention 1/2 delta ba
      Zac = Zac + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    X=1.0; 
    Jac = X;
    mZac = (-1.0)*Zac; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZac;
      kpfac = kpfac /(k+1);
      Jac = Jac + X * kpfac;
    }

    ////////////////////////////
    // Mab
    ////////////////////////////
    Mab = Complex(1.0,0.0);
    Mab = Mab - Jac * Ncb;

    ////////////////////////////
    // det
    ////////////////////////////
    LatticeComplex       det(grid); 
    det = Determinant(Mab);

    ////////////////////////////
    // ln det
    ////////////////////////////
    LatticeComplex       ln_det(grid); 
    ln_det = log(det);

    ////////////////////////////
    // Masked sum
    ////////////////////////////
    ln_det = ln_det * mask;
    Complex result = sum(ln_det);
    return result.real();
  }
public:
  virtual RealD logDetJacobian(void)
  {
    RealD ln_det = 0;
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianLevel(this->get_smeared_conf(ismr-1),ismr);
      }
      ln_det +=logDetJacobianLevel(*(this->ThinLinks),0);

      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationRect: logDetJacobian took " << time << " ms" << std::endl;
    }
    return ln_det;
  }
  // Old top-level, kept for consistency checks; to be deleted together with
  // the old level routines.
  RealD logDetJacobian(int old)
  {
    RealD ln_det = 0;
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
	ln_det+= logDetJacobianLevel(old,this->get_smeared_conf(ismr-1),ismr);
      }
      ln_det +=logDetJacobianLevel(old,*(this->ThinLinks),0);

      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationRect: logDetJacobian(old) took " << time << " ms" << std::endl;
    }
    return ln_det;
  }
  virtual void logDetJacobianForce(GaugeField &force)
  {
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond();

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
	
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);

	// Add back U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = peekLorentz(this->get_smeared_conf(ismr - 1), mu) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
    	
	// Get this levels determinant force
	force_det = Zero();
	logDetJacobianForceLevel(this->get_smeared_conf(ismr-1),force_det,ismr);

	// Sum the contributions
	force = force + force_det;
      }
    
      // remove U in UdSdU...
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = adj(peekLorentz(this->get_smeared_conf(0), mu)) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);

      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(*this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationRect: lnDetJacobianForce took " << time << " ms" << std::endl;  
    }  // if smearingLevels = 0 do nothing
  }

public:
  //====================================================================
  // Override base clas here to mask it
  virtual void fill_smearedSet(GaugeField &U)
  {
    this->ThinLinks = &U;  // attach the smearing routine to the field U

    // check the pointer is not null
    if (this->ThinLinks == NULL)
      std::cout << GridLogError << "[SmearedConfigurationRect] Error in ThinLinks pointer\n";

    if (this->smearingLevels > 0)
    {
      std::cout << GridLogMessage << "[SmearedConfigurationRect] Filling SmearedSet\n";
      GaugeField previous_u(this->ThinLinks->Grid());

      GaugeField smeared_A(this->ThinLinks->Grid());
      GaugeField smeared_B(this->ThinLinks->Grid());
      //std::cout << GridLogDebug << this->smearingLevels <<" "<<this->SmearedSet.size()<<""<<Nsmr_one_step<<std::endl;//DEBUG
      previous_u = *this->ThinLinks;
      double start = usecond();
      for (int smearLvl = 0; smearLvl < this->smearingLevels; smearLvl+=Nsmr_one_step)
	for(int smr=0; smr<Nsmr_one_step; smr++) {
	  int smr_ind = smearLvl/Nsmr_one_step;
	  int flw_knl = mask_types[smr_ind];
	  this->Stouts[smr_ind]->smear(smeared_A, previous_u);
	  ApplyMask(smeared_A,smearLvl+smr);
	  smeared_B = previous_u;
	  ApplyMask(smeared_B,smearLvl+smr);
	  // Replace only the masked portion
	  this->SmearedSet[smearLvl+smr] = previous_u-smeared_B + smeared_A;
	  previous_u = this->SmearedSet[smearLvl+smr];

	  // For debug purposes
	  RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(previous_u);
	  std::cout << GridLogMessage << "[SmearedConfigurationRect] smeared Plaq: " << impl_plaq << std::endl;
	}
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationRect: Link smearing took " << time << " ms" << std::endl;  
    }
  }
  //====================================================================
  // Override base to add masking
  virtual GaugeField AnalyticSmearedForce(const GaugeField& SigmaKPrime,
					  const GaugeField& GaugeK,int level) 
  {
    GridBase* grid = GaugeK.Grid();
    GaugeField SigmaK(grid), iLambda(grid);
    GaugeField SigmaKPrimeA(grid);
    GaugeField SigmaKPrimeB(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);

    int smr_ind   = level/Nsmr_one_step; // which smearer to be used
    int mask_type = mask_types[smr_ind]; // what type of the smearer corresp. to smr_ind
    int sub_level = level%Nsmr_one_step; // Each mask type uses len(red,black)*Nd = 2*4 = Nsmr_one_step steps
    int mmu= (sub_level/2) %Nd;
    int cb= (sub_level%2);
    double rho=0;
    switch(mask_type){
    case 1:
      rho = this->Stouts[smr_ind]->SmearRho[1];
      break;
    case 2:
      rho = ((Rect_Stout<Gimpl> *)this->Stouts[smr_ind])->SmearRhoRs[1];
      break;
    }

    // Can override this to do one direction only.
    SigmaK = Zero();
    iLambda = Zero();

    SigmaKPrimeA = SigmaKPrime;
    ApplyMask(SigmaKPrimeA,level);
    SigmaKPrimeB = SigmaKPrime - SigmaKPrimeA; // not updated
    // Could get away with computing only one polarisation here
    // int mu= (smr/2) %Nd;
    // SigmaKprime_A has only one component
    //    GaugeField C(grid);
    //    this->StoutSmearing->BaseSmear(C, GaugeK);
    //    for (int mu = 0; mu < Nd; mu++)
    int mu =mmu;
    BaseSmear(Cmu, GaugeK,mu,rho,mask_type);

    {
      GaugeKmu = peekLorentz(GaugeK, mu);
      SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
      iQ = Ta(Cmu * adj(GaugeKmu));
      this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu); 
      pokeLorentz(iLambda, iLambda_mu, mu);
      std::cout << " mu "<<mu<<" SigmaKPrime_mu "<<norm2(SigmaKPrime_mu)<< " iLambda_mu " <<norm2(iLambda_mu)<<std::endl;
    }
    //    GaugeField SigmaKcopy(grid);
    //    SigmaKcopy = SigmaK;
    BaseSmearDerivative(SigmaK, iLambda,GaugeK,mu,rho,mask_type);  // derivative of SmearBase
    //    this->StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase
    //    SigmaKcopy = SigmaKcopy - SigmaK;
    //    std::cout << " BaseSmearDerivative fast path error" <<norm2(SigmaKcopy)<<std::endl;
    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;

    return SigmaK;
  }

public:

  /* Standard constructor */
  SmearedConfigurationRect(GridCartesian* _UGrid, unsigned int Nsmear, std::vector<Smear_Stout<Gimpl> *> Stouts, std::vector<int> mask_types={2,1})
    : SmearedConfigurationMasked<Gimpl>(_UGrid, Nsmear,*Stouts[0]), Stouts(Stouts), mask_types(mask_types),
      GhostRect(2,_UGrid)
  {
    assert(Nsmear%(Nsmr_one_step)==0);
    assert(Nsmear/Nsmr_one_step==mask_types.size()); // Nsmr_one_step = #Basic_num_steps = 2*Nd = 8; One step <-> a mask type
    assert(Stouts.size() == mask_types.size());

    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);

    bool has_rect = false;
    for (auto mt : mask_types) if (mt == 2) has_rect = true;
    if (has_rect) {
      // The type-2 mask requires every extent = 0 mod 4: the 2x2 blocks must tile the
      // torus with an even number of blocks per perp. direction for the red-black block
      // colouring to alternate across the periodic boundary, and every direction takes
      // the perp. role at some level.
      Coordinate gdims = _UGrid->GlobalDimensions();
      for (int d = 0; d < Nd; d++) assert(gdims[d] % 4 == 0);

      // Stencils for the optimised rect staple, on the depth-2 padded grid.
      GridBase *ggrid = GhostRect.grids[Nd-1];
      for (int mu = 0; mu < Nd; mu++)
	gStencils_rectsmear.push_back(Rect_Stout<Gimpl>::RectStapleStencilRs(ggrid, mu));
    }

    ///////////////////////////////////////////////////////////////////////
    // Stencils for the force-level PlaqL/PlaqR terms (depth-2 padded grid).
    // Shift lists transcribe the CovShift chains of the old
    // logDetJacobianForceLevel term by term; the entry order must match the
    // kernel read order there. "msk" marks entries read from the padded
    // MASKED link gUtmp (this class cannot use the checkerboard trick of
    // GaugeConfigurationMasked.h: the type-2 mask is not a parity class).
    ///////////////////////////////////////////////////////////////////////
    {
      bool has_plq = false;
      for (auto mt : mask_types) if (mt == 1) has_plq = true;
      GridBase *ggrid = GhostRect.grids[Nd-1];
      std::vector<Coordinate> shifts;

      // plq staple stencil (one per mu; 6 entries per nu!=mu), same
      // geometry as gStencils_smear in GaugeConfigurationMasked.h
      if (has_plq) {
	for(int mu=0;mu<Nd;mu++){
	  Coordinate shift_0(Nd,0);
	  shifts.clear();
	  for(int nu=0;nu<Nd;nu++){
	    if (nu==mu) continue;
	    Coordinate shift_nu(Nd,0);  shift_nu[nu]=1;
	    Coordinate shift_mu(Nd,0);  shift_mu[mu]=1;
	    Coordinate shift_mnu(Nd,0); shift_mnu[nu]=-1;
	    Coordinate shift_pmu_mnu(Nd,0); shift_pmu_mnu[mu]=1; shift_pmu_mnu[nu]=-1;
	    // upper: U_nu(x) U_mu(x+nu) U_nu^d(x+mu)
	    shifts.push_back(shift_0); shifts.push_back(shift_nu); shifts.push_back(shift_mu);
	    // lower: U_nu^d(x-nu) U_mu(x-nu) U_nu(x+mu-nu)
	    shifts.push_back(shift_mnu); shifts.push_back(shift_mnu); shifts.push_back(shift_pmu_mnu);
	  }
	  gStencils_plqsmear.push_back(GeneralLocalStencil(ggrid,shifts));
	}
      }

      for(int mu=0;mu<Nd;mu++){
	for(int nu=0;nu<Nd;nu++){
	  if (nu==mu) continue;
	  auto S = [mu,nu](int amu,int anu){ Coordinate c(Nd,0); c[mu]=amu; c[nu]=anu; return c; };

	  if (has_plq) {
	    // P1 (+nu cw, x=y):    R: U_nu(x) U_mu(x+nu) U_nu^d(x+mu) msk^d(x); L=1
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(1,0)); shifts.push_back(S(0,0));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // P2 (+nu acw, x=y-mu): R: U_nu(x) U_mu^d(x-mu+nu) U_nu^d(x-mu); L: msk^d(x-mu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(-1,1)); shifts.push_back(S(-1,0)); shifts.push_back(S(-1,0));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // P3 (-nu cw, x=y+nu):  L: U_mu(x) U_nu(x+mu) msk^d(x+nu); R = U_nu(x) direct
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(1,0)); shifts.push_back(S(0,1));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // P4 (-nu acw, x=y-mu+nu): L: U_nu(x) msk^d(x-mu+nu); R: U_mu^d(x-mu) U_nu(x-mu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(-1,1)); shifts.push_back(S(-1,0)); shifts.push_back(S(-1,0));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // P5 (mu pol, x=y-nu): L: U_mu(x) U_nu^d(x+mu-nu) msk^d(x-nu); R: U_nu^d(x-nu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(1,-1)); shifts.push_back(S(0,-1)); shifts.push_back(S(0,-1));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // P6 (mu pol, x=y+nu): L: U_mu(x) U_nu(x+mu) msk^d(x+nu); R = U_nu(x) direct
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(1,0)); shifts.push_back(S(0,1));
	    gStencils_plqforce.push_back(GeneralLocalStencil(ggrid,shifts));
	  }
	  if (has_rect) {
	    // R1 (+nu, x=y): R: U_nu(x) U_nu(x+nu) U_mu(x+2nu) U_nu^d(x+mu+nu) U_nu^d(x+mu) msk^d(x); L=1
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(0,2));
	    shifts.push_back(S(1,1)); shifts.push_back(S(1,0)); shifts.push_back(S(0,0));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R2 (x=y-mu): R: U_nu(x) U_nu(x+nu) U_mu^d(x-mu+2nu) U_nu^d(x-mu+nu) U_nu^d(x-mu); L: msk^d(x-mu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(-1,2));
	    shifts.push_back(S(-1,1)); shifts.push_back(S(-1,0)); shifts.push_back(S(-1,0));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R3 (x=y-nu): R: U_nu(x) U_mu(x+nu) U_nu^d(x+mu) U_nu^d(x+mu-nu) msk^d(x-nu); L: U_nu^d(x-nu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(1,0));
	    shifts.push_back(S(1,-1)); shifts.push_back(S(0,-1)); shifts.push_back(S(0,-1));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R4 (x=y-mu-nu): R: U_nu(x) U_mu^d(x-mu+nu) U_nu^d(x-mu) U_nu^d(x-mu-nu); L: U_nu^d(x-nu) msk^d(x-mu-nu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(-1,1)); shifts.push_back(S(-1,0));
	    shifts.push_back(S(-1,-1)); shifts.push_back(S(0,-1)); shifts.push_back(S(-1,-1));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R5 (-nu, x=y+nu): L: U_nu^d(x-nu) U_mu(x-nu) U_nu(x-nu+mu) U_nu(x+mu) msk^d(x+nu); R: U_nu(x)
	    shifts.clear();
	    shifts.push_back(S(0,-1)); shifts.push_back(S(0,-1)); shifts.push_back(S(1,-1));
	    shifts.push_back(S(1,0)); shifts.push_back(S(0,1)); shifts.push_back(S(0,0));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R6 (x=y-mu+nu): L: U_nu^d(x-nu) U_mu^d(x-mu-nu) U_nu(x-mu-nu) U_nu(x-mu); R: U_nu(x) msk^d(x-mu+nu)
	    shifts.clear();
	    shifts.push_back(S(0,-1)); shifts.push_back(S(-1,-1)); shifts.push_back(S(-1,-1));
	    shifts.push_back(S(-1,0)); shifts.push_back(S(0,0)); shifts.push_back(S(-1,1));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R7 (x=y+2nu): L: U_mu(x) U_nu(x+mu) U_nu(x+mu+nu) msk^d(x+2nu); R: U_nu(x) U_nu(x+nu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(1,0)); shifts.push_back(S(1,1));
	    shifts.push_back(S(0,2)); shifts.push_back(S(0,0)); shifts.push_back(S(0,1));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R8 (x=y+2nu-mu): L: U_mu^d(x-mu) U_nu(x-mu) U_nu(x-mu+nu); R: U_nu(x) U_nu(x+nu) msk^d(x-mu+2nu)
	    shifts.clear();
	    shifts.push_back(S(-1,0)); shifts.push_back(S(-1,0)); shifts.push_back(S(-1,1));
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(-1,2));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R9 (mu pol, x=y-2nu): L: U_nu^d(x-nu) U_nu^d(x-2nu); R: U_mu(x) U_nu^d(x+mu-nu) U_nu^d(x+mu-2nu) msk^d(x-2nu)
	    shifts.clear();
	    shifts.push_back(S(0,-1)); shifts.push_back(S(0,-2)); shifts.push_back(S(0,0));
	    shifts.push_back(S(1,-1)); shifts.push_back(S(1,-2)); shifts.push_back(S(0,-2));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	    // R10 (mu pol, x=y+2nu): L: U_nu(x) U_nu(x+nu); R: U_mu(x) U_nu(x+mu) U_nu(x+mu+nu) msk^d(x+2nu)
	    shifts.clear();
	    shifts.push_back(S(0,0)); shifts.push_back(S(0,1)); shifts.push_back(S(0,0));
	    shifts.push_back(S(1,0)); shifts.push_back(S(1,1)); shifts.push_back(S(0,2));
	    gStencils_rectforce.push_back(GeneralLocalStencil(ggrid,shifts));
	  }
	}
      }
    }

    ////////////////////
    // Setup the mask
    ////////////////////
    GridRedBlackCartesian * UrbGrid;
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    //std::vector<Lattice<iScalar<vInteger> > > xs(3,_UGrid);
    LatticeComplex zeros(_UGrid); zeros = Zero();
    LatticeComplex ones(_UGrid); ones = ComplexD(1.0,0.0); 
    LatticeComplex tmp(_UGrid);
    
    for (unsigned int i = 0; i < this->smearingLevels; i+=Nsmr_one_step) {
      int mask_type = mask_types[i/Nsmr_one_step];
      for(unsigned int j = 0; j < Nsmr_one_step; ++j) {
	masks.push_back(*(new LatticeLorentzComplex(_UGrid)));
	
	masks[i+j]=Zero();
	tmp = Zero();
	
	int mu= (j/2) %Nd;
	int cb= (j%2);
	
	switch (mask_type) {
	case 1:
	  {
	    LatticeComplex tmpcb(UrbGrid);
	    
	    pickCheckerboard(cb,tmpcb,ones);
	    setCheckerboard(tmp,tmpcb);
	    break;
	  }
	case 2:
	  {
	    LatticeInteger  coor_nu(_UGrid),coor_sum(_UGrid); coor_sum = Zero();

	    for(int nu=0; nu < Nd; nu++){
	      LatticeCoordinate(coor_nu,nu);
	      if( nu != mu) coor_nu = div(coor_nu,mask_type);
	      coor_sum = coor_sum + coor_nu;
	    }
	    tmp = where( mod(coor_sum,2)==(Integer)(cb), ones, zeros);

#if 0 //DEBUG
	    std::cout << GridLogMessage <<"norm mask "<<norm2(tmp)<<std::endl;
	    
	    LatticeComplex tmp2(_UGrid), tmp3(_UGrid);
	    tmp2 = Cshift(tmp,(mu+1+Nd)%Nd,2);
	    tmp3 = Cshift(tmp,(mu+1+Nd)%Nd,1);
	    Coordinate point1(Nd), point2(Nd), point3(Nd);
	    std::cout << GridLogMessage <<"norm mask "<<norm2(tmp)<<" "<<norm2(ones)<<" "<<norm2(tmp+tmp2)<<" "<<norm2(tmp+tmp3)<<std::endl;
	    for(int nu=0;nu<Nd;nu++) {
	      point1[nu] = 2;
	      point3[nu] = 0;
	    }
	    
	    point2 = point1; point2[mu] = 3;
	    Complex c;
	    peekSite(c,tmp,point1);
	    std::cout << GridLogMessage <<"mask value:  mu="<<mu<<"cb= "<<cb<<" "<<point1<<" "<<c;
	    peekSite(c,tmp,point2); std::cout << point2<< " "<<c;
	    peekSite(c,tmp,point3); std::cout << point3<<" "<<c<<std::endl;
#endif	  
	  break;
	}
	}
	PokeIndex<LorentzIndex>(masks[i+j], tmp, mu);
      }
    }
    delete UrbGrid;
  }
  
  virtual void smeared_force(GaugeField &SigmaTilde) 
  {
    if (this->smearingLevels > 0)
    {
      double start = usecond();
      GaugeField force = SigmaTilde; // actually = U*SigmaTilde
      GaugeLinkField tmp_mu(SigmaTilde.Grid());

      // Remove U from UdSdU
      for (int mu = 0; mu < Nd; mu++)
      {
        // to get just SigmaTilde
        tmp_mu = adj(peekLorentz(this->SmearedSet[this->smearingLevels - 1], mu)) * peekLorentz(force, mu);
        pokeLorentz(force, tmp_mu, mu);
      }

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1),ismr);
      }
      
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);

      // Add U to UdSdU
      for (int mu = 0; mu < Nd; mu++)
      {
        tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
        pokeLorentz(SigmaTilde, tmp_mu, mu);
      }


      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << " GaugeConfigurationRect: Smeared Force chain rule took " << time << " ms" << std::endl;

    }  // if smearingLevels = 0 do nothing
    SigmaTilde=Gimpl::projectForce(SigmaTilde); // Ta
  }

};

NAMESPACE_END(Grid);

