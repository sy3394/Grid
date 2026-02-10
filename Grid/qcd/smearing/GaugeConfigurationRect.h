
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
template <class Gimpl>
class SmearedConfigurationRect : public SmearedConfigurationMasked<Gimpl> //from Masked better???
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

private:
  
  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;

  // These live in base class
  //  const unsigned int smearingLevels;
  //  Smear_Stout<Gimpl> *StoutSmearing;
  //  std::vector<GaugeField> SmearedSet;

  // Conventions:
  //   - flow kernel: 1 =:= Wilson,    2 =:= short-side rectangle
  //   -  mask types: 1 =:= red-black, 2 =:= 2x2 red-black in the plane perp. to \mu but identical in mu-dir
  
  int Nsmr_one_step = 2*Nd; // = #filterings(even colour, odd colour) x #dirs of smearing
  std::vector<int> mask_types;
  std::vector<Smear_Stout<Gimpl> *> Stouts;
  std::vector<LatticeLorentzComplex> masks; // should we turn this to poiners?????????

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
    case 1:{
      GaugeLinkField tmp_stpl(grid);
	  
      std::cout << GridLogMessage <<"BaseSmear: 1" << std::endl;//DEBUG
      WL.Staple(Cmu, U, mu);  //nb staple conventions of IroIro and Grid differ by a dagger
      Cmu = adj(rho * Cmu);
#if 0 // produced the exactly same numbers
      for(int nu=0; nu<Nd; ++nu){
	if (nu != mu) {
	// get the staple in direction mu, nu
        WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dagger
        Cmu += adj(tmp_stpl*rho);
	}
      }
#endif
      break;
    }
    case 2:
      // TODO: prepare optimized version
      std::cout << GridLogMessage <<"BaseSmear: 2" << std::endl;//DEBUG     
      WL.RectStapleUnoptimisedRs(Cmu, U, mu);
          Cmu = adj(rho * Cmu);
      break;
    }

    //    Cmu = adj(rho * Cmu);
    std::cout << GridLogMessage << "BaseSmear: " << mu<<" "<<rho<<" "<<flow_kernel<<" "<<norm2(Cmu) << std::endl;//DEBUG

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

    //Rect_Stout Smearer(grid);
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
      Fdet_pol=Fdet_pol + ci*tmp*te; // Fdet_pol + ci*tmp*te
    }
    pokeLorentz(Fdet, Fdet_pol, nu);
  }

  // tmp comment: no extra factor 
  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, tb);
      tb = 2.0 * ci * tb; // - ci * tb; in Lucher's convention but multiplied the missing factor from below
      Nx = Ta( adj(PlaqL)*tb * PlaqR );
      SU3::LieAlgebraProject(NxAd,Nx,b); /// needs to multiply -2 to get to Lucher's convention
    }
  }

  // tmp comment: orig. extra factor of (-2)*(-2)/(-2) = -2 <- multiplied the result by 2 but still deviation from Luscher by -1
  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
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
    std::cout << GridLogPerformance << " Compute_MpInvJx_dNxxdSy " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
	      << " ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms"<<std::endl;
  }

  void linkTracer(const std::vector<GaugeLinkField> &Umu, const GaugeLinkField &Umskd, const std::vector<int> dirs0, int ind, Real rho, GaugeLinkField &rect){
    // dir in dirs is 1+mu where mu=0,..3 to put sign on dir
    // ind: index of dirs for which masked gauge link field should be used.  The indexing starts from 1, not 0
    //    actually ind is always equal to size of dirs0 => can become just a flag
    GaugeLinkField tmp(Umu[0].Grid()), U(Umu[0].Grid()), tmp1(Umu[0].Grid());
    tmp = Zero(); U = Zero();
    std::vector<int> dirs = dirs0; 
    std::reverse(dirs.begin(), dirs.end());for(auto i :dirs)std::cout<<i<<std::endl;
    for(int i=0; i<dirs.size(); i++){
      int mu = dirs[i];
      //-1: -2 0 0 -1 12288 0
      std::cout << GridLogMessage <<"in_linkTracer "<<mu<<": "<<mu-1<<" "<<-mu-1<<" "<<norm2(tmp)<<" "<<ind<<" "<<norm2(U)<<" "<<std::abs(mu) - 1<<std::endl;
      if(i == 0){
	if ( ind>0) U = Umskd;
	else U = Umu[std::abs(mu) - 1];
	if(mu>0)
	  tmp = Gimpl::CovShiftIdentityForward(U,mu-1);
	else
	  tmp = Gimpl::CovShiftIdentityBackward(U,-mu-1);
      }else{
	if(mu>0)
	  tmp1 = Gimpl::CovShiftForward(Umu[std::abs(mu) - 1],mu-1,tmp);
	else
	  tmp1 = Gimpl::CovShiftBackward(Umu[std::abs(mu) - 1],-mu-1,tmp);
	tmp = tmp1;
      }
    }
    std::cout << GridLogMessage <<"linkTracer "<<dirs.size()<<" "<<rho<<" "<<norm2(tmp)<<std::endl;
    rect = rho*tmp;
  }
    
public:

  void logDetJacobianForceLevel(const GaugeField &U, GaugeField &force ,int smr)
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

#if 1  //testing combined Horner's rule
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
    
    // 2. The Horner Loop (Unchanged)
    for (int j = 12; j >= 1; --j) {
      double kfac = 1.0 / (j + 1);

      for(int b=0; b<8; b++) {
        // dJdX[b] correctly accumulates the negative contributions
        XB_dJdX[b] = (TRb_s[b] * XB_t3 + X * XB_dJdX[b] )* kfac;
      }

      XB_t3 = (XB_t2 * kfac) + aunit;
      XB_t2 = X * XB_t3;
    }
    XB_JxAd = XB_t3;
    std::cout << GridLogMessage << "DEBUG: Horner's method JxAd"<<norm2(XB_JxAd-JxAd)<< " djdx ";
    for(int i =0;i<dJdX.size();i++) std::cout<<norm2(XB_dJdX[i] - dJdX[i])<<" ";
    std::cout <<std::endl;
#endif
    
    time+=usecond();
    std::cout << GridLogMessage << "dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(PlaqL,PlaqR,NxxAd);
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

    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
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
#if 0
    std::vector<GaugeLinkField> Umu_tmp(Nd,grid);
    Umu_tmp = Umu;
    Umu_tmp[mu] = Utmp;//vector of gauge links in each dir where Umu_tmp[mu] is masked
    for(auto ttt: Umu_tmp) std::cout << GridLogMessage << "DEBUG: Umu_tmp[mu] "<< norm2(ttt)<<" "<<std::endl;
#endif
    
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
	case 1:{std::cout << GridLogMessage << "DEBUG: entered loops"<<std::endl;

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

#if 1 //DEBUG
	GaugeLinkField PlaqR2(grid);
	dirs = {(nu+1),mu+1,-(nu+1),-(mu+1)};
	linkTracer(Umu, Utmp, dirs, 4, -rho, PlaqR2);
	std::cout << GridLogMessage << "DEBUG: PlaqR linkT: "<<norm2(PlaqR2-PlaqR)<<std::endl;
#endif
	
	time=-usecond();
	dJdXe_nMpInv_y =   dJdXe_nMpInv;
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond();
	std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu = FdetV;
	time+=usecond();
	std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//     __
	//    |  :
	//    x==y    // nu polarisation -- anticlockwise

	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));
#if 1 //DEBUG
        //GaugeLinkField PlaqR2(grid);
        dirs = {(nu+1),-(mu+1),-(nu+1)};
        linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR2);
        std::cout << GridLogMessage << "DEBUG: PlaqR linkT: "<<norm2(PlaqR2-PlaqR)<<std::endl;
#endif
	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);

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
	
	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
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
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		      Gimpl::CovShiftBackward(Umu[nu], nu,
   		        Gimpl::CovShiftIdentityBackward(Utmp, mu)));

	PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;

	// x==
	// |  |
	// y..          // mu polarisation

	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;

	break;
	}
	case 2:{
	  ///////////////// +ve nu /////////////////
	  //     ->
	  //    |  |
	  //    :  |
	  //    x==    
	  // x = y : Computes contr. from this type to the force for U_nu(y) and U_nu(y+nu) <- similar effect happens in all calc. below
	  
	  time=-usecond();
	  PlaqL=Ident;

	  dirs.clear();
	  dirs = {nu+1,nu+1,mu+1,-(nu+1),-(nu+1),-(mu+1) };
	  linkTracer(Umu, Utmp, dirs, 6, -rho, PlaqR); //Umu_tmp
	  time+=usecond();
	  std::cout << GridLogMessage << "PlaqLR took "<<time<< " us Rect"<<norm2(PlaqR)<<" mu "<<mu<<" nu "<<nu<<std::endl;

	  time=-usecond();
	  dJdXe_nMpInv_y =   dJdXe_nMpInv;
	  ComputeNxy(PlaqL,PlaqR,Nxy);
	  Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	  time+=usecond();
	  std::cout << GridLogMessage << "ComputeNxy (occurs 10x) took "<<time<< " us"<<std::endl;

	  time=-usecond();
	  //PlaqR=(-1.0)*PlaqR; //commented out
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx,FdetV);//PlaqR,PlaqL,MpInvJx,FdetV);//PlaqL,PlaqR,MpInvJx,FdetV);
	  Fdet2_nu = FdetV;
	  time+=usecond();
	  std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 10x) took "<<time<< " us"<<std::endl;

	  //     <-
	  //    |  |
	  //    |  :
	  //    x==y    
	  // x = y - mu 

	  dirs.clear();
	  dirs = {nu+1,nu+1,-(mu+1),-(nu+1),-(nu+1)};
	  //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	  PlaqR = rho * Gimpl::CovShiftForward(Umu[nu],nu,
					       Gimpl::CovShiftForward(Umu[nu],nu,
								      Gimpl::CovShiftBackward(Umu[mu],mu,
											      Gimpl::CovShiftBackward(Umu[nu],nu,
														      Gimpl::CovShiftIdentityBackward(Umu[nu],nu)))));
	  PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu); // Note: adj(PlaqL) is used in ComputeNx & Compute_MpInvJx_dNxxdSy

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

	  dirs.clear();
	  PlaqL = Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	  dirs = {nu+1,(mu+1),-(nu+1),-(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 5,-rho, PlaqR);
	  PlaqR = (-rho) * Gimpl::CovShiftForward(Umu[nu],nu,
						  Gimpl::CovShiftForward(Umu[mu],mu,
									 Gimpl::CovShiftBackward(Umu[nu],nu,
												 Gimpl::CovShiftBackward(Umu[nu],nu,
															 Gimpl::CovShiftIdentityBackward(Utmp,mu)))));
												 
												 
	  
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);
	  ComputeNxy(PlaqL, PlaqR,Nxy);
	  Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	  
	  MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);//PlaqL,PlaqR,MpInvJx_nu,FdetV);
	  Fdet2_nu = Fdet2_nu+FdetV;

	
	  //     <-:
	  //    |  y
	  //    |  |
	  //    x==    
	  // x = y - mu - nu

	  dirs.clear();
	  dirs = {-(nu+1),-(mu+1)};
	  //linkTracer(Umu,Utmp, dirs, 2, 1.0,PlaqL);
	  PlaqL = Gimpl::CovShiftBackward(Umu[nu],nu,
					   Gimpl::CovShiftIdentityBackward(Utmp,mu));
	  dirs.clear();
	  dirs = {nu+1,-(mu+1),-(nu+1),-(nu+1)};
	  //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	  PlaqR = rho * Gimpl::CovShiftForward(Umu[nu],nu,
					       Gimpl::CovShiftBackward(Umu[mu],mu,
									Gimpl::CovShiftBackward(Umu[nu],nu,
												 Gimpl::CovShiftIdentityBackward(Umu[nu],nu))));
					       


	  //dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,-1);
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,-1);
	  ComputeNxy(PlaqL, PlaqR,Nxy);
	  Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	  
	  //MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,-1);
	  MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	  MpInvJx_nu = Cshift(MpInvJx_nu,nu,-1);
	  Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	  Fdet2_nu = Fdet2_nu+FdetV;

	
	  ///////////////// -ve nu /////////////////

	  //    x==
	  //    :  |
	  //    y
	  //    |  |
	  //     <-
	  // x = y + nu: Computes contr. from this type to the force for U_nu(y) and U_nu(y+nu) <- similar effect happens in all calc. below

	  dirs.clear();
	  dirs = {-(nu+1),(mu+1),(nu+1),(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 5, rho, PlaqL);
	  PlaqL = rho* Gimpl::CovShiftBackward(Umu[nu],nu,
						 Gimpl::CovShiftForward(Umu[mu],mu,
									Gimpl::CovShiftForward(Umu[nu],nu,
											       Gimpl::CovShiftForward(Umu[nu],nu,
														      Gimpl::CovShiftIdentityBackward(Utmp,mu)))));
														      
	  PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);

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

	  dirs.clear();
	  dirs = {-(nu+1),-(mu+1),(nu+1),(nu+1)};
	  //linkTracer(Umu, Utmp, dirs, -1, -rho, PlaqL);
	  PlaqL = (-rho)*Gimpl::CovShiftBackward(Umu[nu],nu,
						  Gimpl::CovShiftBackward(Umu[mu],mu,
									   Gimpl::CovShiftForward(Umu[nu],nu,
												  Gimpl::CovShiftIdentityForward(Umu[nu], nu))));
	  dirs.clear();
	  dirs = {(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 2, 1.0,PlaqR);
	  PlaqR = Gimpl::CovShiftForward(Umu[nu],nu,
					 Gimpl::CovShiftIdentityBackward(Utmp,mu));
	  
	  //dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,1);
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
          dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);
	  ComputeNxy(PlaqL, PlaqR,Nxy);
	  Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;

	  //MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,1);
	  MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	  MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);//PlaqL,PlaqR,MpInvJx_nu,FdetV);
	  Fdet2_nu = Fdet2_nu+FdetV;
	  

	  //    x==
	  //    |  |
	  //    :  |
	  //    y->
	  // x = y + 2nu
	  dirs.clear();
	  dirs = {nu+1,nu+1};
	  //linkTracer(Umu, Utmp, dirs, -1, rho, PlaqR);
	  PlaqR = rho*Gimpl::CovShiftForward(Umu[nu],nu,
					     Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	  dirs.clear();
	  dirs = {mu+1,-(nu+1),-(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 4, 1.0, PlaqL);
	  PlaqL = Gimpl::CovShiftForward(Umu[mu],mu,
					 Gimpl::CovShiftBackward(Umu[nu],nu,
								  Gimpl::CovShiftBackward(Umu[nu],nu,
											   Gimpl::CovShiftIdentityBackward(Utmp,mu))));
	
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

	  dirs.clear();
	  dirs = {nu+1,nu+1,-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 3, -rho, PlaqR);
	  PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[nu],nu,
						Gimpl::CovShiftForward(Umu[nu],nu,
								       Gimpl::CovShiftIdentityBackward(Utmp,mu)));
	  dirs.clear();
	  dirs = {-(mu+1),nu+1,nu+1};
          //linkTracer(Umu, Utmp, dirs, -1, 1.0, PlaqL);
	  PlaqL = Gimpl::CovShiftBackward(Umu[mu],mu,
					   Gimpl::CovShiftForward(Umu[nu],nu,
								  Gimpl::CovShiftIdentityForward(Umu[nu], nu)));

	  //dJdXe_nMpInv_y = Cshift(Cshift(dJdXe_nMpInv,mu,-1),nu,2);
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
          dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,2);
	  ComputeNxy(PlaqL, PlaqR,Nxy);
	  Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	  
	  //MpInvJx_nu = Cshift(Cshift(MpInvJx,mu,-1),nu,2);
	  MpInvJx_nu = Cshift(MpInvJx,mu,-1);
          MpInvJx_nu = Cshift(MpInvJx_nu,nu,2);
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);//PlaqL,PlaqR,MpInvJx_nu,FdetV);
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
	  dirs.clear();
	  dirs = {(mu+1),-(nu+1),-(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 4, -rho,PlaqR);
	  PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[mu],mu,
						Gimpl::CovShiftBackward(Umu[nu],nu,
									Gimpl::CovShiftBackward(Umu[nu],nu,
												Gimpl::CovShiftIdentityBackward(Utmp,mu))));
	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-2);
	  
	  ComputeNxy(PlaqL,PlaqR,Nxy);
	  Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	  
	  MpInvJx_nu = Cshift(MpInvJx,nu,-2);
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);//PlaqL,PlaqR,MpInvJx_nu,FdetV);
	  Fdet2_mu = Fdet2_mu+FdetV;
	  
	
	  //    x<=
	  //    |  |
	  //    |  |
	  //    y..
	  // x = y + 2nu

	  dirs.clear();
	  PlaqL=Gimpl::CovShiftForward(Umu[nu],nu,Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	  dirs = {(mu+1),(nu+1),(nu+1),-(mu+1)};
	  //linkTracer(Umu, Utmp, dirs, 4, -rho,PlaqR);
	  PlaqR = (-rho)*Gimpl::CovShiftForward(Umu[mu],mu,
						Gimpl::CovShiftForward(Umu[nu],nu,
								       Gimpl::CovShiftForward(Umu[nu],nu,
											      Gimpl::CovShiftIdentityBackward(Utmp,mu))));

	  dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,2);

	  ComputeNxy(PlaqL,PlaqR,Nxy);
	  Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	  
	  MpInvJx_nu = Cshift(MpInvJx,nu,2);
	  Compute_MpInvJx_dNxxdSy(PlaqR,PlaqL,MpInvJx_nu,FdetV);//PlaqL,PlaqR,MpInvJx_nu,FdetV);
	  Fdet2_mu = Fdet2_mu+FdetV;
	  break;
	}
	default:
	  {
	    assert(1!=1 && " not entered force calculation");
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
  
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
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
  RealD logDetJacobian(void)
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
  void logDetJacobianForce(GaugeField &force)
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
      std::cout << GridLogDebug << this->smearingLevels <<" "<<this->SmearedSet.size()<<""<<Nsmr_one_step<<std::endl;//DEBUG
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
        //7 0 1 7 0.1 120790.48252695
    std::cout << GridLogMessage <<"AnalyticSmearedForce REct "<<level<<" "<<smr_ind<<" "<<mask_type<<" "<<sub_level<<" "<<rho<<" "<<norm2(SigmaKPrimeB)<<" "<<norm2(SigmaKPrimeA)<<" "<<norm2(Cmu)<<" "<<norm2(GaugeK)<<std::endl;//DEBUG

    {
      GaugeKmu = peekLorentz(GaugeK, mu); std::cout << GridLogMessage <<"AnalyticSmearedForce REct gaugek "<<norm2(GaugeKmu)<<std::endl;//DEBUG   
      SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);  std::cout << GridLogMessage <<"AnalyticSmearedForce REct SigmaKPrime_mu "<<norm2(SigmaKPrime_mu)<<std::endl;//DEBUG       
      iQ = Ta(Cmu * adj(GaugeKmu));  std::cout << GridLogMessage <<"AnalyticSmearedForce REct iQ "<<norm2(iQ)<<std::endl;//DEBUG                                   
      this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);  std::cout << GridLogMessage <<"AnalyticSmearedForce REct iLambda_mu "<<norm2(iLambda_mu)<<std::endl;//DEBUG
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu); std::cout << GridLogMessage <<"AnalyticSmearedForce REct SigmaK "<<norm2(SigmaK)<<std::endl;//DEBUG 
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
    : SmearedConfigurationMasked<Gimpl>(_UGrid, Nsmear,*Stouts[0]), Stouts(Stouts), mask_types(mask_types)
  {
    assert(Nsmear%(Nsmr_one_step)==0); 
    assert(Nsmear/Nsmr_one_step==mask_types.size()); // Nsmr_one_step = #Basic_num_steps = 2*Nd = 8; One step <-> a mask type
    assert(Stouts.size() == mask_types.size());
    
    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);

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
	case 1:{
	  LatticeComplex tmpcb(UrbGrid);

	  pickCheckerboard(cb,tmpcb,ones);
	  setCheckerboard(tmp,tmpcb);
	  break;
	}
	case 2:{
	  std::cout << GridLogMessage <<" mask type 2 in seting up Mask"<<i<<" "<<j<<" "<<mu<<" "<<cb<<std::endl;
	  //Lattice<iScalar<vInteger> > coor_nu(_UGrid),coor_sum(_UGrid); coor_sum = Zero();
	  LatticeInteger  coor_nu(_UGrid),coor_sum(_UGrid), int_cb(_UGrid); coor_sum = Zero(); int_cb = LatticeInteger::scalar_type(cb);//Integer(cb);
	  // mu direction is made trivial to reduce coding in LogDetJacobianForceLevel routine
	  //for(int nu=0, c=0; nu < Nd; nu++)
	  for(int nu=0; nu < Nd; nu++)
	    if( nu != mu){
	      LatticeCoordinate(coor_nu,nu);//xs[c],nu);
	      coor_nu = div(coor_nu,mask_type);
	      coor_sum = coor_sum + coor_nu;
	      //xs[c] = div(xs[c],mask_type);
	      //c++;
	    }
	  coor_sum = coor_sum +int_cb;
	  coor_sum = mod(coor_sum,2);
	  //Lattice<iScalar<vInteger>> pred(_UGrid); pred = Zero();
	  //for(int nu=0; nu<3; nu++) pred = pred + xs[nu];
	  
	  //tmp = where( mod(coor_sum,2)==(Integer)(cb), ones, zeros); std::cout << GridLogMessage <<"norm mask "<<norm2(tmp)<<std::endl;
	  tmp = where( coor_sum, ones, zeros);

	  LatticeComplex tmp2(_UGrid), tmp3(_UGrid);
	  tmp2 = Cshift(tmp,(mu+1+Nd)%Nd,2);
	  tmp3 = Cshift(tmp,(mu+1+Nd)%Nd,1);
	  Coordinate point1(Nd), point2(Nd), point3(Nd);
	  std::cout << GridLogMessage <<"norm mask "<<norm2(tmp)<<" "<<norm2(ones)<<" "<<norm2(tmp+tmp2)<<" "<<norm2(tmp+tmp3)<<std::endl;
	  for(int nu=0;nu<Nd;nu++) {
	    point1[nu] = 2;
	    point3[nu] = 0;
	  }
	  
	  point2 = point1; point2[mu] = 4;
	  Complex c;
	  peekSite(c,tmp,point1);
	  std::cout << GridLogMessage <<"mask value:  mu="<<mu<<"cb= "<<cb<<" "<<point1<<" "<<c;
	  peekSite(c,tmp,point2); std::cout << point2<< " "<<c;
	  peekSite(c,tmp,point3); std::cout << point3<<" "<<c<<std::endl;
	  
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

