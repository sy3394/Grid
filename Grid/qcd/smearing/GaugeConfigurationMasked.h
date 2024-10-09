
/*!
  @file GaugeConfiguration.h
  @brief Declares the GaugeConfiguration class
*/
#pragma once

NAMESPACE_BEGIN(Grid);


template<class T> void Dump(const Lattice<T> & lat,
			    std::string s,
			    Coordinate site = Coordinate({0,0,0,0}))
{
  typename T::scalar_object tmp;
  if (lat.Checkerboard() != lat.Grid()->CheckerBoard(site))
    site = Coordinate({0,0,0,1});
  peekSite(tmp,lat,site);
  std::cout << " Dump "<<s<<" "<<tmp<<std::endl;
}


/*!
  @brief Smeared configuration masked container
  Modified for a multi-subset smearing (aka Luscher Flowed HMC)
*/
template <class Gimpl>
class SmearedConfigurationMasked : public SmearedConfiguration<Gimpl>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

private:
  // These live in base class
  //  const unsigned int smearingLevels;
  //  Smear_Stout<Gimpl> *StoutSmearing;
  //  std::vector<GaugeField> SmearedSet;

  GridRedBlackCartesian * UrbGrid; // keep a copy of the redblack grid for life of object
  std::vector<LatticeLorentzComplex> masks;
  std::vector<int> cbs;

  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
  typedef typename SU3::vAlgebraMatrix vAlgebraMatrix;
  
  
  // Assume: lat = full lattice
  template<class T>  void printCheckerboards2norm(T &lat, int cb=-1)
  {
    T lat_0(UrbGrid); lat_0 = Zero();
    T lat_1(UrbGrid); lat_1 = Zero();
    pickCheckerboard(0,lat_0,lat);
    pickCheckerboard(1,lat_1,lat);

    std::string parity = (cb==0)? "Even" : "Odd";
    std::cout << GridLogMessage << " printCheckerboards2norm for " << parity << cb << ": Even Part: " << norm2(lat_0) << " Odd Part: " << norm2(lat_1) << std::endl;
  }

  void BaseSmearDerivative(GaugeField& SigmaTerm,
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

    RealD t = 0;
    t-=usecond();
    
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

    t+=usecond();
    std::cout << GridLogMessage << " BaseSmearDerivative " << t/1e3 << " ms " << std::endl;
    
  }
  
  void BaseSmear(GaugeLinkField& Cup, const GaugeField& U,int mu,RealD rho) {
    GridBase *grid = U.Grid();
    GaugeLinkField tmp_stpl(grid);
    WilsonLoops<Gimpl> WL;
    RealD t = 0;

    t-=usecond();
    Cup = Zero();
    for(int nu=0; nu<Nd; ++nu){
      if (nu != mu) {
	// get the staple in direction mu, nu
	WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dagger
	Cup += adj(tmp_stpl*rho);
      }
    }
    t+=usecond();
    std::cout << GridLogMessage << " BaseSmear " << t/1e3 << " ms " << std::endl;
  }
  void BaseSmear_cb(GaugeLinkField& Cup, const GaugeField& U,int mu,RealD rho) {
    GridBase *grid = U.Grid();
    GridBase *hgrid = Cup.Grid();
    GaugeLinkField tmp_stpl(grid);
    GaugeLinkField tmp_stpl_eo(hgrid);
    WilsonLoops<Gimpl> WL;
    int cb = Cup.Checkerboard();
    RealD t = 0;

    t-=usecond();
    Cup = Zero();
    for(int nu=0; nu<Nd; ++nu){
      if (nu != mu) {
        // get the staple in direction mu, nu
        WL.Staple(tmp_stpl, U, mu, nu);  //nb staple conventions of IroIro and Grid differ by a dagger
	pickCheckerboard(cb,tmp_stpl_eo,tmp_stpl); // ideally, compute tmp_stpl only on the current checkerboard
        Cup += adj(tmp_stpl_eo*rho);
      }
    }
    t+=usecond();
    std::cout << GridLogMessage << " BaseSmear " << t/1e3 << " ms " << std::endl;
  }

  // Adjoint vector to GaugeField force
  void InsertForce(GaugeField &Fdet,AdjVectorField &Fdet_nu,int nu)
  {
    Complex ci(0,1);
    GaugeLinkField Fdet_pol(Fdet.Grid());
    RealD t = 0;

    t-=usecond();
    Fdet_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generator(e, te);
      auto tmp=peekColour(Fdet_nu,e);
      Fdet_pol=Fdet_pol + ci*tmp*te; // but norm of te is different.. why?
    }
    pokeLorentz(Fdet, Fdet_pol, nu);

    t+=usecond();
    std::cout << GridLogMessage << " InsertForce " << t/1e3 << " ms " << std::endl;
  }
  void Compute_MpInvJx_dNxxdSy(int cb,
			       const GaugeLinkField &PlaqL,
                               const GaugeLinkField &PlaqR,
                               AdjMatrixField MpInvJx,
                               AdjVectorField &Fdet2 )
  {
    RealD time = -usecond();
    Fdet2 = Zero();
    GaugeLinkField PlaqLeo(UrbGrid);
    GaugeLinkField PlaqReo(UrbGrid);
    AdjMatrixField MpInvJxeo(UrbGrid);
    AdjVectorField Fdet2eo(UrbGrid);
    pickCheckerboard(cb,PlaqLeo,PlaqL);
    pickCheckerboard(cb,PlaqReo,PlaqR);
    pickCheckerboard(cb,MpInvJxeo,MpInvJx);
    Fdet2eo.Checkerboard()=cb;
    time+=usecond();
    Compute_MpInvJx_dNxxdSy(PlaqLeo,PlaqReo,MpInvJxeo,Fdet2eo);
    time-=usecond();
    setCheckerboard(Fdet2,Fdet2eo);
    time+=usecond();
    std::cout << GridLogMessage << " Checkerboarding_MpInvJx_dNxxdSy " << time/1e3 << " ms " << std::endl;
  }
  void Compute_MpInvJx_dNxxdSy_fused(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy_fused");
    int cb = PlaqL.Checkerboard();
    GridBase *grid = PlaqL.Grid();
    const int Ngen = SU3Adjoint::Dimension;
    //AdjMatrix Dbc;
    Complex ci(0,1);
    ColourMatrix ta;
    RealD t;

    t=-usecond();
    autoView(Fdet2_v,Fdet2,AcceleratorWrite);
    autoView(PlaqL_v,PlaqL,AcceleratorRead);
    autoView(PlaqR_v,PlaqR,AcceleratorRead);
    autoView(MpInvJx_v,MpInvJx,AcceleratorRead);
    const int nsimd = vAlgebraMatrix::Nsimd();
    accelerator_for(ss,grid->oSites(),nsimd,{
	typedef decltype(coalescedRead(MpInvJx_v[0])) adj_mat;
	adj_mat Dbc;
	for(int a=0;a<Ngen;a++) {
	  // Qlat Tb = 2i Tb^Grid
	  SU3::generator(a, ta);
	  ta = 2.0 * ci * ta;
	  auto UtaU = adj(PlaqL_v(ss))*ta*PlaqR_v(ss); // 6ms
	  SU3::LieAlgebraProject(Dbc,UtaU);
	  SU3::trace_product(Fdet2_v[ss],MpInvJx_v(ss),Dbc,a);
	}
      });
    t+=usecond();
    std::cout << GridLogMessage << " Compute_MpInvJx_dNxxdSy_fused " << t/1e3 <<" ms"<<std::endl;
  }
  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GRID_TRACE("Compute_MpInvJx_dNxxdSy");
    int cb = PlaqL.Checkerboard();
    GaugeLinkField UtaU(PlaqL.Grid());         UtaU.Checkerboard() = cb;
    GaugeLinkField D(PlaqL.Grid());            D.Checkerboard() = cb;
    AdjMatrixField Dbc(PlaqL.Grid());          Dbc.Checkerboard() = cb;
    AdjMatrixField Dbc_opt(PlaqL.Grid());      Dbc_opt.Checkerboard() = cb;
    LatticeComplex tmp(PlaqL.Grid());          tmp.Checkerboard() = cb;
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    RealD t=0, time;
    RealD tp=0, tpl=0;
    RealD tta=0;
    RealD tpk=0;
    t-=usecond();
    for(int a=0;a<Ngen;a++) {
      tta-=usecond();
      SU3::generator(a, ta);
      ta = 2.0 * ci * ta;
      // Qlat Tb = 2i Tb^Grid
      {
	GRID_TRACE("UtaU");
	UtaU= adj(PlaqL)*ta*PlaqR; // 6ms
      }
      tta+=usecond();
      ////////////////////////////////////////////
      // Could add this entire C-loop to a projection routine
      // for performance. Could also pick checkerboard on UtaU
      // and set checkerboard on result for 2x perf
      ////////////////////////////////////////////
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	tc = 2.0*ci*tc;
	tp-=usecond(); 
	{
	  GRID_TRACE("D = Ta tc*UtaU");
	  D = Ta( tc *UtaU); // 2ms
	}
#if 1
	{
	  GRID_TRACE("LieAlgebraProject");
	  time=-usecond();
	  SU3::LieAlgebraProject(Dbc_opt,D,c); // 5.5ms
	  time+=usecond();
	  tpl+=time;
	  //std::cout << GridLogMessage << " LieAlgebraProject_in_Compute_MpInvJx " << a << " " << c << " "<<time/1e3<<" ms"<<std::endl;
	}
#else
	for(int b=0;b<Ngen;b++){
	  SU3::generator(b, tb);
	  tmp =-trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // Adjoint rep
	}
#endif
	tp+=usecond();
      }
      //Dump(Dbc_opt,"Dbc_opt");
      //Dump(Dbc,"Dbc");
      tpk-=usecond();
      {
	GRID_TRACE("traceMpDbc");
	tmp = trace(MpInvJx * Dbc_opt);
      }
#if 0
      if (a==0){
	GRID_TRACE("traceMpDbc2");
	LatticeComplex tmp2(PlaqL.Grid()); tmp2.Checkerboard() = cb;
	SU3::trace_product(tmp2,MpInvJx,Dbc_opt);
	std::cout << GridLogMessage <<  "DEBUG: Compute_MpInvJx_dNxxdSy_fused trace " << norm2(tmp2-tmp)<<std::endl;
      }
	
#endif	
      {
	GRID_TRACE("pokeIndecx");
	PokeIndex<ColourIndex>(Fdet2,tmp,a);
      }
      tpk+=usecond();
    }
    t+=usecond();
    std::cout << GridLogMessage << " Compute_MpInvJx_dNxxdSy " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
	      << " ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms LieAlgebraProject "<<tpl/1e3<<" ms"<<std::endl;
  }
  
  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GRID_TRACE("ComputeNxy");
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   tb;
    ColourMatrix   tc;
    RealD t = 0, tta = 0, tp = 0, tgen = 0;

    t-=usecond();
    for(int b=0;b<Ngen;b++) {
      tgen-=usecond();
      SU3::generator(b, tb);
      tgen+=usecond();
      tb = 2.0 * ci * tb;
      tta-=usecond();
      {
        GRID_TRACE("UtaU");
	Nx = Ta( adj(PlaqL)*tb * PlaqR );
      }
      tta+=usecond();
      tp-=usecond();
#if 1
      {
	GRID_TRACE("LieAlgebraProject");
	SU3::LieAlgebraProject(NxAd,Nx,b);
      }
#else
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	auto tmp =closure( -trace(ci*tc*Nx));
	{
	  GRID_TRACE("pokeIndecx");
	  PokeIndex<ColourIndex>(NxAd,tmp,c,b);
	}
      }
#endif
      tp+=usecond();
    }
    t+=usecond();
    std::cout << GridLogMessage << " ComputeNxy " << t/1e3 << " ms  proj "<<tp/1e3<< " ms"
              << " ta "<<tta/1e3<<" ms tgen "<< tgen/1e3 << std::endl;
  }
  
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
public:

  void logDetJacobianForceLevel(const GaugeField &U, GaugeField &force ,int smr)
  {
    GRID_TRACE("logDetJacobianForceLevel");
    GridBase* grid = U.Grid();
    GridBase* hgrid = UrbGrid; // For now, assume masking is based on red-black checkerboarding
    ColourMatrix   tb;
    ColourMatrix   tc;
    ColourMatrix   ta;
    GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid);
    GaugeLinkField Cmu(hgrid); // U and staple; C contains factor of epsilon
    GaugeLinkField Zx(hgrid);  // U times Staple, contains factor of epsilon
    GaugeLinkField Utmp(grid);
    GaugeLinkField Ueo(hgrid);
    GaugeLinkField PlaqL(hgrid);
    GaugeLinkField PlaqR(hgrid);
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    ColourMatrix Ident;
    LatticeComplex  cplx(hgrid); 
    
    AdjVectorField  dJdXe_nMpInv(hgrid); 
    AdjVectorField  dJdXe_nMpInv_y(hgrid); 
    AdjMatrixField  MpAd(hgrid);    // Mprime luchang's notes
    AdjMatrixField  MpAdInv(hgrid); // Mprime inverse
    AdjMatrixField  NxxAd(hgrid);    // Nxx in adjoint space
    AdjMatrixField  JxAd(hgrid);     
    AdjMatrixField  ZxAd(hgrid);
    AdjMatrixField  mZxAd(hgrid);
    AdjMatrixField  X(hgrid);
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
    int cb = cbs[smr];
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    Umsk = U;
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);
    pickCheckerboard(cb,Ueo,Utmp);

    Cmu.Checkerboard() = cb;
    cplx.Checkerboard() = cb;
    Zx.Checkerboard() = cb;
    Ueo.Checkerboard() = cb;
    PlaqL.Checkerboard() = cb;
    PlaqR.Checkerboard() = cb;
    MpAd.Checkerboard() = cb;
    MpAdInv.Checkerboard() = cb;
    dJdXe_nMpInv.Checkerboard() = cb;
    dJdXe_nMpInv_y.Checkerboard() = cb;
    NxxAd.Checkerboard() = cb;
    JxAd.Checkerboard() = cb;
    ZxAd.Checkerboard() = cb;
    mZxAd.Checkerboard() = cb;
    X.Checkerboard() = cb;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    double rho=this->StoutSmearing->SmearRho[1];
    int idx=0;
    for(int mu=0;mu<4;mu++){
    for(int nu=0;nu<4;nu++){
      if ( mu!=nu) assert(this->StoutSmearing->SmearRho[idx]==rho);
      else         assert(this->StoutSmearing->SmearRho[idx]==0.0);
      idx++;
    }}
    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only and do it here
    RealD time;
    time=-usecond();
    BaseSmear_cb(Cmu, U, mu, rho);
    
    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Ueo));
    time+=usecond();
    std::cout << GridLogMessage << "Z took "<<time<< " us"<<std::endl;
    
    time=-usecond();
    // Move Z to the Adjoint Rep == make_adjoint_representation
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      GRID_TRACE("ZxAd_b");
      // Adj group sets traceless antihermitian T's -- Guido, really????
      SU3::generator(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba
      ZxAd = ZxAd + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }
    time+=usecond();
    std::cout << GridLogMessage << "ZxAd took "<<time<< " us"<<std::endl;
    
    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    time=-usecond();
    X=1.0; 
    JxAd = X;
    mZxAd = (-1.0)*ZxAd; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      GRID_TRACE("JxAd_k");
      X=X*mZxAd;
      kpfac = kpfac /(k+1);
      JxAd = JxAd + X * kpfac;
    }
    time+=usecond();
    std::cout << GridLogMessage << "Jx took "<<time<< " us"<<std::endl;
    
    //////////////////////////////////////
    // dJ(x)/dxe
    //////////////////////////////////////
    time=-usecond();
#if 1
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,hgrid); for(auto &M : dJdX) M.Checkerboard() = cb; 
    std::vector<AdjMatrix> TRb_s; TRb_s.resize(8);              
    AdjMatrixField tbXn(hgrid);                                 tbXn.Checkerboard() = cb;
    AdjMatrixField sumXtbX(hgrid);                              sumXtbX.Checkerboard() = cb;
    AdjMatrixField t2(hgrid);                                   t2.Checkerboard() = cb;
    AdjMatrixField dt2(hgrid);                                  dt2.Checkerboard() = cb;
    AdjMatrixField t3(hgrid);                                   t3.Checkerboard() = cb;
    AdjMatrixField dt3(hgrid);                                  dt3.Checkerboard() = cb;
    AdjMatrixField aunit(hgrid);                                aunit.Checkerboard() = cb;
    {
      GRID_TRACE("dJx");
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
    for(int b=0;b<8;b++){
      dJdX[b] = -dJdX[b];
    }
    }
#else
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);
    for(int b=0;b<8;b++){
      aunit = ComplexD(1.0);
      SU3Adjoint::generator(b, TRb); //dt2

      X  = (-1.0)*ZxAd; 
      t2 = X;
      dt2 = TRb;
      for (int j = 12; j > 1; --j) {
	t3  = t2*(1.0 / (j + 1))  + aunit;
	dt3 = dt2*(1.0 / (j + 1));
	t2 = X * t3;
	dt2 = TRb * t3 + X * dt3;
      }
      dJdX[b] = -dt2; 
    }
#endif  
    time+=usecond();
    std::cout << GridLogMessage << "dJx took "<<time<< " us"<<std::endl;
    
    /////////////////////////////////////////////////////////////////
    // NxxAd
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Ueo*adj(Cmu);
    ComputeNxy(PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "ComputeNxy took "<<time<< " us"<<std::endl;
    
    ////////////////////////////
    // Mab
    ////////////////////////////
    MpAd = Complex(1.0,0.0);std::cout << GridLogMessage <<"after MpAd def"<<std::endl;Coordinate lcoor,lcoor1;int site=0;hgrid->LocalIndexToLocalCoor(site, lcoor);int site1=1;hgrid->LocalIndexToLocalCoor(site1, lcoor1);
    MpAd = MpAd - JxAd * NxxAd; std::cout << GridLogMessage <<"after MpAd apxy "<< MpAd.Checkerboard()<<" " <<cb<< " "<<hgrid->CheckerBoard(lcoor) << " "<<hgrid->CheckerBoard(lcoor1)<<std::endl;

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    time=-usecond();
    //temporary
    AdjMatrixField tmp(grid);setCheckerboard(tmp,MpAd);
    tmp = Inverse(tmp);
    pickCheckerboard(cb,MpAdInv,tmp);
    time+=usecond();
    std::cout << GridLogMessage << "MpAdInv took "<<time<< " us"<<std::endl;
    
    RealD t3a = usecond();
    /////////////////////////////////////////////////////////////////
    // Nxx Mp^-1
    /////////////////////////////////////////////////////////////////
    AdjVectorField  FdetV(hgrid);     FdetV.Checkerboard() = cb;
    AdjVectorField  Fdet1_nu(grid);
    AdjVectorField  Fdet2_nu(grid);
    AdjVectorField  Fdet2_mu(grid);
    AdjVectorField  Fdet1_mu(grid);

    AdjVectorField  Fdet1_nu_eo(hgrid); Fdet1_nu_eo.Checkerboard() = cb;       
    AdjVectorField  Fdet1_nu_oe(hgrid); Fdet1_nu_oe.Checkerboard() = (cb+1)%2; 
    AdjVectorField  Fdet2_nu_eo(hgrid); Fdet2_nu_eo.Checkerboard() = cb;       
    AdjVectorField  Fdet2_nu_oe(hgrid); Fdet1_nu_oe.Checkerboard() = (cb+1)%2;
    
    AdjVectorField  Fdet1_mu_oe(hgrid); Fdet1_mu_oe = Zero(); Fdet1_mu_oe.Checkerboard() = (cb+1)%2;
    AdjVectorField  Fdet2_mu_oe(hgrid); Fdet2_mu_oe = Zero(); Fdet2_mu_oe.Checkerboard() = (cb+1)%2;
    
    AdjMatrixField nMpInv(hgrid);       nMpInv.Checkerboard() = cb;
    nMpInv= NxxAd *MpAdInv;

    AdjMatrixField MpInvJx(hgrid);      MpInvJx.Checkerboard() = cb; 
    AdjMatrixField MpInvJx_nu(hgrid);   MpInvJx_nu.Checkerboard() = cb;
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor

    LatticeComplexD tr(hgrid); tr.Checkerboard() = cb;
    for(int e =0 ; e<8 ; e++){
      //      ColourMatrix te;
      //      SU3::generator(e, te);
      tr = trace(dJdX[e] * nMpInv);
      pokeColour(dJdXe_nMpInv,tr,e);
    }

    setCheckerboard(Fdet1_mu, (AdjVectorField) (transpose(NxxAd)*dJdXe_nMpInv)); 
    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
    setCheckerboard(Fdet2_mu,FdetV);
#if 1
    AdjVectorField  FdetV2(hgrid);     FdetV2.Checkerboard() = cb; //FdetV2=Zero();
    Compute_MpInvJx_dNxxdSy_fused(PlaqL,PlaqR,MpInvJx,FdetV2);
    std::cout << GridLogMessage << " DEBUG: logDetJacobianForce_F_detVdiff " <<smr<<" "<<mu<<" "<<cb<<" "<<" simd "<<AdjMatrix::Nsimd()<<" "<<norm2(FdetV-FdetV2)<<" "<<norm2(FdetV)<<" " <<norm2(FdetV2)<<std::endl;
#endif

    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1)
    //       Nxy_mu one site forward  in each nu direction (3)
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.
    AdjMatrixField Nxy(hgrid);

    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);

    RealD t4 = usecond(), tLR = 0, tNxy = 0, tMJx = 0;
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	GRID_TRACE("MuNuLoopBody");
	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise

	time=-usecond(); tLR -= usecond();
	PlaqL=Ident;
	{
	  GRID_TRACE("Staple");
	pickCheckerboard(cb,PlaqR,(GaugeLinkField) ((-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
				          Gimpl::CovShiftForward(Umu[mu], mu,
				           Gimpl::CovShiftBackward(Umu[nu], nu,
					    Gimpl::CovShiftIdentityBackward(Utmp, mu))))));
	}
	time+=usecond(); tLR += usecond();
	std::cout << GridLogMessage << "PlaqLR took "<<time<< " us"<<std::endl;

	time=-usecond(); tNxy -= usecond();
	PlaqL.Checkerboard() = cb; Nxy.Checkerboard() = cb; FdetV.Checkerboard() = cb;
	
	dJdXe_nMpInv_y = dJdXe_nMpInv;
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu_eo = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond(); tNxy += usecond();
	std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond(); tMJx -= usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu_eo = FdetV;
	time+=usecond(); tMJx += usecond();
	std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	pickCheckerboard((cb+1)%2,PlaqR,(GaugeLinkField) ((rho)*Gimpl::CovShiftForward(Umu[nu], nu,
					       Gimpl::CovShiftBackward(Umu[mu], mu,
								       Gimpl::CovShiftIdentityBackward(Umu[nu], nu)))));

	pickCheckerboard((cb+1)%2,PlaqL, (GaugeLinkField) (Gimpl::CovShiftIdentityBackward(Utmp, mu)));
	}
	tLR += usecond();

	tNxy -= usecond();
	Nxy.Checkerboard() = (cb+1)%2; 	FdetV.Checkerboard() = (cb+1)%2;
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	ComputeNxy(PlaqL, PlaqR,Nxy);
	Fdet1_nu_oe = transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	
	Fdet2_nu_oe = FdetV;
	tMJx += usecond();
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	pickCheckerboard((cb+1)%2,PlaqL,(GaugeLinkField) ((rho)* Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftForward(Umu[nu], nu,
								       Gimpl::CovShiftIdentityBackward(Utmp, mu)))));

        pickCheckerboard((cb+1)%2,PlaqR, (GaugeLinkField) (Gimpl::CovShiftIdentityForward(Umu[nu], nu)));
	}
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu_oe = Fdet1_nu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu_oe = Fdet2_nu_oe+FdetV;
	tMJx += usecond();
	
	// x==
	// |  |
	// |__|         // nu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	pickCheckerboard(cb,PlaqL,(GaugeLinkField) ((-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
										  Gimpl::CovShiftIdentityBackward(Utmp, mu))));

	pickCheckerboard(cb,PlaqR,(GaugeLinkField) (Gimpl::CovShiftBackward(Umu[mu], mu,
									    Gimpl::CovShiftIdentityForward(Umu[nu], nu))));
	}
	tLR += usecond();

	tNxy -= usecond();
	Nxy.Checkerboard() = cb;  FdetV.Checkerboard() = cb;
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu_eo = Fdet1_nu_eo + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu_eo = Fdet2_nu_eo+FdetV;
	tMJx += usecond();
	
	/////////////////////////////////////////////////////////////////////
	// Set up the determinant force contribution in 3x3 algebra basis
	/////////////////////////////////////////////////////////////////////
	setCheckerboard(Fdet1_nu, Fdet1_nu_eo);
	setCheckerboard(Fdet1_nu, Fdet1_nu_oe);
	InsertForce(Fdet1,Fdet1_nu,nu);
	setCheckerboard(Fdet2_nu, Fdet2_nu_eo);
        setCheckerboard(Fdet2_nu, Fdet2_nu_oe);
	InsertForce(Fdet2,Fdet2_nu,nu);
	
	//////////////////////////////////////////////////
	// Parallel direction terms
	//////////////////////////////////////////////////

	Nxy.Checkerboard() = (cb+1)%2; FdetV.Checkerboard() = (cb+1)%2;
	
        //     __
	//    |  "
	//    |__"x    // mu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	pickCheckerboard((cb+1)%2,PlaqL,(GaugeLinkField) ((-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftBackward(Umu[nu], nu,
									Gimpl::CovShiftIdentityBackward(Utmp, mu)))));

	pickCheckerboard((cb+1)%2,PlaqR,(GaugeLinkField) (Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));
	}
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu_oe = Fdet1_mu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,-1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu_oe = Fdet2_mu_oe+FdetV;
	tMJx += usecond();

	//  __
	// "  |
	// x__|          // mu polarisation
	tLR -= usecond();
	{
	  GRID_TRACE("Staple");
	pickCheckerboard((cb+1)%2,PlaqL,(GaugeLinkField) ((-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
						Gimpl::CovShiftForward(Umu[nu], nu,
								       Gimpl::CovShiftIdentityBackward(Utmp, mu)))));

	pickCheckerboard((cb+1)%2,PlaqR,(GaugeLinkField) (Gimpl::CovShiftIdentityForward(Umu[nu], nu)));
	}
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu_oe = Fdet1_mu_oe + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu_oe = Fdet2_mu_oe+FdetV;
	tMJx += usecond();
	
      }
    }
    RealD t5 = usecond();
    setCheckerboard(Fdet1_mu, Fdet1_mu_oe);
    setCheckerboard(Fdet2_mu, Fdet2_mu_oe);
    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force= (-0.5)*( Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t4-t3 dJdXe_nMpInv "<<t4-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us Plaq "
	      <<tLR/1e3<<" ms Nxy "<<tNxy/1e3<<" ms MpInvJx_dNxxdSy "<<tMJx/1e3<<" ms"<<std::endl;
    std::cout << GridLogMessage << " logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl; // turn adj vec to SU3 force
    std::cout << GridLogMessage << " logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }
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

    ////////////////////////////////////////////////////////////////////////////////
    // Retrieve the eps/rho parameter(s) -- could allow all different but not so far
    ////////////////////////////////////////////////////////////////////////////////
    double rho=this->StoutSmearing->SmearRho[1];
    int idx=0;
    for(int mu=0;mu<4;mu++){
    for(int nu=0;nu<4;nu++){
      if ( mu!=nu) assert(this->StoutSmearing->SmearRho[idx]==rho);
      else         assert(this->StoutSmearing->SmearRho[idx]==0.0);
      idx++;
    }}
    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only and do it here
    RealD time;
    time=-usecond();
    BaseSmear(Cmu, U,mu,rho);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra
    Zx  = Ta(Cmu * adj(Umu[mu]));
    time+=usecond();
    std::cout << GridLogMessage << "Full: Z took "<<time<< " us"<<std::endl;

    time=-usecond();
    // Move Z to the Adjoint Rep == make_adjoint_representation
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      // Adj group sets traceless antihermitian T's -- Guido, really????
      SU3::generator(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba
      ZxAd = ZxAd + cplx * TRb; // is this right? YES - Guido used Anti herm Ta's and with bloody wrong sign.
    }
    time+=usecond();
    std::cout << GridLogMessage << "Full: ZxAd took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
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
    std::cout << GridLogMessage << "Full: Jx took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // dJ(x)/dxe
    //////////////////////////////////////
    time=-usecond();
#if 1
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    std::vector<AdjMatrix> TRb_s; TRb_s.resize(8);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);

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
    for(int b=0;b<8;b++){
      dJdX[b] = -dJdX[b];
    }
#else
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid);
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);
    for(int b=0;b<8;b++){
      aunit = ComplexD(1.0);
      SU3Adjoint::generator(b, TRb); //dt2

      X  = (-1.0)*ZxAd; 
      t2 = X;
      dt2 = TRb;
      for (int j = 12; j > 1; --j) {
	t3  = t2*(1.0 / (j + 1))  + aunit;
	dt3 = dt2*(1.0 / (j + 1));
	t2 = X * t3;
	dt2 = TRb * t3 + X * dt3;
      }
      dJdX[b] = -dt2; 
    }
#endif  
    time+=usecond();
    std::cout << GridLogMessage << "Full: dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident;
    PlaqR = Utmp*adj(Cmu);
    ComputeNxy(PlaqL,PlaqR,NxxAd);
    time+=usecond();
    std::cout << GridLogMessage << "Full: ComputeNxy took "<<time<< " us"<<std::endl;
    
    ////////////////////////////
    // Mab
    ////////////////////////////
    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxxAd;

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    time=-usecond();
    MpAdInv = Inverse(MpAd);
    time+=usecond();
    std::cout << GridLogMessage << "Full: MpAdInv took "<<time<< " us"<<std::endl;
    
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

    RealD t4 = usecond(), tLR = 0, tNxy = 0, tMJx = 0;
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	///////////////// +ve nu /////////////////
	//     __
	//    |  |
	//    x==    // nu polarisation -- clockwise

	time=-usecond(); tLR -= usecond();
	PlaqL=Ident;

	PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	       Gimpl::CovShiftForward(Umu[mu], mu,
	         Gimpl::CovShiftBackward(Umu[nu], nu,
		   Gimpl::CovShiftIdentityBackward(Utmp, mu))));
	time+=usecond(); tLR += usecond();
	std::cout << GridLogMessage << "Full: PlaqLR took "<<time<< " us"<<std::endl;

	time=-usecond(); tNxy -= usecond();
	dJdXe_nMpInv_y =   dJdXe_nMpInv;
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y;
	time+=usecond(); tNxy += usecond();
	std::cout << GridLogMessage << "Full: ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond(); tMJx -= usecond();
	PlaqR=(-1.0)*PlaqR;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu = FdetV;
	time+=usecond(); tMJx += usecond();
	std::cout << GridLogMessage << "Full: Compute_MpInvJx_dNxxSy (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    x==
	//    |  |
	//    .__|    // nu polarisation -- anticlockwise

	tLR -= usecond();
	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));

	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	ComputeNxy(PlaqL, PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	///////////////// -ve nu /////////////////
	//  __
	// |  |
	// x==          // nu polarisation -- clockwise

	tLR -= usecond();
	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	// x==
	// |  |
	// |__|         // nu polarisation

	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	        Gimpl::CovShiftIdentityBackward(Utmp, mu));

	PlaqR=Gimpl::CovShiftBackward(Umu[mu], mu,
	        Gimpl::CovShiftIdentityForward(Umu[nu], nu));
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1);
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv_y,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,mu,-1);
	MpInvJx_nu = Cshift(MpInvJx_nu,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	tMJx += usecond();
	
	/////////////////////////////////////////////////////////////////////
	// Set up the determinant force contribution in 3x3 algebra basis
	/////////////////////////////////////////////////////////////////////
	InsertForce(Fdet1,Fdet1_nu,nu);
	InsertForce(Fdet2,Fdet2_nu,nu);
	
	//////////////////////////////////////////////////
	// Parallel direction terms
	//////////////////////////////////////////////////

        //     __
	//    |  "
	//    |__"x    // mu polarisation
	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		      Gimpl::CovShiftBackward(Umu[nu], nu,
   		        Gimpl::CovShiftIdentityBackward(Utmp, mu)));

	PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	tMJx += usecond();

	//  __
	// "  |
	// x__|          // mu polarisation
	tLR -= usecond();
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));

        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);
	tLR += usecond();

	tNxy -= usecond();
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;
	tNxy += usecond();

	tMJx -= usecond();
	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	tMJx += usecond();
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force= (-0.5)*( Fdet1 + Fdet2);
    RealD t1 = usecond();
    std::cout << GridLogMessage << " Full: logDetJacobianForce t3-t0 "<<t3a-t0<<" us "<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t4-t3 dJdXe_nMpInv "<<t4-t3a<<" us "<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t5-t4 mu nu loop "<<t5-t4<<" us Plaq "
	      <<tLR/1e3<<" ms Nxy "<<tNxy/1e3<<" ms MpInvJx_dNxxdSy "<<tMJx/1e3<<" ms"<<std::endl;
    std::cout << GridLogMessage << " Full: logDetJacobianForce t1-t5 "<<t1-t5<<" us "<<std::endl; // turn adj vec to SU3 force
    std::cout << GridLogMessage << " Full: logDetJacobianForce level took "<<t1-t0<<" us "<<std::endl;
  }
  RealD logDetJacobianLevel(const GaugeField &U,int smr)
  {
    GridBase* grid = U.Grid();
    GaugeField C(grid);
    GaugeLinkField Nb(grid), Nb_opt(grid);
    GaugeLinkField Z(grid);
    GaugeLinkField Umu(grid), Cmu(grid);
    ColourMatrix   Tb;
    ColourMatrix   Tc;
    typedef typename SU3Adjoint::AMatrix AdjMatrix;
    typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
    typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;
    const int Ngen = SU3Adjoint::Dimension;
    AdjMatrix TRb;
    LatticeComplex       cplx(grid); 
    AdjVectorField  AlgV(grid); 
    AdjMatrixField  Mab(grid);
    AdjMatrixField  Ncb(grid);
    AdjMatrixField  Jac(grid);
    AdjMatrixField  Zac(grid);
    AdjMatrixField  mZac(grid);
    AdjMatrixField  X(grid);

    RealD time=0, tta=0, tpk=0, tN=0, tZ=0, tJ=0, tlnDetM=0;
    int mu= (smr/2) %Nd;

    time -= usecond();
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    //////////////////////////////////////////////////////////////////
    // Assemble the N matrix
    //////////////////////////////////////////////////////////////////

    tN -= usecond();
    double rho=this->StoutSmearing->SmearRho[1];
    BaseSmear(Cmu, U,mu,rho);

    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, Tb);
      // Qlat Tb = 2i Tb^Grid
      tta -= usecond();
      Nb = (2.0)*Ta( ci*Tb * Umu * adj(Cmu));
      tta += usecond();
      // FIXME -- replace this with LieAlgebraProject
#if 1
      SU3::LieAlgebraProject(Ncb_opt,Nb,b);
#else
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, Tc);
	auto tmp = -trace(ci*Tc*Nb); // Luchang's norm: (2Tc) (2Td) N^db = -2 delta cd N^db // - was important
	tpk -= usecond();
	PokeIndex<ColourIndex>(Ncb,tmp,c,b);
	tpk += usecond();
      }
#endif
    }
    Dump(Ncb_opt,"Ncb_opt");
    Dump(Ncb,"Ncb");
    tN += usecond();
    
    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    tZ -= usecond();
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
    tZ += usecond();
    
    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zac)^k/(k+1)!
    //////////////////////////////////////
    tJ -= usecond();
    X=1.0; 
    Jac = X;
    mZac = (-1.0)*Zac; 
    RealD kpfac = 1;
    for(int k=1;k<12;k++){
      X=X*mZac;
      kpfac = kpfac /(k+1);
      Jac = Jac + X * kpfac;
    }
    tJ += usecond();

    tlnDetM -= usecond();
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
    tlnDetM += usecond();
    time += usecond();
    std::cout << GridLogMessage << " logDetJacobianLevel " << time/1e3 << " ms ta "<<tta/1e3<<" ms" << " poke "<<tpk/1e3<< " ms"
	      <<" N "<<tN/1e3<<" ms Z "<<tZ/1e3<<" ms J " <<tJ/1e3<<" ms lnDetM "<<tlnDetM/1e3<<" ms" <<std::endl;
    
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
      std::cout << GridLogMessage << "GaugeConfigurationMasked: logDetJacobian took " << time << " ms" << std::endl;  
    }
    return ln_det;
  }
  void logDetJacobianForce(int old, GaugeField &force)
  {
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond(), tas = 0;

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

	tas -= usecond();
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);
	tas += usecond();
	
	// Add back U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = peekLorentz(this->get_smeared_conf(ismr - 1), mu) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}
    	
	// Get this levels determinant force
	force_det = Zero();
	logDetJacobianForceLevel(old, this->get_smeared_conf(ismr-1),force_det,ismr);

	// Sum the contributions
	force = force + force_det;
      }
    
      // remove U in UdSdU...
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = adj(peekLorentz(this->get_smeared_conf(0), mu)) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      tas -= usecond();
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);
      tas += usecond();
      
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(old, *this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
    }  // if smearingLevels = 0 do nothing
    std::cout << GridLogMessage << " DEBUG: logDetJacobianForce Full " << std::endl;
  }
  
  void logDetJacobianForce(GaugeField &force)
  {
    force =Zero();
    GaugeField force_det(force.Grid());

    if (this->smearingLevels > 0)
    {
      double start = usecond(), tas = 0;

      GaugeLinkField tmp_mu(force.Grid());

      for (int ismr = this->smearingLevels - 1; ismr > 0; --ismr) {

	// remove U in UdSdU...
	for (int mu = 0; mu < Nd; mu++) {
	  tmp_mu = adj(peekLorentz(this->get_smeared_conf(ismr), mu)) * peekLorentz(force, mu);
	  pokeLorentz(force, tmp_mu, mu);
	}

	tas -= usecond();
      	// Propagate existing force
        force = this->AnalyticSmearedForce(force, this->get_smeared_conf(ismr - 1), ismr);
	tas += usecond();
	
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

      tas -= usecond();
      force = this->AnalyticSmearedForce(force, *this->ThinLinks,0);
      tas += usecond();
      
      for (int mu = 0; mu < Nd; mu++) {
	tmp_mu = peekLorentz(*this->ThinLinks, mu) * peekLorentz(force, mu);
	pokeLorentz(force, tmp_mu, mu);
      }

      force_det = Zero();

      logDetJacobianForceLevel(*this->ThinLinks,force_det,0);

      force = force + force_det;

      force=Ta(force); // Ta
      
#if 1 // debug
      GaugeField force_debug(force.Grid()); 
      logDetJacobianForce(1,force_debug);
      std::cout << GridLogMessage << " DEBUG: logDetJacobianForce_diff " << norm2(force-force_debug) << std::endl;
#endif
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: AnalyticSmearedForce in lnDetJacobianForce took " << tas/1e3 << " ms" << std::endl;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: lnDetJacobianForce took " << time << " ms" << std::endl;  
    }  // if smearingLevels = 0 do nothing
  }

private:
  //====================================================================
  // Override base clas here to mask it
  virtual void fill_smearedSet(GaugeField &U)
  {
    this->ThinLinks = &U;  // attach the smearing routine to the field U

    // check the pointer is not null
    if (this->ThinLinks == NULL)
      std::cout << GridLogError << "[SmearedConfigurationMasked] Error in ThinLinks pointer\n";

    if (this->smearingLevels > 0)
    {
      std::cout << GridLogMessage << "[SmearedConfigurationMasked] Filling SmearedSet\n";
      GaugeField previous_u(this->ThinLinks->Grid());

      GaugeField smeared_A(this->ThinLinks->Grid());
      GaugeField smeared_B(this->ThinLinks->Grid());

      previous_u = *this->ThinLinks;
      double start = usecond();
      for (int smearLvl = 0; smearLvl < this->smearingLevels; ++smearLvl)
      {
        this->StoutSmearing->smear(smeared_A, previous_u);
	ApplyMask(smeared_A,smearLvl);
	smeared_B = previous_u;
	ApplyMask(smeared_B,smearLvl);
	// Replace only the masked portion
	this->SmearedSet[smearLvl] = previous_u-smeared_B + smeared_A;
        previous_u = this->SmearedSet[smearLvl];

        // For debug purposes
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(previous_u);
        std::cout << GridLogMessage << "[SmearedConfigurationMasked] smeared Plaq: " << impl_plaq << std::endl;
      }
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << "GaugeConfigurationMasked: Link smearing took " << time << " ms" << std::endl;  
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

    int mmu= (level/2) %Nd;
    int cb= (level%2);
    double rho=this->StoutSmearing->SmearRho[1];

    RealD time = 0;

    time -= usecond();
    // Can override this to do one direction only.
    SigmaK = Zero();
    iLambda = Zero();

    SigmaKPrimeA = SigmaKPrime;
    ApplyMask(SigmaKPrimeA,level);
    SigmaKPrimeB = SigmaKPrime - SigmaKPrimeA;
    
    // Could get away with computing only one polarisation here
    // int mu= (smr/2) %Nd;
    // SigmaKprime_A has only one component
#if 0
    BaseSmear(Cmu, GaugeK,mu,rho);
    GaugeKmu = peekLorentz(GaugeK, mu);
    SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
    iQ = Ta(Cmu * adj(GaugeKmu));
    this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
    pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
    pokeLorentz(iLambda, iLambda_mu, mu);
    BaseSmearDerivative(SigmaK, iLambda,GaugeK,mu,rho);  // derivative of SmearBase
#else
    //    GaugeField C(grid);
    //    this->StoutSmearing->BaseSmear(C, GaugeK);
    //    for (int mu = 0; mu < Nd; mu++)
    int mu =mmu;
    BaseSmear(Cmu, GaugeK,mu,rho);
    {
      // Cmu = peekLorentz(C, mu);
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
    BaseSmearDerivative(SigmaK, iLambda,GaugeK,mu,rho);  // derivative of SmearBase
    //    this->StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase
    //    SigmaKcopy = SigmaKcopy - SigmaK;
    //    std::cout << " BaseSmearDerivative fast path error" <<norm2(SigmaKcopy)<<std::endl;
#endif
    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;
    time += usecond();

    std::cout << GridLogMessage << "GaugeConfigurationMasked: Analytic smearing force  took " << time/1e3 << " ms" << std::endl;
    return SigmaK;
  }

public:

  /* Standard constructor */
  virtual ~SmearedConfigurationMasked()
  {
    delete UrbGrid;
  }
  SmearedConfigurationMasked(GridCartesian* _UGrid, unsigned int Nsmear, Smear_Stout<Gimpl>& Stout)
    : SmearedConfiguration<Gimpl>(_UGrid, Nsmear,Stout)
  {
    assert(Nsmear%(2*Nd)==0); // Or multiply by 8??

    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);
    
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    LatticeComplex one(_UGrid); one = ComplexD(1.0,0.0);
    LatticeComplex tmp(_UGrid);

    for (unsigned int i = 0; i < this->smearingLevels; ++i) {

      masks.push_back(*(new LatticeLorentzComplex(_UGrid)));

      int mu= (i/2) %Nd;
      int cb= (i%2);
      LatticeComplex tmpcb(UrbGrid);

      cbs.push_back(cb);
      
      masks[i]=Zero();
      ////////////////////
      // Setup the mask
      ////////////////////
      tmp = Zero();
      pickCheckerboard(cb,tmpcb,one);
      setCheckerboard(tmp,tmpcb);
      PokeIndex<LorentzIndex>(masks[i],tmp, mu);
	
    }
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
      std::cout << GridLogMessage << " GaugeConfigurationMasked: Smeared Force chain rule took " << time << " ms" << std::endl;

    }  // if smearingLevels = 0 do nothing
    SigmaTilde=Gimpl::projectForce(SigmaTilde); // Ta
  }

};

NAMESPACE_END(Grid);

