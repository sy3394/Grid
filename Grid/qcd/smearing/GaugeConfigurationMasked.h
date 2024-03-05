/*!
  @file GaugeConfiguration.h
  @brief Declares the GaugeConfiguration class
  
  Convention:
    - follow Luscher's normalization convention for su(3) internally
    - Assume: T^a = -i*t^a where T^a:Luscher; t^a:Grid
*/
#pragma once

NAMESPACE_BEGIN(Grid);

/*!
  @brief Smeared configuration masked container
  Modified for a multi-subset smearing (aka Luscher Flowed HMC)
*/
template <class Gimpl>
class SmearedConfigurationMasked : public SmearedConfiguration<Gimpl>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
  /*
#define INHERIT_GIMPL_TYPES(GImpl)                  \
  typedef typename GImpl::Simd Simd;                \
  typedef typename GImpl::Scalar Scalar;            \
  typedef typename GImpl::LinkField GaugeLinkField; \
  typedef typename GImpl::Field GaugeField;         \
  typedef typename GImpl::ComplexField ComplexField;\
  typedef typename GImpl::SiteField SiteGaugeField; \
  typedef typename GImpl::SiteComplex SiteComplex;  \
  typedef typename GImpl::SiteLink SiteGaugeLink;

  typedef Lattice<SiteLink>    LinkField;
  typedef iImplGaugeLink<Simd>  SiteLink;
  template <typename vtype> using iImplGaugeLink  = iScalar<iScalar<iMatrix<vtype, Nrepresentation> > >; <- a[1][1][Nrep x Nrep][len(vtype)]

  => i) LinkField has no Lorentz index, just a field of SU(3) matrix in some rep; U(x)
     ii) GaugeField has Lorentz index; U_\mu(x)

  /Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/qcd/action/gauge/GaugeImplTypes.h  
   */
private:
  // These live in base class
  //  const unsigned int smearingLevels;
  //  Smear_Stout<Gimpl> *StoutSmearing;
  //  std::vector<GaugeField> SmearedSet;

  /*
    template<typename vtype> using iLorentzComplex            = iVector<iScalar<iScalar<vtype> >, Nd > ; <- arr[len(vtype)][1][1][4]
    typedef iLorentzComplex<vComplex > vLorentzComplex; <- vtype = vComplex; = internal dof
    typedef Lattice<vLorentzComplex>  LatticeLorentzComplex;

    ~/BNL/src/Grid_sy3394/Grid/qcd/QCD.h
   */
  std::vector<LatticeLorentzComplex> masks;

  typedef typename SU3Adjoint::AMatrix AdjMatrix;
  typedef typename SU3Adjoint::LatticeAdjMatrix  AdjMatrixField;
  typedef typename SU3Adjoint::LatticeAdjVector  AdjVectorField;

  // Adjoint vector in Luscher' converntion; If we use T^a_Luscher for expansion, the result su(3) field should be basis indep.
  void InsertForce(GaugeField &Fdet,AdjVectorField &Fdet_nu,int nu)
  {
    Complex ci(0,1);
    GaugeLinkField Fdet_pol(Fdet.Grid());
    Fdet_pol=Zero();
    for(int e=0;e<8;e++){
      ColourMatrix te;
      SU3::generator(e, te);
      auto tmp=peekColour(Fdet_nu,e); 
      Fdet_pol=Fdet_pol + (-1.0)*ci*tmp*te;
    }
    pokeLorentz(Fdet, Fdet_pol, nu);
  }
  // 2nd term of Eq. (53) from Luchang's note: M'^{-1} J(X) d N_xx / ds(y) -> Fdet2
  //  N^{bc} = [ <P(T^c PlaqL^\dagger T^a PlaqR), T^b> = d N_bc /ds^a ] pre-multiplied by M'^{-1} J(X)
  void Compute_MpInvJx_dNxxdSy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR, AdjMatrixField MpInvJx,AdjVectorField &Fdet2 )
  {
    GaugeLinkField UtaU(PlaqL.Grid()); // <- N; contains T^a
    GaugeLinkField D(PlaqL.Grid());    // <- d N/ds; contains T^a as well as T^c from s(y,\nu)^c
    AdjMatrixField Dbc(PlaqL.Grid());
    LatticeComplex tmp(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta,tb,tc;
    
    for(int a=0;a<Ngen;a++) {
      SU3::generator(a, ta);
      UtaU= (-1.0)*ci*adj(PlaqL)*ta*PlaqR; //adj(PlaqL) probably b/c P(U^\dagger C_\mu) is used for projection to su(3)
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, tc);
	// /Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/qcd/utils/GaugeFix.h confirms Ta(U) is proj of U to tracelss anti-hermitian 
	D = Ta( (-1.0)*ci*tc *UtaU); // tc is part of def of N, c.f. Eq. (17, 18); c <-> b there 
	for(int b=0;b<Ngen;b++){
	  SU3::generator(b, tb);
	  tmp =(2.0)*trace(ci*tb*D); 
	  PokeIndex<ColourIndex>(Dbc,tmp,b,c);  // d N_{bc}/ ds^a <- Note: D contains t^a from ds^a
	}
      }
      tmp = trace(MpInvJx * Dbc); // c.f., 2nd term of Eq. (53)
      PokeIndex<ColourIndex>(Fdet2,tmp,a);
    }
  }

  /*
    Note:
      - derivative inserts T^a in plaquettes for Wilson flow
    @param PlaqL: left portion of plaquette to the inserted T^a
    @param PlaqR: right portion of plaquette to the inserted T^a
    @brief Compute N(x,\mu; y,\nu)^{a,b} in Luchang's note
    Comments:
      - The name 'adj' simply b/c it is an 8x8 matrix in {T^a} basis of transformation on su(3)
      - not to be confused with adj rep of some su(3) element
   */
  void ComputeNxy(const GaugeLinkField &PlaqL,const GaugeLinkField &PlaqR,AdjMatrixField &NxAd)
  {
    GaugeLinkField Nx(PlaqL.Grid());
    const int Ngen = SU3Adjoint::Dimension;
    Complex ci(0,1);
    ColourMatrix   ta;
    ColourMatrix   tb;
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, tb);
      Nx = Ta( (-1.0)*adj(PlaqL)*ci*tb * PlaqR ); 
      for(int a=0;a<Ngen;a++) {
	SU3::generator(a, ta);
	///Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/lattice/Lattice_ET.h; forces to eval the expr
	auto tmp =closure( trace((-1.0)*ci*ta*Nx));
	PokeIndex<ColourIndex>(NxAd,tmp,a,b); 
      }
    }
  }
  // U <- U*masks[smr]
  /*
    Apply the same masking for each mu of U_\mu
    Waring: mask is defined in the following, diff from masks[smr]
      masks[smr] is either 0 or 1 for each (x,\mu)
    Recall
      - F(x,mu;y,nu) is force at (y,nu) updated with Z(x,mu)
      - The complete update consists of repeated updates of a single update with Z(x,mu) for all (x,mu) for a number of smearing steps
      - Updates with diferent (x,mu) can be done concurrently as long as they do not depend on each other
      - e.g., U_\mu(x) depends on links in the staples of the to-be-updated link => links of even sites or odd sites
   */
  void ApplyMask(GaugeField &U,int smr)
  {
    LatticeComplex tmp(U.Grid());
    GaugeLinkField Umu(U.Grid());
    for(int mu=0;mu<Nd;mu++){

      // PeekIndex allocates new memory
      ///Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/lattice/Lattice_peekpoke.h
      Umu=PeekIndex<LorentzIndex>(U,mu);
      tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
      Umu=Umu*tmp;
      PokeIndex<LorentzIndex>(U, Umu, mu);
    }
  }
public:

  void logDetJacobianForceLevel(const GaugeField &U, GaugeField &force ,int smr)
  {
    /*
      Lattice(GridBase *grid,ViewMode mode=AcceleratorWriteDiscard) {
      this->_grid = grid;
      resize(this->_grid->oSites());
      assert((((uint64_t)&this->_odata[0])&0xF) ==0);
      this->checkerboard=0;
      SetViewMode(mode);
  }
     */
    GridBase* grid = U.Grid();
    ColourMatrix   tb;
    ColourMatrix   tc;
    ColourMatrix   ta;
    GaugeField C(grid);
    GaugeField Umsk(grid);
    std::vector<GaugeLinkField> Umu(Nd,grid); // <- compiler cast this to  std::vector<GaugeLinkField> Umu(Nd,GaugeLinkField(grid));
    // https://stackoverflow.com/questions/6142830/how-do-i-initialize-a-stl-vector-of-objects-who-themselves-have-non-trivial-cons
    // (3) of https://en.cppreference.com/w/cpp/container/vector/vector
    // This is because of conversion constructor https://www.geeksforgeeks.org/g-fact-35/
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
    int mu= (smr/2) %Nd;// smr needs to go up to 2*Nd to cover every dir

    ////////////////////////////////////////////////////////////////////////////////
    // Mask the gauge field
    ////////////////////////////////////////////////////////////////////////////////
    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the mask for this level of smearing; std::vector<LatticeLorentzComplex> masks;

    // As Umsk is already created and initialized during the construction, (overloaded) assignment op. is called
    //https://www.geeksforgeeks.org/copy-constructor-vs-assignment-operator-in-c/
    Umsk = U;//deep copy  c.f. assignemt temp /Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/lattice/Lattice_base.h
    ApplyMask(Umsk,smr);
    Utmp = peekLorentz(Umsk,mu);// masked gauge field in \mu dir

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
    this->StoutSmearing->BaseSmear(C, U);
    Cmu = peekLorentz(C, mu);

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta projects SU(3) element to su(3) =>  Z lives in Lie algabra
    // /Users/CoffeeBreak/BNL/src/Grid_sy3394/Grid/qcd/utils/GaugeFix.h confirms Ta(U) is proj of U to tracelss anti-hermitian
    Zx  = Ta(Cmu * adj(Umu[mu])); // Z = -P(M) = P(M^\dagger); M = U C^\dagger; M^\dagger = Cmu*Umu
    time+=usecond();
    std::cout << GridLogMessage << "Z took "<<time<< " us"<<std::endl;

    time=-usecond();
    // Move Z to the Adjoint Rep == make_adjoint_representation
    /*
      Note:
        - ZxAd = ZxAd + cplx * TRb
	- TRb: adj rep in Luscher's normalization x (-1)
	    c.f. utils/SUnAdjoint.h
	=> Z is traceless anti-hermitian
      Goal: Compute ad_Z = [Z, . ] = \sum_c Z^c [T^c, . ]= \sum_c Z^c (T^c)_adj
      Note: In Luscher's convention
        - [(T^c)_adj]^{ab} = -f^{abc} where f^{abc} is what shows up for \lambda^a's
	- T^a = -(i/2)\lambda^a
      Here, c.f. SUn.impl.h, SUnAdjoint.h
        - ta = (1/2) \lambda^a = i T^a, i.e., T^a = -it^a
     */
    ZxAd = Zero();
    for(int b=0;b<8;b++) {
      SU3::generator(b, tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb); 
      TRb=-TRb; // iT_adj is in -T^a basis, c.f., SUnAdjoint.h => iT_adj = -T^{Luscher}_adj \therefore we multp. by -1
      cplx = 2.0*trace(ci*tb*Zx); // my convention 1/2 delta ba <- ci = 0 + 1*i
      /*
	Luscher's Convention:
          <A,B>_ta = -2 tr[AB] 
	=> ZxAd = ZxAd + <A,B>_ta (T_b)_adj^Luscher = ZxAd - 2tr[Zx (-i*tb)](T_b)_adj^Luscher = ZxAd +  2tr[Zx i*tb](T_b)_adj^Luscher
       */
      ZxAd = ZxAd + cplx * TRb; 
    }
    time+=usecond();
    std::cout << GridLogMessage << "ZxAd took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // J(x) = 1 + Sum_k=1..N (-Zad)^k/(k+1)! <- If N finite, this is approx
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
    std::cout << GridLogMessage << "Jx took "<<time<< " us"<<std::endl;

    //////////////////////////////////////
    // dJ(x)/dxe <- x = \epsilon Z here; \rho = \epsilon
    //////////////////////////////////////
    /*
      Note:
        - J(x) = 1 + Sum_k=1..N (-Zad)^k/(k+1)!
	- d X/ dX^e = T^e_adj
     */
    time=-usecond();
    std::vector<AdjMatrixField>  dJdX;    dJdX.resize(8,grid); // for each T^a
    AdjMatrixField tbXn(grid);
    AdjMatrixField sumXtbX(grid);
    AdjMatrixField t2(grid);
    AdjMatrixField dt2(grid);
    AdjMatrixField t3(grid);
    AdjMatrixField dt3(grid);
    AdjMatrixField aunit(grid);
    for(int b=0;b<8;b++){
      aunit = ComplexD(1.0);
      SU3Adjoint::generator(b, TRb); //TRb = - (T^b)_adj^Luscher

      X  = (-1.0)*ZxAd; 
      t2 = X;
      dt2 = TRb; // = -ad_Tb^Luscher; - (d ad_X/ dX_e) = - ad_{T^e} <- (-1) comes from -Z_ad in (-Z_ad)^k in J(X)
      for (int j = 20; j > 1; --j) {// t2<->t3, dt2->coeff*dt2 = dt3, dt3+t3->dt2
	t3 = t2*(1.0 / (j + 1))  + aunit; 
	dt3 = dt2*(1.0 / (j + 1)); // rescaling of dt2 by 1/(j+1)
	t2 = X * t3; // 2\sum_j  X^(j)/ (j +1)!
	dt2 = TRb * t3 + X * dt3; // = -T_adj^b*t3 + X*dt2/(j+1)
      }
      dJdX[b] = 0.5*dt2; // the above expr comes with extra factor of 2
    }
    time+=usecond();
    std::cout << GridLogMessage << "dJx took "<<time<< " us"<<std::endl;
    /////////////////////////////////////////////////////////////////
    // Mask Umu for this link
    /////////////////////////////////////////////////////////////////
    time=-usecond();
    PlaqL = Ident; // = 1^{3x3}<-unit color matrix
    PlaqR = Utmp*adj(Cmu); //Utmp = peekLorentz(Umsk,mu); -proj(PlaqR) onto su(3) = Z
    ComputeNxy(PlaqL,PlaqR,NxxAd); //Computed: P(T^b U C_\mu^\dagger) = N_xx in {T^a} basis
    time+=usecond();
    std::cout << GridLogMessage << "ComputeNxy took "<<time<< " us"<<std::endl;
    
    ////////////////////////////
    // Mab: Mp = M' = 1 - \epsilon J(X)N
    ////////////////////////////
    MpAd = Complex(1.0,0.0);
    MpAd = MpAd - JxAd * NxxAd; 

    /////////////////////////
    // invert the 8x8
    /////////////////////////
    time=-usecond();
    MpAdInv = Inverse(MpAd);
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
    nMpInv= NxxAd *MpAdInv; // = N_xx M'^{-1}

    ///////////////////////////////
    // Case: (y,nu) = (x,mu)
    ///////////////////////////////
    AdjMatrixField MpInvJx(grid);
    AdjMatrixField MpInvJx_nu(grid);
    MpInvJx = (-1.0)*MpAdInv * JxAd;// rho is on the plaq factor; // accounts for minus sign of the 2nd term of Eq. (53)
    /*
      PlaqL = Ident = 1^{3x3}<-unit color matrix
      PlaqR = Utmp*adj(Cmu); //Utmp = peekLorentz(Umsk,mu);
         <- If (x,mu) is not among the set of (x,mu)'s used for parallel updates, the contribution from such Z_(x,mu) is 0
     */
    Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
    Fdet2_mu=FdetV;
    Fdet1_mu=Zero();
    
    for(int e =0 ; e<8 ; e++){
      LatticeComplexD tr(grid);
      tr = trace(dJdX[e] * nMpInv); //c.f., 1st term of Eq. 53. prepare adj vec with index e
      pokeColour(dJdXe_nMpInv,tr,e);
    }
    ///////////////////////////////
    // Mask it off
    ///////////////////////////////
    auto tmp=PeekIndex<LorentzIndex>(masks[smr],mu);
    dJdXe_nMpInv = dJdXe_nMpInv*tmp;
    
    //    dJdXe_nMpInv needs to multiply:
    //       Nxx_mu (site local)                           (1) <- 2nd term of Eq. (53) from this contr. is computed above
    //       Nxy_mu one site forward  in each nu direction (3) 
    //       Nxy_mu one site backward in each nu direction (3)
    //       Nxy_nu 0,0  ; +mu,0; 0,-nu; +mu-nu   [ 3x4 = 12]
    // 19 terms.
    /*
      Recall: In F(x,\mu;y,\nu),
        - (x,\mu) is the index on Z and U to be updated
	- [- \sum_{x,\u} F(x,\mu;y,\nu) ] is the 2nd cont. to the seq. updated force (c.f. Eq. 6.5 of Luscher)
	     <- (y,\nu) is the index of force after update
	 => F(x,\mu;y,\nu) is summed over (x,\mu) of parallel updates (done via masking)
      Note:
       - We fix the updating dir to mu, but several x are considered in parallel
     */
    // s.y.) F(x,\mu;y,\nu) <- for a given (y,\nu), need to sum over (x,\mu) <- ds(y)_\nu picks up U_nu(y) in the masked sum
    //         2nd case: take derivative of U(y,mu) with y 
    //         3rd case: U_\nu(y) we take derivative w/r/t/ via ds(y,\nu) shows up in the backward path in plaq starting at (y-\nu,\mu)
    //                     ->
    //                    |  :  take derivative of plaq starting at such (x,\mu) that U_\mu^\dagger(y) shows up, indicated by :
    //                     <- 
    
    AdjMatrixField Nxy(grid);

    GaugeField Fdet1(grid);
    GaugeField Fdet2(grid);
    GaugeLinkField Fdet_pol(grid); // one polarisation

    RealD t4 = usecond();
    for(int nu=0;nu<Nd;nu++){

      if (nu!=mu) {
	// ds(y,\nu) for fixed mu but different x contributing to force propagation
	
	///////////////// +ve nu /////////////////
	//     __
	//    :  |
	//    x==    // nu polarisation -- clockwise
	// It's actually counterclockwise; turned to clockwise by daggering w/ change in sign for N(x,y)
	
	// Case: x = y; \mu \neq \nu
	time=-usecond();
	PlaqL=Ident;

	// Recall: U_{-\mu}(x+\mu) = U_\mu^\dagger(x) => PlaqR = plaq(\mu,\nu) in the clockwise manner
	//  This is why Utmp is used.  It is 0, and thus PlaqR=0 so that the contribution to the force from such Z(x,\mu)

	// N(x,\mu; x,\nu) = d P( U(x,\mu)...U(x,\nu)^\dagger )/ds(x,\nu)^b = P( U(x,\mu)...[T^b U(x,\nu)]^\dagger ) = -P(T^b  U(x,\nu)...U(x,\mu)^\dagger)
	PlaqR=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu, // U_nu(x) U_\mu(x+nu) U_\nu^\dagger(x-nu+mu+nu) Utmp^\dagger(x-nu+nu)_\mu 
	       Gimpl::CovShiftForward(Umu[mu], mu, // U_\mu(x) U_\nu^\dagger(x-nu+mu) Utmp^\dagger(x-mu-nu+mu)_\mu 
		 Gimpl::CovShiftBackward(Umu[nu], nu, // U_\nu^\dagger(x-\nu) Utmp^\dagger(x-mu-nu)_\mu
		   Gimpl::CovShiftIdentityBackward(Utmp, mu)))); // Utmp^\dagger(x-mu)_\mu
	time+=usecond();
	std::cout << GridLogMessage << "PlaqLR took "<<time<< " us"<<std::endl;

	time=-usecond();
	dJdXe_nMpInv_y = dJdXe_nMpInv; //<- a field of vectors of 8 entries
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = transpose(Nxy)*dJdXe_nMpInv_y; //Note: N^{ea} appears in Eq. (53); we take transpose so that correct indicies are contracted
	//  In the above, probably contraction is also done over x to leave only mu,y,nu dependence
	time+=usecond();
	std::cout << GridLogMessage << "ComputeNxy (occurs 6x) took "<<time<< " us"<<std::endl;

	time=-usecond();
	PlaqL=(-1.0)*adj(PlaqR);
	PlaqR=Ident;
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx,FdetV);
	Fdet2_nu = FdetV;
	time+=usecond();
	std::cout << GridLogMessage << "Compute_MpInvJx_dNxxSy s.y. (occurs 6x) took "<<time<< " us"<<std::endl;
	
	//    ___
	//    |  :
	//    x==y    // nu polarisation -- anticlockwise
	// Case: x = y - mu <- Z(y-mu, mu) has non-zero contribution to F(y,nu)
	//       we shift the entire lattice by -mu to bring y into the 'center'

	// Unu(y) Umu^dagger(y-mu+nu)Unu^\dagger(y-mu)
	PlaqR=(rho)*Gimpl::CovShiftForward(Umu[nu], nu,
		      Gimpl::CovShiftBackward(Umu[mu], mu,
    	 	        Gimpl::CovShiftIdentityBackward(Umu[nu], nu)));
	// Umu(y-mu)
	PlaqL=Gimpl::CovShiftIdentityBackward(Utmp, mu); //<- to be daggered

	/*
	  Note:
	    - x in N(x,y) refers to a site near the current site y
	    - x is to be contracted with dJdXe_nMpInv(x)
	    - So we need to bring the value of dJdXe_nMpInv at different site
	 */
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,mu,-1); // to match x index to be contracted; 
	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_nu = Fdet1_nu+transpose(Nxy)*dJdXe_nMpInv_y;
	

	MpInvJx_nu = Cshift(MpInvJx,mu,-1); // to match x index to be contracted
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	// in the presence of the side effects of overloaded method, += is avoided
	// https://stackoverflow.com/questions/34465848/operator-in-c~
	
	///////////////// -ve nu /////////////////
	// x==
	// :  |
	// y__|          // nu polarisation -- clockwise
	// Case: x = y + nu

	// Umu(y)Unu(y+mu+nu)Umu(y+nu)^dag <- to be daggered
	PlaqL=(rho)* Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
			 Gimpl::CovShiftIdentityBackward(Utmp, mu))); 
	// Unu(y) <- take deriv. w/r/t/ this
        PlaqR = Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);
	ComputeNxy(PlaqL,PlaqR,Nxy); // Nxy contains derivative of plaq, thus comes with minus sign
	Fdet1_nu = Fdet1_nu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);
	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_nu = Fdet2_nu+FdetV;
	
	// x==
	// |  :
	// |__y         // nu polarisation
	// Case: x = y - mu + nu

	// Unu(y)Umu(y-mu+nu)^dag <- to be daggered; (-1) comes from [T^a]^dag b/c deriv. w/r/t/ Unu(y)^dag
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[nu], nu,
 	        Gimpl::CovShiftIdentityBackward(Utmp, mu));
	// Umu(y-mu)^dag Unu(y-mu)
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

        //    y""
	//    |  |
	//    x==    // mu polarisation
	// Case: x = y - mu

	// Umu(y)Unu^\dag(y-nu+mu)Utmp_mu^dag(y-nu) <- to be daggered; (-1) comes from [T^a]^dag <- plaq contains Umu(y)^dag
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		      Gimpl::CovShiftBackward(Umu[nu], nu,
   		        Gimpl::CovShiftIdentityBackward(Utmp, mu)));
	// Unu(y)^dag
	PlaqR=Gimpl::CovShiftIdentityBackward(Umu[nu], nu);
	
	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,-1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,-1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;

	// x==
	// |  |
	// y"""           // mu polarisation
	// Case: x = y + nu

	// Umu(y)Unu(y+mu)Utmp_mu(y-mu)^dag <- to be daggered; (-1) comes from [T^a]^dag <- plaq contains Umu(y)^dag
	PlaqL=(-rho)*Gimpl::CovShiftForward(Umu[mu], mu,
		       Gimpl::CovShiftForward(Umu[nu], nu,
		 	 Gimpl::CovShiftIdentityBackward(Utmp, mu)));
	// Unu(y)
        PlaqR=Gimpl::CovShiftIdentityForward(Umu[nu], nu);

	dJdXe_nMpInv_y = Cshift(dJdXe_nMpInv,nu,1);

	ComputeNxy(PlaqL,PlaqR,Nxy);
	Fdet1_mu = Fdet1_mu + transpose(Nxy)*dJdXe_nMpInv_y;

	MpInvJx_nu = Cshift(MpInvJx,nu,1);

	Compute_MpInvJx_dNxxdSy(PlaqL,PlaqR,MpInvJx_nu,FdetV);
	Fdet2_mu = Fdet2_mu+FdetV;
	
      }
    }
    RealD t5 = usecond();

    Fdet1_mu = Fdet1_mu + transpose(NxxAd)*dJdXe_nMpInv;

    InsertForce(Fdet1,Fdet1_mu,mu);
    InsertForce(Fdet2,Fdet2_mu,mu);

    force=  (-1.0)*(Fdet1 + Fdet2);
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
    LatticeComplex       cplx(grid); 
    AdjVectorField  AlgV(grid); 
    AdjMatrixField  Mab(grid);
    AdjMatrixField  Ncb(grid);
    AdjMatrixField  Jac(grid);
    AdjMatrixField  Zac(grid);
    AdjMatrixField  mZac(grid);
    AdjMatrixField  X(grid);

    int mu= (smr/2) %Nd;

    auto mask=PeekIndex<LorentzIndex>(masks[smr],mu); // the cb mask

    /*
      Goal: Compute ln det M'
      Note:
       - det is basis-indep as long as we use the same convention consistently
     */
    //////////////////////////////////////////////////////////////////
    // Assemble the Nxx matrix
    //////////////////////////////////////////////////////////////////
    // Computes ALL the staples -- could compute one only here
    this->StoutSmearing->BaseSmear(C, U);
    Cmu = peekLorentz(C, mu);
    Umu = peekLorentz(U, mu);
    Complex ci(0,1);
    for(int b=0;b<Ngen;b++) {
      SU3::generator(b, Tb);
      Nb = (-1.0)*Ta( ci*Tb * Umu * adj(Cmu));
      for(int c=0;c<Ngen;c++) {
	SU3::generator(c, Tc);
	auto tmp = (2.0)*trace(ci*Tc*Nb);
	PokeIndex<ColourIndex>(Ncb,tmp,c,b); 
      }
    }      

    //////////////////////////////////////////////////////////////////
    // Assemble Luscher exp diff map J matrix 
    //////////////////////////////////////////////////////////////////
    // Ta so Z lives in Lie algabra; Z = -P(Umu*adj(Cmu)) = P(Cmu*adj(Umu))
    Z  = Ta(Cmu * adj(Umu));

    // Move Z to the Adjoint Rep == make_adjoint_representation
    Zac = Zero();
    for(int b=0;b<8;b++) {
      // Is the mapping of these the same? Same structure constants
      // Might never have been checked.
      SU3::generator(b, Tb);         // Fund group sets traceless hermitian T's
      SU3Adjoint::generator(b,TRb);
      TRb=-TRb;
      cplx = 2.0*trace(ci*Tb*Z);
      Zac = Zac + cplx * TRb; 
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
      std::cout << GridLogMessage << "GaugeConfigurationMasked: logDetJacobian took " << time << " ms" << std::endl;  
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
    GaugeField C(grid), SigmaK(grid), iLambda(grid);
    GaugeField SigmaKPrimeA(grid);
    GaugeField SigmaKPrimeB(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);
    
    this->StoutSmearing->BaseSmear(C, GaugeK);
    SigmaK = Zero();
    iLambda = Zero();

    SigmaKPrimeA = SigmaKPrime;
    ApplyMask(SigmaKPrimeA,level);
    SigmaKPrimeB = SigmaKPrime - SigmaKPrimeA;
    
    // Could get away with computing only one polarisation here
    // int mu= (smr/2) %Nd;
    // SigmaKprime_A has only one component
    for (int mu = 0; mu < Nd; mu++)
    {
      Cmu = peekLorentz(C, mu);
      GaugeKmu = peekLorentz(GaugeK, mu);
      SigmaKPrime_mu = peekLorentz(SigmaKPrimeA, mu);
      iQ = Ta(Cmu * adj(GaugeKmu));
      this->set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
      pokeLorentz(iLambda, iLambda_mu, mu);
    }
    this->StoutSmearing->derivative(SigmaK, iLambda,GaugeK);  // derivative of SmearBase

    ////////////////////////////////////////////////////////////////////////////////////
    // propagate the rest of the force as identity map, just add back
    ////////////////////////////////////////////////////////////////////////////////////
    SigmaK = SigmaK+SigmaKPrimeB;

    return SigmaK;
  }

public:

  /* Standard constructor */
  SmearedConfigurationMasked(GridCartesian* _UGrid, unsigned int Nsmear, Smear_Stout<Gimpl>& Stout)
    : SmearedConfiguration<Gimpl>(_UGrid, Nsmear,Stout)
  {
    assert(Nsmear%(2*Nd)==0); // Or multiply by 8??

    // was resized in base class
    assert(this->SmearedSet.size()==Nsmear);
    
    GridRedBlackCartesian * UrbGrid;
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(_UGrid);
    LatticeComplex one(_UGrid); one = ComplexD(1.0,0.0);
    LatticeComplex tmp(_UGrid);

    for (unsigned int i = 0; i < this->smearingLevels; ++i) {

      masks.push_back(*(new LatticeLorentzComplex(_UGrid)));

      int mu= (i/2) %Nd;
      int cb= (i%2);
      LatticeComplex tmpcb(UrbGrid);
	
      masks[i]=Zero();
      ////////////////////
      // Setup the mask
      ////////////////////
      tmp = Zero();
      pickCheckerboard(cb,tmpcb,one);
      setCheckerboard(tmp,tmpcb);
      PokeIndex<LorentzIndex>(masks[i],tmp, mu);
	
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
      std::cout << GridLogMessage << " GaugeConfigurationMasked: Smeared Force chain rule took " << time << " ms" << std::endl;

    }  // if smearingLevels = 0 do nothing
    SigmaTilde=Gimpl::projectForce(SigmaTilde); // Ta
  }

};

NAMESPACE_END(Grid);

