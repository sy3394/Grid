/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/smearing/StoutSmearing.h
 
 Copyright (C) 2019
 
 Author: unknown
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
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
/*
  @file RectSmearing.h
  @brief Declares Stout smearing class with rectangular kernel
*/
#pragma once

NAMESPACE_BEGIN(Grid);

/*!  @brief Stout smearing of link variable using (1x1, 1x2, 2x1 Wilson loops)
     @details kernel: plaquette, 2x1 rectangle, 1x2 rectangle <- mu-dir in the x-dir
              Rs: 2x1 rectangle: U_\mu is the shorter side of the rectangle
	      Rl: 1x2 rectangle: U_\mu is part of the longer side of the rectangle

!!! Fow now, smearing is done in the order plq -> Rs -> Rl
     However, the order matters, and this complicate force calculation.  Also, it is diff 
     from the one with its kernel sum of the three, i.e., C in Morningstar and Peardon
      => the user should specify which kernel (Rs, Rl) at the time of construction
  
*/
template <class Gimpl>
class Rect_Stout : public Smear_Stout<Gimpl> {
protected:
  int OrthogDimRs = -1;
  int OrthogDimRl = -1; // TODO: implement this vertical rectangle smearing
  
public:
  INHERIT_GIMPL_TYPES(Gimpl)
    
  const std::vector<double> SmearRhoRs, SmearRhoRl;

  // to be used in analytic smearing and derivative
  // mu: dir along which force is computed, i.e., dir of updating link; nu: dir of link smeared in a given step
  static void RectStapleUnoptimisedRsUpper(GaugeLinkField &Stap, const GaugeField &Umu,
					   int mu, int nu) {
    GridBase *grid = Umu.Grid();

    std::vector<GaugeLinkField> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    Stap = Zero();

    if (nu != mu) {

      //  ^ mu  
      //  |     --> nu
      //      ->-
      //      |  |
      //         |
      //      x-<-
      // shift by mu Umu(x)Unu(x+mu)Udag_mu(x-mu+nu+mu)Udag_mu(x-mu-mu+nu+mu)Udag_nu(x-nu-2mu+nu+mu)
      Stap += Gimpl::ShiftStaple(Gimpl::CovShiftForward( U[mu], mu,
				       Gimpl::CovShiftForward( U[nu], nu,
							       Gimpl::CovShiftBackward( U[mu], mu,
										       Gimpl::CovShiftBackward( U[mu], mu,
													       Gimpl::CovShiftIdentityBackward(U[nu], nu)))))
				 ,mu);

    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Optimised short-side rectangle staple on a padded (ghost) grid.
  //
  // RectStapleStencilRs builds, for direction mu, the stencil of the ten
  // link reads per nu!=mu consumed by RectStaplePaddedRs; the two must be
  // kept in step. The stencil lives on the padded grid, which needs
  // depth>=2 because the paths reach two hops in nu.
  //
  // RectStaplePaddedRs writes Cup = rho * adj(Stap) on the padded grid,
  // where Stap is the staple of WilsonLoops::RectStapleUnoptimisedRs;
  // the caller extracts the interior (PaddedCell::Extract).
  //////////////////////////////////////////////////////////////////////////
  static GeneralLocalStencil RectStapleStencilRs(GridBase *ggrid, int mu) {
    std::vector<Coordinate> shifts;
    for (int nu = 0; nu < Nd; nu++) {
      if (nu == mu) continue;
      auto shift = [mu, nu](int smu, int snu) {
	Coordinate s(Nd, 0); s[mu] = smu; s[nu] = snu; return s;
      };
      // upper: U_nu(x) U_nu(x+nu) U_mu(x+2nu) U_nu^dag(x+mu+nu) U_nu^dag(x+mu)
      shifts.push_back(shift(0, 0));
      shifts.push_back(shift(0, 1));
      shifts.push_back(shift(0, 2));
      shifts.push_back(shift(1, 1));
      shifts.push_back(shift(1, 0));
      // lower: U_nu^dag(x-nu) U_nu^dag(x-2nu) U_mu(x-2nu) U_nu(x+mu-2nu) U_nu(x+mu-nu)
      shifts.push_back(shift(0, -1));
      shifts.push_back(shift(0, -2));
      shifts.push_back(shift(0, -2));
      shifts.push_back(shift(1, -2));
      shifts.push_back(shift(1, -1));
    }
    return GeneralLocalStencil(ggrid, shifts);
  }

  static void RectStaplePaddedRs(GaugeLinkField &Cup, const GaugeField &gU,
				 const GeneralLocalStencil &gStencil, int mu, RealD rho) {
    GRID_TRACE("RectStaplePaddedRs");
    GridBase *ggrid = gU.Grid();
    conformable(ggrid, Cup.Grid());

    autoView( Cup_v , Cup, AcceleratorWrite);
    autoView( gU_v , gU, AcceleratorRead);
    autoView( gStencil_v, gStencil, AcceleratorRead);
    accelerator_for(ss, ggrid->oSites(), ggrid->Nsimd(), {
	typedef decltype(coalescedRead(Cup_v[0])) LinkMat;

	LinkMat tmp = Zero();
	for (int nu = 0; nu < Nd; nu++) {
	  int inc = 10*(nu - (mu<=nu));
	  if (nu != mu) {
	    GeneralStencilEntry const* e;
	    e = gStencil_v.GetEntry(0+inc,ss);
	    auto U_nu_x          =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (nu)();
	    e = gStencil_v.GetEntry(1+inc,ss);
	    auto U_nu_xpnu       =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (nu)();
	    e = gStencil_v.GetEntry(2+inc,ss);
	    auto U_mu_xp2nu      =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (mu)();
	    e = gStencil_v.GetEntry(3+inc,ss);
	    auto Udag_nu_xpmupnu = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
	    e = gStencil_v.GetEntry(4+inc,ss);
	    auto Udag_nu_xpmu    = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();

	    tmp()() = tmp()() + U_nu_x * U_nu_xpnu * U_mu_xp2nu * Udag_nu_xpmupnu * Udag_nu_xpmu;

	    e = gStencil_v.GetEntry(5+inc,ss);
	    auto Udag_nu_xmnu    = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
	    e = gStencil_v.GetEntry(6+inc,ss);
	    auto Udag_nu_xm2nu   = adj(coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd))(nu)();
	    e = gStencil_v.GetEntry(7+inc,ss);
	    auto U_mu_xm2nu      =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (mu)();
	    e = gStencil_v.GetEntry(8+inc,ss);
	    auto U_nu_xpmum2nu   =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (nu)();
	    e = gStencil_v.GetEntry(9+inc,ss);
	    auto U_nu_xpmumnu    =     coalescedReadGeneralPermute(gU_v[e->_offset], e->_permute, Nd) (nu)();

	    tmp()() = tmp()() + Udag_nu_xmnu * Udag_nu_xm2nu * U_mu_xm2nu * U_nu_xpmum2nu * U_nu_xpmumnu;
	  }
	}
	coalescedWrite(Cup_v[ss], rho*tmp);
      });
  }


protected:

  // Assume: SmearRhoRs is set
  void rectStapleRs(GaugeField &C, const std::vector<GaugeLinkField> &U, const std::vector<GaugeLinkField> &U2) const{
    
    GaugeLinkField Stap(C.Grid()), tmp(C.Grid());

    for (int mu = 0; mu < Nd; mu++) {
      Stap = Zero();
      for (int nu = 0; nu < Nd; nu++) {
	if (nu != mu) {
	  //      -<--
	  //      |  |
	  //
	  //      |  |
	  //      x
	  
	  // U_nu(x+mu)... => Umu(x)*Stap is a rectangle staple
	  tmp = Gimpl::CshiftLink(adj(U2[nu]), nu, -2);
	  tmp = Gimpl::CovShiftBackward(U[mu], mu, tmp);
	  tmp = U2[nu] * Gimpl::CshiftLink(tmp, nu, 2);
	  Stap += SmearRhoRs[mu + Nd * nu]*Gimpl::CshiftLink(tmp, mu, 1);
	  
	  //      |  |
	  //
	  //      |  |
	  //      -<--
	  
	  tmp = Gimpl::CovShiftBackward(U[mu], mu, U2[nu]);
	  tmp = adj(U2[nu]) * tmp;
	  tmp = Gimpl::CshiftLink(tmp, nu, -2);
	  Stap += SmearRhoRs[mu + Nd * nu]*Gimpl::CshiftLink(tmp, mu, 1);

	}
      }
      pokeLorentz(C, adj(Stap), mu); // C[mu] = Cup^dag   see conventions for Staple
    }
  }

  void rectStapleRl(GaugeField &C, const std::vector<GaugeLinkField> U, const std::vector<GaugeLinkField> &U2) const{
    
    GaugeLinkField Stap(C.Grid()), Staple2x1(C.Grid()), tmp(C.Grid());

    for (int mu = 0; mu < Nd; mu++) {
      Stap = Zero();
      for (int nu = 0; nu < Nd; nu++) {
        if (nu != mu) {
	  // Up staple    ___ ___
	  //             |       |
	  tmp = Gimpl::CshiftLink(adj(U[nu]), nu, -1);
	  tmp = adj(U2[mu]) * tmp;
	  tmp = Gimpl::CshiftLink(tmp, mu, -2);
	  
	  Staple2x1 = Gimpl::CovShiftForward(U[nu], nu, tmp);
	  
	  // Down staple
	  //             |___ ___|
	  //
	  tmp = adj(U2[mu]) * U[nu];
	  Staple2x1 += Gimpl::CovShiftBackward(U[nu], nu, Gimpl::CshiftLink(tmp, mu, -2));
	  
	  //              ___ ___
	  //             |    ___|
	  //             |___ ___|
	  //
	  
	  Stap += SmearRhoRl[mu + Nd * nu]*Gimpl::CshiftLink(Gimpl::CovShiftForward(U[mu], mu, Staple2x1), mu, 1);
	  
	  //              ___ ___
	  //             |___    |
	  //             |___ ___|
	  //
	  
	  Stap += SmearRhoRl[mu + Nd * nu]*Gimpl::CshiftLink(Staple2x1, mu, 1) * Gimpl::CshiftLink(U[mu], mu, -1);
	}
      }
      pokeLorentz(C, adj(Stap), mu); // C[mu] = Cup^dag   see conventions for Staple  
    }
  }
  
public:

  /*! Stout smearing with base explicitly specified */
  /* disable it; cannot access to rho
  Rect_Stout(Smear<Gimpl>* base) : SmearBase{base} {
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
  }
  */
  
  /*! Construct stout smearing object from explicitly specified rho matrix; Assume: SmearRhoP = SmearRhoRs = SmearRhoRl */
  /* disable it; this object should perform only one type of smearing.  Otherwise, deriv becomes complicated
  Rect_Stout(const std::vector<double>& rho_)
  : Smear_Stout<Gimpl>(rho_), SmearRhoRs(rho_),SmearRhoRl(rho_) {
  }
  */
  /*! Default constructor: rho is constant in all directions, optionally except for orthogonal dimension */
  Rect_Stout(double rho = 0.0, double rho_s = 1.0, double rho_l = 0.0, int orthogdim = -1, int orthogdim_s = -1, int orthogdim_l = -1)
    : Smear_Stout<Gimpl>(rho, orthogdim),
      OrthogDimRs{orthogdim_s}, SmearRhoRs{ this->rho3D(rho_s,orthogdim_s) },
      OrthogDimRl{orthogdim_l}, SmearRhoRl{ this->rho3D(rho_l,orthogdim_l) }{
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
    assert( (rho+rho_s == 0.0 || rho+rho_l == 0.0 || rho_s+rho_l == 0.0) && "Only one of the input rho values should be non-zero");
  }

  ~Rect_Stout() {}  // delete SmearBase...

  // Return: stout link = e^(iQ)U
  void smear(GaugeField& u_smr, const GaugeField& U) const {
    GaugeField C(U.Grid());
    GaugeLinkField tmp(U.Grid()), iq_mu(U.Grid()), Umu(U.Grid());

    std::cout << GridLogMessage << "Rect Stout smearing started\n";

    // C contains the staples multiplied by some rho
    C = Zero();
    if(this->SmearRho[1] > 0) { // If we specialize this class to only rectangler smearing, this case can be deleted
      this->SmearBase->smear(C, U); // Assume: SmearBase = Smear_APE 
      std::cout << GridLogMessage << "BaseSmearREct: " <<norm2(C)<<" "<<this->SmearRho[1]<<std::endl; // Assume: SmearBase = Smear_APE
    }
    else {
      WilsonLoops<Gimpl> WL;
      std::vector<GaugeLinkField> Us(Nd, U.Grid()),U2s(Nd, U.Grid());

      for (int mu = 0; mu < Nd; mu++) {
	Us[mu] = PeekIndex<LorentzIndex>(U, mu);
	WL.RectStapleDouble(U2s[mu], Us[mu], mu);
      }
      
      if(SmearRhoRs[1]>0) {
	rectStapleRs(C, Us, U2s);
	//std::cout << GridLogMessage << "BaseSmearREc t rs: " <<norm2(C)<<SmearRhoRs[1]<<std::endl;
      }
      else if(SmearRhoRl[1]>0) {
	rectStapleRl(C, Us, U2s);
	//std::cout << GridLogMessage << "BaseSmearREct rl: " <<norm2(C)<<" "<<SmearRhoRl[1]<<std::endl;
      }
    }
    
    u_smr = U; // set the smeared field to the current gauge field
    for (int mu = 0; mu < Nd; mu++) {
      if( mu == this->OrthogDim || mu == OrthogDimRs || mu == OrthogDimRl) continue ;
      // u_smr = exp(iQ_mu)*U_mu apart from Orthogdim
      Umu = peekLorentz(U, mu);
      tmp = peekLorentz(C, mu);
      iq_mu = Ta( tmp * adj(Umu));  
      this->exponentiate_iQ(tmp, iq_mu);
      pokeLorentz(u_smr, tmp * Umu, mu);
    }

  };

  void derivative(GaugeField& SigmaTerm, const GaugeField& iLambda,
                  const GaugeField& U) const {
    GridBase *grid = U.Grid();

    WilsonLoops<Gimpl> WL;
    GaugeLinkField staple(grid), u_tmp(grid);
    GaugeLinkField iLambda_mu(grid), iLambda_nu(grid);
    GaugeLinkField U_mu(grid), U_nu(grid);
    GaugeLinkField sh_field(grid), temp_Sigma(grid);
    Real rho_munu, rho_numu;

    if(this->SmearRho[1] > 0)
      this->SmearBase->derivative(SigmaTerm, iLambda, U);
    else if(SmearRhoRs[1]>0) {
      // Force from C = Rs in stout smearing
      for(int mu = 0; mu < Nd; ++mu){
	U_mu       = peekLorentz(      U, mu);
	iLambda_mu = peekLorentz(iLambda, mu);
      
	for(int nu = 0; nu < Nd; ++nu){
	  if(nu==mu) continue;
	  
	  U_nu       = peekLorentz(U, nu);
	  iLambda_nu = peekLorentz(iLambda, nu);
	  
	  rho_munu = SmearRhoRs[mu + Nd * nu];
	  rho_numu = SmearRhoRs[nu + Nd * mu];
	  
	  RectStapleUnoptimisedRsUpper(staple, U, mu, nu);
	  
	  // 1st
	  temp_Sigma = -rho_numu*staple*iLambda_nu;
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	  
	  // 2nd
	  sh_field = adj(U_mu)*Cshift(temp_Sigma*U_mu, mu, -1);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);

	  // 3rd
	  sh_field = Cshift(iLambda_nu, mu, 1);
	  temp_Sigma = rho_numu*sh_field*adj(U_mu)*Cshift(staple*U_mu,mu,-1);
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	  
	  // 4th
	  sh_field = Cshift(U_mu*temp_Sigma, mu,1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);
	  
	  // 5th
	  temp_Sigma = rho_numu*adj(U_mu)*Cshift(adj(staple*U_nu)*adj(U_mu)*iLambda_nu*U_nu,nu,-1);
	  Gimpl::AddLink(SigmaTerm, temp_Sigma, mu);
	  
	  // 6th
	  sh_field = adj(U_mu)*Cshift(temp_Sigma*U_mu,mu,-1);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);
	
	  // 7th
	  temp_Sigma = -rho_numu*Cshift(Cshift(adj(U_nu)*iLambda_nu*adj(adj(U_nu)*Cshift(adj(U_mu)*Cshift(staple*U_mu,mu,-1)*U_mu,mu,-1)),mu,1),nu,-1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm,temp_Sigma, mu);
	  
	  // 8th
	  sh_field = Cshift(U_mu*temp_Sigma,mu,1)*adj(U_mu);
	  Gimpl::AddLink(SigmaTerm, sh_field, mu);
	  
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
    } else if (SmearRhoRl[1]>0) {
      // TODO: implement for the case C = Rl
      assert(1 != 1 && "Longer rectangle flow is not implemeted as yet");
    }
  
  };

  // Retrun: fat link with staples
  void BaseSmear(GaugeField& C, const GaugeField& U) const {
    GaugeLinkField tmp(U.Grid());
    std::vector<GaugeLinkField> Us(Nd, U.Grid()),U2s(Nd, U.Grid());
    WilsonLoops<Gimpl> WL;

    //std::cout << GridLogMessage << "Rect Stout smearing base \n";
    for (int mu = 0; mu < Nd; mu++) {
      Us[mu] = PeekIndex<LorentzIndex>(U, mu);
      WL.RectStapleDouble(U2s[mu], Us[mu], mu);
    }

    if(this->SmearRho[1] > 0) this->SmearBase->smear(C, U);
    else if (SmearRhoRs[1]>0) rectStapleRs(C, Us, U2s);
    else if (SmearRhoRl[1]>0) rectStapleRl(C, Us, U2s);
  };
  

};

NAMESPACE_END(Grid);
