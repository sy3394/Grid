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
  
*/
template <class Gimpl>
class Rect_Stout : public Smear_Stout<Gimpl> {
protected:
  int OrthogDimRs = -1;
  int OrthogDimRl = -1; // TODO: implement this vertical rectangle smearing
  
public:
  INHERIT_GIMPL_TYPES(Gimpl)
    
  const std::vector<double> SmearRhoRs, SmearRhoRl;
  
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
	  
	  tmp = Gimpl::CshiftLink(adj(U2[nu]), nu, -2);
	  tmp = Gimpl::CovShiftBackward(U[mu], mu, tmp);
	  tmp = U2[nu] * Gimpl::CshiftLink(tmp, nu, 2);
	  Stap += Gimpl::CshiftLink(tmp, mu, 1);
	  
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
      pokeLorentz(C, Stap, mu);
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
	  
	  Stap += Gimpl::CshiftLink(Gimpl::CovShiftForward(U[mu], mu, Staple2x1), mu, 1);
	  
	  //              ___ ___
	  //             |___    |
	  //             |___ ___|
	  //
	  
	  Stap += SmearRhoRl[mu + Nd * nu]*Gimpl::CshiftLink(Staple2x1, mu, 1) * Gimpl::CshiftLink(U[mu], mu, -1);
	}
      }
      pokeLorentz(C, Stap, mu);
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
  Rect_Stout(const std::vector<double>& rho_)
  : Smear_Stout<Gimpl>(rho_), SmearRhoRs(rho_),SmearRhoRl(rho_) {
  }

  /*! Default constructor: rho is constant in all directions, optionally except for orthogonal dimension */
  Rect_Stout(double rho = 1.0, double rho_s = 1.0, double rho_l = 1.0, int orthogdim = -1, int orthogdim_s = -1, int orthogdim_l = -1)
    : Smear_Stout<Gimpl>(rho, orthogdim),
      OrthogDimRs{orthogdim_s}, SmearRhoRs{ this->rho3D(rho_s,orthogdim_s) },
      OrthogDimRl{orthogdim_l}, SmearRhoRl{ this->rho3D(rho_l,orthogdim_l) }{
    assert(Nc == 3 && "Stout smearing currently implemented only for Nc==3");
  }

  ~Rect_Stout() {}  // delete SmearBase...

  // Return: stout link = e^(iQ)U
  void smear(GaugeField& u_smr, const GaugeField& U) const {
    GaugeField C(U.Grid()), C_tmp(U.Grid());
    GaugeLinkField tmp(U.Grid()), iq_mu(U.Grid()), Umu(U.Grid());
    std::vector<GaugeLinkField> Us(Nd, U.Grid()),U2s(Nd, U.Grid());
    WilsonLoops<Gimpl> WL;

    for (int mu = 0; mu < Nd; mu++) {
      Us[mu] = PeekIndex<LorentzIndex>(U, mu);
      WL.RectStapleDouble(U2s[mu], Us[mu], mu);
    }
    
    std::cout << GridLogDebug << "Rect Stout smearing with Plq + Rect (Rs + Rl) started\n";

    // C contains the staples multiplied by some rho
    u_smr = U; // set the smeared field to the current gauge field
    this->SmearBase->smear(C, U); // Assume: SmearBase = Smear_APE
    rectStapleRs(C_tmp, Us, U2s);
    C = C + C_tmp;
    rectStapleRl(C_tmp, Us, U2s);
    C =	C + C_tmp;
    for (int mu = 0; mu < Nd; mu++) {
      if( mu == this->OrthogDim || mu == OrthogDimRs || mu == OrthogDimRl) continue ;
      // u_smr = exp(iQ_mu)*U_mu apart from Orthogdim
      Umu = peekLorentz(U, mu);
      tmp = peekLorentz(C, mu);
      iq_mu = Ta( tmp * adj(Umu));  
      this->exponentiate_iQ(tmp, iq_mu);
      pokeLorentz(u_smr, tmp * Umu, mu);
    }

    std::cout << GridLogDebug << "Rect Stout smearing with Plq + Rect (Rs + Rl) completed\n";
  };

  // TODO: fix the below to take into account Rs, Rl
  void derivative(GaugeField& SigmaTerm, const GaugeField& iLambda,
                  const GaugeField& Gauge) const {
    this->SmearBase->derivative(SigmaTerm, iLambda, Gauge);
  };

  // Retrun: fat link with staples
  void BaseSmear(GaugeField& C, const GaugeField& U) const {
    GaugeLinkField tmp(U.Grid());
    std::vector<GaugeLinkField> Us(Nd, U.Grid()),U2s(Nd, U.Grid());
    WilsonLoops<Gimpl> WL;

    for (int mu = 0; mu < Nd; mu++) {
      Us[mu] = PeekIndex<LorentzIndex>(U, mu);
      WL.RectStapleDouble(U2s[mu], Us[mu], mu);
    }

    this->SmearBase->smear(C, U);
    rectStapleRs(tmp, Us, U2s);
    C = C + tmp;
    rectStapleRl(tmp, Us, U2s);
    C = C + tmp;    
  };
  
  // TODO: implement derivative, BaseSmear for Rs, Rl for completeness <- not yet done as not used for FTHMC

};

NAMESPACE_END(Grid);
