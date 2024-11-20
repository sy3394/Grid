/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_trace.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_LATTICE_TRACE_H
#define GRID_LATTICE_TRACE_H

///////////////////////////////////////////////
// Tracing, transposing, peeking, poking
///////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace
////////////////////////////////////////////////////////////////////////////////////////////////////
/*
template<class vobj>
inline auto trace(const Lattice<vobj> &lhs)  -> Lattice<decltype(trace(vobj()))>
{
  Lattice<decltype(trace(vobj()))> ret(lhs.Grid());
  autoView(ret_v , ret, AcceleratorWrite);
  autoView(lhs_v , lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), vobj::Nsimd(), {
    coalescedWrite(ret_v[ss], trace(lhs_v(ss)));
  });
  return ret;
};
*/
    
////////////////////////////////////////////////////////////////////////////////////////////////////
// Trace Index level dependent operation
////////////////////////////////////////////////////////////////////////////////////////////////////
template<int Index,class vobj>
inline auto TraceIndex(const Lattice<vobj> &lhs) -> Lattice<decltype(traceIndex<Index>(vobj()))>
{
  Lattice<decltype(traceIndex<Index>(vobj()))> ret(lhs.Grid());
  autoView( ret_v , ret, AcceleratorWrite);
  autoView( lhs_v , lhs, AcceleratorRead);
  accelerator_for( ss, lhs_v.size(), vobj::Nsimd(), {
    coalescedWrite(ret_v[ss], traceIndex<Index>(lhs_v(ss)));
  });
  return ret;
};

template<int N, class Vec>
Lattice<iScalar<iScalar<iScalar<Vec> > > > Determinant(const Lattice<iScalar<iScalar<iMatrix<Vec, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto lvol = grid->lSites();
  Lattice<iScalar<iScalar<iScalar<Vec> > > > ret(grid);
  typedef typename Vec::scalar_type scalar;
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,lvol,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);
    Coordinate lcoor;
    grid->LocalIndexToLocalCoor(site, lcoor);
    iScalar<iScalar<iMatrix<scalar, N> > > Us;
    peekLocalSite(Us, Umu_v, lcoor);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	scalar tmp= Us()()(i,j);
	ComplexD ztmp(real(tmp),imag(tmp));
	EigenU(i,j)=ztmp;
      }}
    ComplexD detD  = EigenU.determinant();
    typename Vec::scalar_type det(detD.real(),detD.imag());
    pokeLocalSite(det,ret_v,lcoor);
  });
  return ret;
}

template<int N, class Vec>
Lattice<iScalar<iScalar<iScalar<Vec> > > > Determinant_real(const Lattice<iScalar<iScalar<iMatrix<Vec, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iScalar<Vec> > > > ret(grid);
  ret.Checkerboard() = Umu.Checkerboard();
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
      Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);
      
      iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
      iScalar<iScalar<iScalar<ComplexD> > > Ud;
      for(int lane=0;lane<Nsimd;lane++){
	Us = extractLane(lane,Umu_v[site]);
	for(int i=0;i<N;i++){
	  for(int j=0;j<N;j++){
	    EigenU(i,j)=real(Us()()(i,j));
	  }}
	Ud()()() = EigenU.determinant();
	insertLane(lane,ret_v[site],Ud);
      }
    });
  return ret;
}

template<int N>
Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > Inverse(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
{
#if 0
  GridBase *grid=Umu.Grid();
  auto lvol = grid->lSites();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
  
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,lvol,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);
    Coordinate lcoor;
    grid->LocalIndexToLocalCoor(site, lcoor);
    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;
    peekLocalSite(Us, Umu_v, lcoor);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	EigenU(i,j) = Us()()(i,j);
      }}
    Eigen::MatrixXcd EigenUinv = EigenU.inverse();
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	Ui()()(i,j) = EigenUinv(i,j);
      }}
    pokeLocalSite(Ui,ret_v,lcoor);
  });
#else
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
    Eigen::MatrixXcd EigenU = Eigen::MatrixXcd::Zero(N,N);

    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;

    for(int lane=0;lane<Nsimd;lane++){
      Us = extractLane(lane,Umu_v[site]);
      for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
	  EigenU(i,j) = Us()()(i,j);
	}}
      Eigen::MatrixXcd EigenUinv = EigenU.inverse();
      for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
	  Ui()()(i,j) = EigenUinv(i,j);
	}}
      insertLane(lane,ret_v[site],Ui);
    }
  });
#endif  
  return ret;
}

template<int N>
Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > Inverse_RealPart(const Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > &Umu)
{
  GridBase *grid=Umu.Grid();
  auto osites = grid->oSites();
  const int Nsimd=grid->Nsimd();
  Lattice<iScalar<iScalar<iMatrix<vComplexD, N> > > > ret(grid);
#if 1
  autoView(Umu_v,Umu,CpuRead);
  autoView(ret_v,ret,CpuWrite);
  thread_for(site,osites,{
    Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);

    iScalar<iScalar<iMatrix<ComplexD, N> > > Us;
    iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;

    for(int lane=0;lane<Nsimd;lane++){
      Us = extractLane(lane,Umu_v[site]);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          EigenU(i,j) = real(Us()()(i,j));
        }}
      Eigen::MatrixXd EigenUinv = EigenU.inverse();
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          Ui()()(i,j) = EigenUinv(i,j);
        }}
      insertLane(lane,ret_v[site],Ui);
    }
  });
#else // Eigen supports inversion on GPU's only for matrices of size < 5
  iScalar<iScalar<iMatrix<ComplexD, N> > > Ui;
  autoView(Umu_v,Umu,AcceleratorRead);
  autoView(ret_v,ret,AcceleratorWrite);
  accelerator_for(ss,grid->oSites(),Nsimd,{
      Eigen::MatrixXd EigenU = Eigen::MatrixXd::Zero(N,N);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          EigenU(i,j) = real(Umu_v(ss)()()(i,j));
        }}
      //Linker error occurs w/r/t Eigen when combining the below two lines into a one liner
      Eigen::MatrixXd EigenUinv; 
      EigenUinv= EigenU.inverse();
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          Ui()()(i,j) = EigenUinv(i,j);
        }}
      coalescedWrite(ret_v[ss],Ui);
    });
#endif
 return ret;
}

NAMESPACE_END(Grid);
#endif

