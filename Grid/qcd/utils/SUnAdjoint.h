#ifndef QCD_UTIL_SUNADJOINT_H
#define QCD_UTIL_SUNADJOINT_H

////////////////////////////////////////////////////////////////////////
//
// * Adjoint representation generators
//
// * Normalisation for the fundamental generators: 
//   trace ta tb = 1/2 delta_ab = T_F delta_ab
//   T_F = 1/2  for SU(N) groups
//    - s.y.) diff from the convention in https://en.wikipedia.org/wiki/Gell-Mann_matrices
//    - s.y.) tb = i T^b of Luscher (\because tb = \lambda/2 from SUn.impl.h)
//
//   base for NxN hermitian traceless matrices
//   normalized to 1:
//
//   (e_Adj)^a = t^a / sqrt(T_F)
//
//   then the real, antisymmetric generators for the adjoint representations
//   are computed ( shortcut: e^a == (e_Adj)^a )
//     - s.y.) If basis matrices of adj. rep. are real & anti-symmetric, it cannot be Hermitian
//     - s.y.) If T_adj is purely imaginary & antisymmetric, then  Hermitian
//     - s.y.) Here, it computes: i*(T^c)_adj^{Grid} = - (T^c)_adj^{Lushcer}   
//
//   (iT_adj)^d_ba = i tr[e^a t^d e^b - t^d e^a e^b]
//
////////////////////////////////////////////////////////////////////////

NAMESPACE_BEGIN(Grid);

template <int ncolour>
class SU_Adjoint : public SU<ncolour> {
public:
  static const int Dimension = ncolour * ncolour - 1;

  template <typename vtype>
  using iSUnAdjointMatrix =
    iScalar<iScalar<iMatrix<vtype, Dimension > > >;

  // Actually the adjoint matrices are real...
  // Consider this overhead... FIXME
  typedef iSUnAdjointMatrix<Complex> AMatrix;
  typedef iSUnAdjointMatrix<ComplexF> AMatrixF;
  typedef iSUnAdjointMatrix<ComplexD> AMatrixD;

  typedef iSUnAdjointMatrix<vComplex>  vAMatrix;
  typedef iSUnAdjointMatrix<vComplexF> vAMatrixF;
  typedef iSUnAdjointMatrix<vComplexD> vAMatrixD;

  typedef Lattice<vAMatrix>  LatticeAdjMatrix;
  typedef Lattice<vAMatrixF> LatticeAdjMatrixF;
  typedef Lattice<vAMatrixD> LatticeAdjMatrixD;

  typedef Lattice<iVector<iScalar<iMatrix<vComplex, Dimension> >, Nd> >  LatticeAdjField;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexF, Dimension> >, Nd> > LatticeAdjFieldF;
  typedef Lattice<iVector<iScalar<iMatrix<vComplexD, Dimension> >, Nd> > LatticeAdjFieldD;


  template <typename vtype>
  using iSUnMatrix = iScalar<iScalar<iMatrix<vtype, ncolour> > >; // SU_N = [1][1][ncolour x ncolour][len(vtype)]

  typedef Lattice<iScalar<iScalar<iVector<vComplex, Dimension> > > >  LatticeAdjVector;

  template <class cplx>
  static void generator(int Index, iSUnAdjointMatrix<cplx> &iAdjTa) {
    // returns i(T_Adj)^index necessary for the projectors
    // s.y.:
    //  - T_Adj^a is the adj rep of t^a in the basis t^a's
    //   => i(T_Adj) is -(adj rep) in Luscher's basis (T^a = -it^a) => i(T_Adj) is adj rep in the basis -T^a
    // see definitions above
    iAdjTa = Zero();
    Vector<iSUnMatrix<cplx> > ta(ncolour * ncolour - 1);
    iSUnMatrix<cplx> tmp;

    // FIXME not very efficient to get all the generators everytime
    for (int a = 0; a < Dimension; a++) SU<ncolour>::generator(a, ta[a]);

    for (int a = 0; a < Dimension; a++) {
      tmp = ta[a] * ta[Index] - ta[Index] * ta[a]; // = -[t^c, t^a] = -(T^c)_adj[T^a] <- [(T^c)_adj]^{ab} = - [(T^c)_adj]^{ba} = - <(T^c)_adj T^a, T^b> 
      for (int b = 0; b < (ncolour * ncolour - 1); b++) { // ncolour * ncolour - 1 = #basis in su(N)
        iSUnMatrix<cplx> tmp1 = 2.0 * tmp * ta[b];  // 2.0 from the normalization; <A, B> = -2 Tr[AB]_Luscher = 2 Tr[AB]_here <- tb = i T^b of Luscher
	// s.y.) Up to this point, basis of su(3) follow Grid's convention
        Complex iTr = TensorRemove(timesI(trace(tmp1)));
        //iAdjTa()()(b, a) = iTr;<-?
        iAdjTa()()(a, b) = iTr; // = i*(T^c)_adj^{Grid} = - (T^c)_adj^{Lushcer}
      }
    }
  }

  static void printGenerators(void) {
    for (int gen = 0; gen < Dimension; gen++) {
      AMatrix ta;
      generator(gen, ta);
      std::cout << GridLogMessage << "Nc = " << ncolour << " t_" << gen
                << std::endl;
      std::cout << GridLogMessage << ta << std::endl;
    }
  }

  static void testGenerators(void) {
    AMatrix adjTa;
    std::cout << GridLogMessage << "Adjoint - Checking if real" << std::endl;
    for (int a = 0; a < Dimension; a++) {
      generator(a, adjTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adjTa - conjugate(adjTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;

    std::cout << GridLogMessage << "Adjoint - Checking if antisymmetric"
              << std::endl;
    for (int a = 0; a < Dimension; a++) {
      generator(a, adjTa);
      std::cout << GridLogMessage << a << std::endl;
      assert(norm2(adjTa + transpose(adjTa)) < 1.0e-6);
    }
    std::cout << GridLogMessage << std::endl;
  }

  static void AdjointLieAlgebraMatrix(
				      const typename SU<ncolour>::LatticeAlgebraVector &h,
				      LatticeAdjMatrix &out, Real scale = 1.0) {
    conformable(h, out);
    GridBase *grid = out.Grid();
    LatticeAdjMatrix la(grid);
    AMatrix iTa;

    out = Zero();
    for (int a = 0; a < Dimension; a++) {
      generator(a, iTa);
      la = peekColour(h, a) * iTa;
      out += la;
    }
    out *= scale;
  }

  // Projects the algebra components a lattice matrix (of dimension ncol*ncol -1 )
  static void projectOnAlgebra(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeAdjMatrix &in, Real scale = 1.0) 
  {
    
    conformable(h_out, in);
    h_out = Zero();
    AMatrix iTa;
    Real coefficient = - 1.0/(ncolour) * scale;// 1/Nc for the normalization of the trace in the adj rep

    for (int a = 0; a < Dimension; a++) {
      generator(a, iTa);
      pokeColour(h_out, real(trace(iTa * in)) * coefficient, a);
    }
  }

  // a projector that keeps the generators stored to avoid the overhead of recomputing them 
  static void projector(typename SU<ncolour>::LatticeAlgebraVector &h_out, const LatticeAdjMatrix &in, Real scale = 1.0) {
    conformable(h_out, in);
    static std::vector<AMatrix> iTa(Dimension);  // to store the generators
    h_out = Zero();
    static bool precalculated = false; 
    if (!precalculated){
      precalculated = true;
      for (int a = 0; a < Dimension; a++) generator(a, iTa[a]);
    }

    Real coefficient = -1.0 / (ncolour) * scale;  // 1/Nc for the normalization of
    // the trace in the adj rep

    for (int a = 0; a < Dimension; a++) {
      auto tmp = real(trace(iTa[a] * in)) * coefficient; 
      pokeColour(h_out, tmp, a);
    }
  }


};




// Some useful type names

typedef SU_Adjoint<2> SU2Adjoint;
typedef SU_Adjoint<3> SU3Adjoint;
typedef SU_Adjoint<4> SU4Adjoint;
typedef SU_Adjoint<5> SU5Adjoint;

typedef SU_Adjoint<Nc> AdjointMatrices;

NAMESPACE_END(Grid);

#endif
