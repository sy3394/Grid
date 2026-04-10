// Derived from VTK/Examples/Cxx/Medical2.cxx
// The example reads a volume dataset, extracts two isosurfaces that
// represent the skin and bone, and then displays them.
//
// Modified heavily by Peter Boyle to display lattice field theory data as movies and compare multiple files

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkMetaImageReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkStripper.h>
#include <vtkImageData.h>
#include <vtkVersion.h>
#include <vtkCallbackCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#define MPEG
#ifdef MPEG
#include <vtkFFMPEGWriter.h>
#endif

#include <vtkProperty2D.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkWindowToImageFilter.h>

#include <array>
#include <string>
#include <fstream>
#include <sstream>

#include <Grid/Grid.h>
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

int mpeg = 0 ;
int Ls = -1;
int xlate = 0 ;
int take_diff = 0;
int omit_dir = 4;
int dynm_dir = 3;
std::vector<std::string> dynm_labels = {"X", "Y", "Z", "T", "tau"};

// HDF5 reader: dataset "field" shape (Nsites,2) float64, Grid lex order (x fastest)
// Uses the HDF5 C API (hdf5.h) — works even when C++ bindings (H5Cpp.h) are absent.
template <class T> void readFileHDF5(T& out, std::string const fname){
#if defined(HAVE_HDF5)
  typedef typename T::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  GridBase *grid   = out.Grid();
  int64_t   Nsites = grid->_gsites;

  hid_t fid = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); assert(fid >= 0);
  hid_t did = H5Dopen2(fid, "field", H5P_DEFAULT);                  assert(did >= 0);

  std::vector<double> buf(2 * Nsites);
  herr_t err = H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  assert(err >= 0);
  H5Dclose(did); H5Fclose(fid);

  std::vector<sobj> lexbuf(Nsites);
  for(int64_t i = 0; i < Nsites; i++)
    lexbuf[i]()()() = Grid::ComplexD(buf[2*i], buf[2*i+1]);

  Grid::vectorizeFromLexOrdArray(lexbuf, out);
  std::cout << Grid::GridLogMessage << "readFileHDF5: loaded " << fname << std::endl;
#endif
}

static bool isHDF5file(std::string const &fname){
  if(fname.size() > 3 && fname.substr(fname.size()-3) == ".h5")   return true;
  if(fname.size() > 5 && fname.substr(fname.size()-5) == ".hdf5") return true;
  return false;
}

template <class T> void readFile(T& out, std::string const fname){
  if(isHDF5file(fname)){
    readFileHDF5(out, fname);
    return;
  }
#ifdef HAVE_LIME
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
#endif
}
template <class T> void writeFile(T& in, std::string const fname){
#ifdef HAVE_LIME
  // Ref: https://github.com/paboyle/Grid/blob/feature/scidac-wp1/tests/debug/Test_general_coarse_hdcg_phys48.cc#L111
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in,record,0); // Lexico
  WR.close();
#endif
}

using namespace Grid;

int n_dims = Nd; //can remove this if define it in class FrameUpdater based on ext_latt_size  


int main(int argc, char* argv[])
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian    UGrid(latt_size, simd_layout, mpi_layout);
  GridCartesian * grid  = &UGrid, *gridF;
 
  std::cout << argc << " command Line arguments "<<std::endl;
  for(int c=0;c<argc;c++) {
    std::cout << " - "<<argv[c]<<std::endl;
  }

  /*
  For now,
    - xlate: the first non-omitted coor
    - the last coordinate is dynamic, i.e., to be upated
    - default: omit_dir = 4 && dynm_dir = 3
         diff. simulation time corresp. to diff. frame & time is updated
    - if omit_dir < 4;
        - index on input files is treated as a dimension of the frame
	- if dynm_dir = 4; 3 out of 4 latt dim is displayed
	- if dynm_dir < 4; only 2 out of 4 latt dim is displayed + sim. time
 */
  int save_file = 0, PCF = 0, CTCDs = 0, shifting=0, sum_all=0;
  double cut=0;
  std::ofstream data_f;
  std::vector<std::string> file_list1, file_list2;
  double default_contour = 1.0;
  std::vector<int> corner(4,0);
  std::vector<int> b_size(4,0);
  std::vector<int> shift(4,0);
  std::vector<int> n_0modes;
  std::string arg;
  std::vector<std::string> save_fnames;
#ifdef MPEG
  std::string mpeg_fname = "movie.avi";
  if( GridCmdOptionExists(argv,argv+argc,"--mpeg") ){
    mpeg = 1;
    std::string fname = GridCmdOptionPayload(argv,argv+argc,"--mpeg");
    if(!fname.empty()) mpeg_fname = fname;
  }
#endif

  if( GridCmdOptionExists(argv,argv+argc,"--Ls") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--Ls");
    GridCmdOptionInt(arg,Ls);
    assert(Ls !=0 && "Ls needs to be greater than 0 when specified");
    gridF = SpaceTimeGrid::makeFiveDimGrid(Ls, &UGrid);
    latt_size = gridF->GlobalDimensions();
    dynm_labels.insert(dynm_labels.begin(),"Ls");
    //omit_dirs[0] = 5;
    n_dims++;
    std::cout<<gridF->GlobalDimensions()<<" "<<latt_size<<std::endl;
  }
    
  if( GridCmdOptionExists(argv,argv+argc,"--xlate") ){
    xlate = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--omit_dir") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--omit_dir");
    GridCmdOptionInt(arg,omit_dir);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--dynm_dir") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--dynm_dir");
    GridCmdOptionInt(arg,dynm_dir);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--compute_PCF") ){
    PCF = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--compareTCD_defs") ){
    CTCDs = 1;
  }
    
  if( GridCmdOptionExists(argv,argv+argc,"--take_adj_diff") ){
    take_diff = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--sum_all_files") ){
    sum_all = 1;
    arg = GridCmdOptionPayload(argv,argv+argc,"--sum_all_files");
    GridCmdOptionCSL(arg, save_fnames);  // comma-separated: fname1[,fname2]
  }

  if( GridCmdOptionExists(argv,argv+argc,"--shift") ){
    shifting=1;
    arg=GridCmdOptionPayload(argv,argv+argc,"--shift");
    GridCmdOptionIntVector(arg,shift);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--box_corner") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--box_corner");
    GridCmdOptionIntVector(arg,corner);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--box_size") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--box_size");
    GridCmdOptionIntVector(arg,b_size);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--cut") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--cut");
    GridCmdOptionFloat(arg,cut);
  }
    
  if( GridCmdOptionExists(argv,argv+argc,"--isosurface") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--isosurface");
    GridCmdOptionFloat(arg,default_contour);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--files1") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--files1");
    std::cout<<"files1 "<<arg<<std::endl;
    GridCmdOptionCSL(arg,file_list1);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--files2") ){ // For eigen
    arg = GridCmdOptionPayload(argv,argv+argc,"--files2");
    std::cout<<"files2 "<<arg<<std::endl;
    GridCmdOptionCSL(arg,file_list2);
    //assert(file_list1.size() == file_list2.size() && "files1 and files2 must have the same number of files");
  }
  if( GridCmdOptionExists(argv,argv+argc,"--chiral_modes") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--chiral_modes");
    n_0modes.resize(file_list2.size());
    GridCmdOptionIntVector(arg,n_0modes);
  }

  if( GridCmdOptionExists(argv,argv+argc,"--save_data_to") ){
    save_file = 1;
    arg = GridCmdOptionPayload(argv,argv+argc,"--save_data_to");
    if(!arg.empty()) data_f.open(arg,std::ios::trunc);
  }

  // ---------- Topo charge density from eigenvectors ----------
  // evals: eigenvalues of H_DWF = Gamma5R5 * D_DWF(mass).
  //   These come from Compute_DWF_G5R5.cc (eMe[i]) and INCLUDE m_f.
  //   Near-zero modes have |mu_n| ~ m_f; bulk modes have |mu_n| ~ O(1).
  //   --evals accepts either:
  //     (a) a comma-separated list:  --evals 0.012,0.015,-0.011
  //     (b) a path to a text file:   --evals /path/to/evals.txt
  //         File format: one eigenvalue per line (blank lines / '#' comments ignored).
  //   The number of eigenvalues should match the number of --f2 density files;
  //   if fewer are provided the remaining modes are treated as mu_n=0 (sign=+1).
  // mass_f: input quark mass m_f used when generating the eigenvectors.
  // m_gap:  spectral gap m_gap = m_f + m_res (physical quark mass including residual
  //         chiral symmetry breaking).  Used as the weight in Formula A'.
  //         If not supplied via --m_gap it is estimated as min_n|mu_n| from --evals.
  // topo_out: output file prefix for the four q_top fields.
  double mass_f = 0.0;
  double m_gap  = -1.0;   // sentinel: negative means "not yet set"
  std::vector<double> evals;
  std::string topo_out = "topo_evec";

  if( GridCmdOptionExists(argv,argv+argc,"--mass") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mass");
    GridCmdOptionFloat(arg, mass_f);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--evals") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--evals");
    // Detect file vs inline list: try parsing the first token as a double.
    // If that fails (or arg contains a '/' or '.txt'), treat as a filename.
    bool is_file = false;
    {
      // Heuristic: if arg contains a path separator or cannot be parsed as a
      // number at all, assume it is a file path.
      std::istringstream probe(arg);
      double v; probe >> v;
      is_file = probe.fail() || (arg.find('/') != std::string::npos)
                              || (arg.find('\\') != std::string::npos);
    }
    if(is_file){
      // Read one eigenvalue per line; skip blank lines and '#' comments.
      std::ifstream fin(arg);
      if(!fin) { std::cerr << "ERROR: cannot open evals file: " << arg << std::endl; exit(1); }
      std::string line;
      while(std::getline(fin, line)){
        // Strip inline comments
        auto pos = line.find('#');
        if(pos != std::string::npos) line = line.substr(0, pos);
        std::istringstream iss(line);
        double v;
        while(iss >> v) evals.push_back(v);
      }
      std::cout << "Loaded " << evals.size() << " eigenvalues from file: " << arg << std::endl;
    } else {
      // Comma-separated inline list, e.g. --evals 0.012,0.015,-0.011
      std::vector<std::string> eval_strs;
      GridCmdOptionCSL(arg, eval_strs);
      for(auto& s: eval_strs) evals.push_back(std::stod(s));
      std::cout << "Loaded " << evals.size() << " eigenvalues (inline)" << std::endl;
    }
  }
  if( GridCmdOptionExists(argv,argv+argc,"--topo_out") ){
    topo_out = GridCmdOptionPayload(argv,argv+argc,"--topo_out");
  }
  if( GridCmdOptionExists(argv,argv+argc,"--m_gap") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--m_gap");
    GridCmdOptionFloat(arg, m_gap);
    std::cout << "m_gap set explicitly: " << m_gap << std::endl;
  }
  // If m_gap was not supplied, estimate it as min_n |mu_n|.
  // This approximates m_f + m_res = m_phys (the spectral gap of H_DWF).
  if(m_gap < 0.0 && !evals.empty()){
    m_gap = std::fabs(evals[0]);
    for(auto& v : evals) if(std::fabs(v) < m_gap) m_gap = std::fabs(v);
    std::cout << "m_gap auto-estimated as min|mu_n| = " << m_gap << std::endl;
  } else if(m_gap < 0.0){
    m_gap = mass_f;  // last-resort fallback
    std::cout << "m_gap fallback to mass_f = " << m_gap << std::endl;
  }
  // ------------------------------------------------------------

  assert( omit_dir != dynm_dir && "The omitted dir cannot be the same as updated dimension" );

  /**************   Read in Files   ************************************/
  /*
    For now,
       data1: 4D fields
       data2: initially 5D fields -> summed over Ls -> ends up as 4D fields
   */
  // 4D fields: data1
  FieldMetaData header;
  std::vector<LatticeComplexD> data1(file_list1.size()-take_diff,grid);
  for(int c=0;c<data1.size();c++) {
    std::cout << "Reading file1: "<<file_list1[c]<<std::endl;
    readFile(data1[c],file_list1[c]);
    std::cout<<"Sum "<<c<<" "<<real(TensorRemove(sum(data1[c])))<<std::endl;
  }
  if(take_diff){
    for(int c=0;c<data1.size()-1;c++)
      data1[c] = data1[c+1] - data1[c];
    LatticeComplexD tmp(data1[0].Grid());
    std::cout << "Reading file1: "<<file_list1.back()<<std::endl;
    readFile(tmp,file_list1.back());
    std::cout<<"Sum last "<<real(TensorRemove(sum(tmp)))<<std::endl;
    data1.back() = tmp - data1.back();
  }

  // 5D fields: data2
  //   data2[c]: plain sum over s (existing behaviour, used for visualisation)
  //   q_eps, q_bdy, q_mid: three q_top formulas accumulated over eigenvectors
  bool compute_topo = (Ls > 0) && !file_list2.empty();
  LatticeComplexD q_eps(grid), q_bdy(grid), q_mid(grid), q_prime(grid);
  if(compute_topo){ q_eps = Zero(); q_bdy = Zero(); q_mid = Zero(); q_prime = Zero(); }

  std::vector<LatticeComplexD> data2(file_list2.size()-take_diff,grid);
  for(int c=0;c<data2.size();c++) {
    std::cout << "Reading file2: "<<file_list2[c]<<std::endl;
    LatticeComplexD tmp(gridF);
    readFile(tmp,file_list2[c]);
    if(Ls > 0){
      // 5D input: sum over the Ls dimension to produce a 4D density
      LatticeComplexD tmp4D(grid); tmp4D = Zero(); data2[c] = Zero();
      for(int i=0; i<Ls;i++){
        ExtractSlice(tmp4D,tmp,i,0);
        data2[c] = data2[c] + tmp4D;
      }
    } else {
      // 4D input: use directly (e.g. pre-computed fermion TCD definition files)
      data2[c] = tmp;
    }
    std::cout<<"Sum "<<c<<" "<<real(TensorRemove(sum(data2[c])))<<std::endl;

    // ==================== Topo charge density (3 formulas) ====================
    //
    // Background (Blum et al. 2004, arXiv:hep-lat/0105006):
    //
    //   D_H = gamma5 * R5 * D_DWF    (Eq. 12, Hermitian DWF operator)
    //   D_H psi_n = mu_n psi_n        (real eigenvalues mu_n)
    //   H_DWF^2 = D_DWF† D_DWF       => D†D evecs phi_k are related to psi_n
    //                                    (Compute_DWF_G5R5.cc rotates phi->psi)
    //
    // Input files (file_list2) contain the 5D scalar density per mode:
    //   rho_n(x,s) = |psi_n(x,s)|^2  (sum over spinor & color indices)
    // as written by the "evec_density" output of Compute_DWF_G5R5.cc.
    //
    // Scalar sign function eps_code(s)  [= -Gamma5_Blum, Blum Eq. 17]:
    //   eps_code(s) = -1  for s <  Ls/2   (left-wall side)
    //   eps_code(s) = +1  for s >= Ls/2   (right-wall side)
    // This sign convention matches G5evec construction in Compute_DWF_G5R5.cc
    // (axpby_ssp applies factor -1 for s < Ls/2, +1 for s >= Ls/2).
    //
    // evals[c] = mu_n = eigenvalue of H_DWF (includes m_f).
    //   Near-zero (topological) modes: |mu_n| ~ m_f  -> sign(mu_n) captures chirality.
    //   Bulk (+/-mu) pairs:            contributions cancel in the sum.
    //
    // WHY sign(mu_n) instead of m_f/mu_n:
    //   The exact formula contains the factor m_f/mu_n, but for near-zero modes
    //   |mu_n| ~ m_f (up to m_res), so m_f/mu_n ~ sign(mu_n).  Using sign(mu_n)
    //   directly avoids explicit m_f dependence and works correctly even at m_f=0.
    //
    // -----------------------------------------------------------------------
    // FORMULA B:  eps_code(s)-chirality density   (Bulk form, from Blum Eq. 17 / tr[Gamma5])
    // -----------------------------------------------------------------------
    //   q_B(x) = sum_n sign(mu_n) * sum_s eps_code(s) * rho_n(x,s)
    //
    //   Derivation:
    //     Q = -(1/2) sum_n sign(mu_n) <psi_n|Gamma5_Blum|psi_n>
    //       =  (1/2) sum_n sign(mu_n) sum_s eps_code(s) |psi_n(x,s)|^2   (summed over x)
    //   Here Gamma5_Blum = delta_{ss'} sgn((Ls-1)/2 - s)  (Blum Eq. 17, SCALAR sign,
    //   not gamma5*R5).  The factor 1/2 is absorbed into the normalization convention;
    //   we accumulate without it to match the other two formulas.
    //   Near-zero modes localized at one wall contribute +/-1; bulk symmetric modes
    //   cancel because eps_code sums to zero over a uniformly distributed mode.
    //
    // -----------------------------------------------------------------------
    // FORMULA C:  boundary projection density      (Boundary form, related to Blum Eq. 8)
    // -----------------------------------------------------------------------
    //   q_C(x) = -sum_n sign(mu_n) * [rho_n(x,Ls-1) - rho_n(x,0)]
    //
    //   Derivation:
    //     q_top(x) = -m_f * tr[gamma5 * S^{4D}(x,x)]
    //              = -m_f * sum_n (1/mu_n) * [psi_n†(x,Ls-1) P_R psi_n(x,Ls-1)
    //                                        + psi_n†(x,0)   P_L psi_n(x,0)]
    //   where P_R = (1+gamma5)/2, P_L = (1-gamma5)/2  with the STANDARD 4D gamma5,
    //   and the PLUS sign between the two wall terms is exact (it arises from the
    //   gamma5 factor inside tr[gamma5 D^{-1}] when D^{-1} is expressed via D_H^{-1}).
    //   Replacing m_f/mu_n -> sign(mu_n) and the spinor-projected densities by the
    //   scalar proxy rho(x,s) = |psi(x,s)|^2:
    //     psi†(x,Ls-1) P_R psi(x,Ls-1) ~ rho(x,Ls-1)  [right-wall mode: P_R ~ 1]
    //     psi†(x,0)    P_L psi(x,0)    ~ rho(x,0)      [left-wall mode:  P_L ~ 1]
    //   The relative minus sign in (rho(Ls-1) - rho(0)) encodes the eps_code
    //   convention (eps(Ls-1)=+1, eps(0)=-1).
    //
    // -----------------------------------------------------------------------
    // FORMULA A:  midpoint density                 (Midpoint form, analog of Blum Eq. 9)
    // -----------------------------------------------------------------------
    //   q_A(x) = -sum_n sign(mu_n) * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
    //
    //   Derivation:
    //     Analogous to Formula C but using the midpoint slices s=Ls/2 and s=Ls/2-1
    //     instead of the physical walls s=Ls-1 and s=0.  Corresponds to the
    //     "midpoint axial current" J^a_{5q,mid}(x) defined in Blum Eq. (9).
    //     The factor 0.5 normalises relative to the wall formula.
    // ===========================================================================
    if(compute_topo) {
      // Eigenvalue of H_DWF = gamma5*R5*D_DWF (includes m_f, passed via --evals).
      double mu_n    = (c < (int)evals.size()) ? evals[c] : 0.0;
      // sign(mu_n): +1 or -1.  For near-zero modes, m_f/mu_n ~ sign(mu_n)*1.
      // Bulk (+/-mu) pairs cancel automatically in the sum over modes.
      double sign_mu = (mu_n >= 0.0) ? 1.0 : -1.0;
      // m_gap/mu_n weight for Formula A' (q_prime).
      // Uses the spectral gap m_gap = m_f + m_res instead of sign(mu_n):
      //   - topological modes: |mu_n| ~ m_gap  =>  w_Ap ~ sign(mu_n)  (same as q_A)
      //   - bulk modes:        |mu_n| ~ Lambda  =>  w_Ap ~ m_gap/Lambda << 1  (suppressed)
      // Guard against division by zero (should never occur with physical evals).
      double w_Ap = (mu_n != 0.0) ? (m_gap / mu_n) : 0.0;

      // --- Formula B:  sign(mu_n)-weighted eps_code chirality sum ---
      // --- Formula B': m_gap/mu_n-weighted eps_code chirality sum  ---
      // Both accumulate sum_s eps_code(s) * rho_n(x,s); they differ only in weight.
      // eps_code(s) = -1 for s < Ls/2  (left),  +1 for s >= Ls/2  (right).
      {
        LatticeComplexD eps_slice(grid);
        for(int s = 0; s < Ls; s++){
          ExtractSlice(eps_slice, tmp, s, 0);
          double eps_s = (s >= Ls/2) ? 1.0 : -1.0;   // eps_code(s)
          q_eps   = q_eps   + (sign_mu * eps_s) * eps_slice;  // Formula B
          q_prime = q_prime + (w_Ap    * eps_s) * eps_slice;  // Formula B'
        }
      }

      // --- Formula C: boundary projection ---
      // q_C(x) += -sign(mu_n) * [rho_n(x,Ls-1) - rho_n(x,0)]
      // s=Ls-1: right wall (eps_code=+1),  s=0: left wall (eps_code=-1).
      // Scalar proxy for the exact spinor boundary formula (see comment above).
      {
        LatticeComplexD bdy_s0(grid), bdy_sLs(grid);
        ExtractSlice(bdy_s0,  tmp, 0,    0);   // left wall,  s=0
        ExtractSlice(bdy_sLs, tmp, Ls-1, 0);   // right wall, s=Ls-1
        q_bdy = q_bdy - sign_mu * (bdy_sLs - bdy_s0);
      }

      // --- Formula A: midpoint density ---
      // q_A(x) += -sign(mu_n) * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
      // Midpoint slices straddle the 5D bulk; analogous to Formula C at mid-plane.
      // Corresponds to Blum Eq. (9) J^a_{5q,mid}(x).
      if(Ls >= 2){
        LatticeComplexD mid_lo(grid), mid_hi(grid);
        ExtractSlice(mid_lo, tmp, Ls/2-1, 0);   // below midpoint (eps_code=-1)
        ExtractSlice(mid_hi, tmp, Ls/2,   0);   // above midpoint (eps_code=+1)
        q_mid = q_mid - sign_mu * 0.5 * (mid_hi - mid_lo);
      }

      std::cout << "TopoContrib evec=" << c
                << " mu_n=" << mu_n << " sign(mu_n)=" << sign_mu
                << " w_Bp=" << w_Ap
                << " Q_B="  << real(TensorRemove(sum(q_eps)))
                << " Q_B'=" << real(TensorRemove(sum(q_prime)))
                << " Q_C="  << real(TensorRemove(sum(q_bdy)))
                << " Q_A="  << real(TensorRemove(sum(q_mid))) << std::endl;
    }
    // ==========================================================================
  }
  if(take_diff){
    for(int c=0;c<data2.size()-1;c++)
      data2[c] = data2[c+1] - data2[c];
    LatticeComplexD tmp(data2[0].Grid());
    std::cout << "Reading file2: "<<file_list2.back()<<std::endl;
    readFile(tmp,file_list2.back());
    std::cout<<"Sum last "<<real(TensorRemove(sum(tmp)))<<std::endl;
    data2.back() = tmp - data2.back();
  }
  if(sum_all && !save_fnames.empty()){
    if(!data1.empty()){
      LatticeComplexD tmp(data1[0].Grid()); tmp = Zero();
      for(int c=0;c<data1.size();c++) tmp = tmp + data1[c];
      writeFile(tmp, save_fnames[0]);
    }
    if(!data2.empty()){
      // Use save_fnames[1] when provided (both datasets); else save_fnames[0] (single-dataset usage)
      std::string fname2 = (save_fnames.size() >= 2) ? save_fnames[1] : save_fnames[0];
      LatticeComplexD tmp(data2[0].Grid()); tmp = Zero();
      for(int c=0;c<data2.size();c++) tmp = tmp + data2[c];
      writeFile(tmp, fname2);
    }
  }
  
  /****** Write topo charge density reconstructed from eigenvectors (4 formulas) *****/
  // Output file naming: pass --topo_out with a {def} placeholder, e.g.
  //   --topo_out /path/Top_dnsty_q_{def}_0_smr.702
  // C++ replaces {def} with A, B, Bp, C to produce the four output files.
  // q_A, q_B, q_C use sign(mu_n); q_Bp uses m_gap/mu_n (requires --evals).
  if(compute_topo && !evals.empty()){
    // Lambda: substitute {def} placeholder in topo_out template string.
    auto fill_def = [](std::string tmpl, const std::string& d) -> std::string {
      auto pos = tmpl.find("{def}");
      if(pos != std::string::npos) tmpl.replace(pos, 5, d);
      return tmpl;
    };

    // Formula B: eps_code(s)-chirality density  [sign(mu_n) weight, bulk form]
    //   q_B(x) = sum_n sign(mu_n) * sum_s eps_code(s) * rho_n(x,s)
    writeFile(q_eps, fill_def(topo_out,"B"));
    std::cout << "Wrote q_B   -> " << fill_def(topo_out,"B")
              << "  Q_B=" << real(TensorRemove(sum(q_eps))) << std::endl;

    // Formula B': m_gap-weighted chirality density  [m_gap/mu_n weight, bulk improved]
    //   q_B'(x) = sum_n (m_gap/mu_n) * sum_s eps_code(s) * rho_n(x,s)
    //   Bulk modes suppressed by m_gap/Lambda_bulk << 1.
    //   Recommended for pointwise comparison with gradient-flowed q^gf(x).
    writeFile(q_prime, fill_def(topo_out,"Bp"));
    std::cout << "Wrote q_B'  -> " << fill_def(topo_out,"Bp")
              << "  Q_B'=" << real(TensorRemove(sum(q_prime)))
              << "  (m_gap=" << m_gap << ")" << std::endl;

    // Formula C: boundary projection density  [Blum Eq. 8 scalar proxy]
    //   q_C(x) = -sum_n sign(mu_n) * [rho_n(x,Ls-1) - rho_n(x,0)]
    writeFile(q_bdy, fill_def(topo_out,"C"));
    std::cout << "Wrote q_C   -> " << fill_def(topo_out,"C")
              << "  Q_C=" << real(TensorRemove(sum(q_bdy))) << std::endl;

    // Formula A: midpoint density  [Blum Eq. 9 analog]
    //   q_A(x) = -sum_n sign(mu_n) * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
    writeFile(q_mid, fill_def(topo_out,"A"));
    std::cout << "Wrote q_A   -> " << fill_def(topo_out,"A")
              << "  Q_A=" << real(TensorRemove(sum(q_mid))) << std::endl;
  }
  /****** IP & Corr of each fermion TCD definition vs gluonic TCD (--topo_compare) *****/
  // Requires: --topo_out (so q_eps/q_prime/q_bdy/q_mid are available),
  //           --files1   (one file per gluonic flow time, e.g. TD_tau=0,4,16),
  //           --conf_id  (integer config number, written to output for notebook parsing).
  // Output per line (pure numeric): TD_tau  tau  conf  Corr  IP
  //   where TD_tau is the index i of files1[i], tau comes from --conf_id context,
  //   and the label "Topo_PCF: q_X" is printed to stderr for diagnostics only.
  // q_A, q_B, q_C use sign(mu_n) (Approximation 1).
  // The 4 output blocks (one per q_def) are labelled by a comment line "# q_X"
  // so the shell can split them into 4 per-def files via grep.
  int conf_id = -1;
  if( GridCmdOptionExists(argv,argv+argc,"--conf_id") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--conf_id");
    GridCmdOptionInt(arg, conf_id);
  }
  int topo_compare = GridCmdOptionExists(argv,argv+argc,"--topo_compare");
  if(compute_topo && !evals.empty() && topo_compare && !data1.empty()){
    typedef typename PeriodicGimplR::ComplexField ComplexField;
    // data1[i] = gluonic TCD at flow time index i  (one file per TD_tau)
    struct QDef { std::string name; LatticeComplexD* field; };
    std::vector<QDef> qdefs = {{"q_A",&q_mid},{"q_B",&q_eps},{"q_Bp",&q_prime},{"q_C",&q_bdy}};
    for(auto& qd : qdefs){
      std::cout << "# " << qd.name << std::endl;  // shell uses this to split output
      for(int i=0; i<(int)data1.size(); i++){
        LatticeComplex X(grid), Y(grid), one(grid); one = ComplexField::scalar_type(1.0,0.0);
        ComplexD avg1 = TensorRemove(sum(data1[i]))/RealD(grid->gSites());
        ComplexD avg2 = TensorRemove(sum(*qd.field))/RealD(grid->gSites());
        X = data1[i] - avg1*one;
        Y = *qd.field  - avg2*one;
        double corr = TensorRemove(sum(X*Y)).real()/sqrt(norm2(X)*norm2(Y));
        X = data1[i]; Y = *qd.field;
        double ip   = TensorRemove(innerProduct(X,Y)).real()/sqrt(norm2(X))/sqrt(norm2(Y));
        // format: TD_tau_index  conf  corr   (one Corr line + one IP line, matching old format)
        std::cout << "Topo PCF Corr: " << i << " " << conf_id << " " << corr << std::endl;
        std::cout << "Topo PCF IP:   " << i << " " << conf_id << " " << ip   << std::endl;
      }
    }
  }
  /****** Compare fermion TCD definitions vs qlat reference fields (--comp_file) ***/
  // --comp_file <file>[,<file2>,...]
  //   SCIDAC files from qlat: each is a separate stochastic DWF TCD estimate
  //   (e.g. topo_field_0.scidac, topo_field_1.scidac from pickle_to_scidac.ipynb).
  //   Each file is compared individually against ALL 4 fermion TCD definitions.
  //   After comparison, every comp field is appended to data1 so it appears as
  //   an additional frame in the visualisation pipeline (FieldDensityAnimateMultiFiles).
  //
  //   Output per line: "CompRef: def comp_idx conf Q_evec Q_ref Corr IP rms_diff"
  std::vector<LatticeComplexD> comp_fields;   // populated below; appended to data1 at end
  if(compute_topo && !evals.empty() && GridCmdOptionExists(argv,argv+argc,"--comp_file")){
    arg = GridCmdOptionPayload(argv,argv+argc,"--comp_file");
    std::vector<std::string> comp_fnames;
    GridCmdOptionCSL(arg, comp_fnames);

    typedef typename PeriodicGimplR::ComplexField ComplexField;
    LatticeComplexD one(grid); one = ComplexField::scalar_type(1.0, 0.0);

    struct QDef { std::string name; LatticeComplexD* field; };
    std::vector<QDef> qdefs = {{"q_A",&q_mid},{"q_B",&q_eps},{"q_Bp",&q_prime},{"q_C",&q_bdy}};

    // Header for human-readable table
    std::cout << GridLogMessage
              << std::left  << std::setw(5)  << "ci"
              << std::setw(6)  << "def"
              << std::right << std::setw(12) << "Q_evec"
                            << std::setw(12) << "Q_ref"
                            << std::setw(12) << "Corr"
                            << std::setw(12) << "IP"
                            << std::setw(12) << "rms_diff" << std::endl;

    for(int ci = 0; ci < (int)comp_fnames.size(); ci++){
      LatticeComplexD ref(grid);
      readFile(ref, comp_fnames[ci]);
      double   Q_ref   = real(TensorRemove(sum(ref)));
      ComplexD avg_ref = TensorRemove(sum(ref)) / RealD(grid->gSites());
      std::cout << GridLogMessage << "comp_file[" << ci << "] " << comp_fnames[ci]
                << "  Q=" << Q_ref << std::endl;

      for(auto& qd : qdefs){
        double   Q_evec   = real(TensorRemove(sum(*qd.field)));
        ComplexD avg_evec = TensorRemove(sum(*qd.field)) / RealD(grid->gSites());

        // Pearson correlation (mean-subtracted)
        LatticeComplexD X(grid), Y(grid);
        X = *qd.field - avg_evec * one;
        Y = ref       - avg_ref  * one;
        double corr = real(TensorRemove(sum(X*Y))) / std::sqrt(norm2(X) * norm2(Y));

        // Normalised inner product
        X = *qd.field;  Y = ref;
        double ip = real(TensorRemove(innerProduct(X,Y))) / std::sqrt(norm2(X)) / std::sqrt(norm2(Y));

        // RMS pointwise difference
        LatticeComplexD diff(grid); diff = *qd.field - ref;
        double rms_diff = std::sqrt(norm2(diff) / RealD(grid->gSites()));

        std::cout << GridLogMessage
                  << std::left  << std::setw(5)  << ci
                  << std::setw(6)  << qd.name
                  << std::right << std::setw(12) << Q_evec
                                << std::setw(12) << Q_ref
                                << std::setw(12) << corr
                                << std::setw(12) << ip
                                << std::setw(12) << rms_diff << std::endl;

        // Machine-readable for shell/notebook: "CompRef: def comp_idx conf Q_evec Q_ref Corr IP rms_diff"
        std::cout << "CompRef: " << qd.name << " " << ci << " " << conf_id << " "
                  << Q_evec << " " << Q_ref << " " << corr << " " << ip << " " << rms_diff << std::endl;
      }

      // Stash for appending to data1 (visualisation frames) below
      comp_fields.push_back(ref);
    }
  }
  /******************************************************************************/

  /****************   Filter Gluonic TCD by Fermionic Chiral Eigen Modes  ******/
  if( CTCDs ){
    std::vector<int> shift_c(4,0);
    for(int x=0;x<shift[0];x++)
      for(int y=0;y<shift[1];y++)
	for(int z=0;z<shift[2];z++)
	  for(int t=0;t<shift[3];t++){
	    shift_c[0]=x;shift_c[1]=y;shift_c[2]=z;shift_c[3]=t;
	    
	    for(int i=0; i<data1.size(); i++) {
	      typedef typename PeriodicGimplR::ComplexField ComplexField;
	      ComplexField filter(grid), ones(grid), zeros(grid); ones = ComplexField::scalar_type(1.0,0.0); zeros = Zero();
	      LatticeReal tmp(grid);
	      tmp = toReal(localNorm2(data2[i]));
	      RealD cut2= cut*maxLocalNorm2(data2[i]);
	      filter = where( tmp > cut2, ones, zeros);
	      for(int mu=0;mu<Nd;mu++) 
		filter = Cshift(filter,mu,shift_c[mu]);
	      //LatticeComplexD tmp2(grid); tmp2 = toComplex(filter);
	      //std::cout<<"Filtered Sum: " << i <<" "<< real(TensorRemove(sum(data1[i]*filter))) << std::endl;
	      std::cout<<"Filtered Sum: " << shift_c <<" "<< i <<" "<<real(innerProduct(data1[i],filter))
		       << " max norm " << std::sqrt(maxLocalNorm2(data2[i]))
		       <<" avg filter "<<real(sum(filter))/RealD(grid->gSites())
		       <<"filter sum  "<<TensorRemove(sum(filter))
		       <<" volume "<<RealD(grid->gSites())<<std::endl;//" "<<maxLocalNorm2(tmp2)<<std::endl;
	    }
    }
  }
  
  /****************   Compute Corr Coeff & IP  ********************************/
  if(PCF){
    assert(data1.size()==data2.size());
    typedef typename PeriodicGimplR::ComplexField ComplexField;
    
    std::vector<ComplexD> avgs1;
    for(int i=0; i<data1.size(); i++) avgs1.push_back(TensorRemove(sum(data1[i]))/RealD(grid->gSites()));
    std::vector<ComplexD> avgs2;
    for(int i=0; i<data2.size(); i++) avgs2.push_back(TensorRemove(sum(data2[i]))/RealD(grid->gSites()));

    for(int i=0; i<data1.size(); i++){
      LatticeComplex X(grid), Y(grid), one(grid); one = ComplexField::scalar_type(1.0,0.0);
      X = data1[i] - avgs1[i]*one;
      Y = data2[i] - avgs2[i]*one;
      
      std::cout<< "Corr Coeff: " << i << " " << TensorRemove(sum(X*Y)).real()/sqrt(norm2(X)*norm2(Y)) << std::endl;

      X = data1[i];
      Y = data2[i];

      std::cout<< "Inner Product: " << i<<" "<<TensorRemove(innerProduct(X,Y)).real()/sqrt(norm2(X))/sqrt(norm2(Y)) << std::endl;
    }
  }

  data_f.close();

  /****** Append comp fields to data1 for visualisation ************************/
  // comp fields are appended after all fixed-size analyses (PCF, CTCDs) so their
  // presence does not disturb assertions on data1.size().
  // In FieldDensityAnimateMultiFiles they appear as additional frames alongside
  // the gluonic TCD fields; labels: "comp_0", "comp_1", ...
  if(!comp_fields.empty()){
    std::cout << GridLogMessage << "Appending " << comp_fields.size()
              << " comp field(s) to data1 for visualisation." << std::endl;
    for(int ci = 0; ci < (int)comp_fields.size(); ci++){
      data1.push_back(comp_fields[ci]);
      std::cout << GridLogMessage << "  data1[" << data1.size()-1 << "] = comp_" << ci
                << "  Q=" << real(TensorRemove(sum(data1.back()))) << std::endl;
    }
  }
  /******************************************************************************/

  Grid_finalize();

  return EXIT_SUCCESS;
}
