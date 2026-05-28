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
  Grid::GridBase *grid   = out.Grid();
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

// Raw binary reader: big-endian complex128, Grid lex order (x fastest), no header.
// Python: arr.T.astype(np.dtype('>c16')).tofile("field.bin")
template <class T> void readFileBinary(T& out, std::string const fname){
  typedef typename T::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  uint32_t nersc_csum=0, scidac_csuma=0, scidac_csumb=0;
  Grid::BinarySimpleMunger<sobj,sobj> munge;
  Grid::BinaryIO::readLatticeObject<vobj,sobj>(out, fname, munge, 0, "IEEE64BIG",
                                               nersc_csum, scidac_csuma, scidac_csumb);
  std::cout << Grid::GridLogMessage << "readFileBinary: loaded " << fname << std::endl;
}

static bool isRawBinaryFile(std::string const &fname){
  if(fname.size() > 4 && fname.substr(fname.size()-4) == ".bin") return true;
  if(fname.size() > 4 && fname.substr(fname.size()-4) == ".raw") return true;
  return false;
}

template <class T> void readFile(T& out, std::string const fname){
  if(isHDF5file(fname)){
    readFileHDF5(out, fname);
    return;
  }
  if(isRawBinaryFile(fname)){
    readFileBinary(out, fname);
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

// H_DWF_EvalRecord: embedded in each evec_density SCIDAC file by Compute_DWF_G5R5.cc
namespace Grid {
  struct H_DWF_EvalRecord : Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(H_DWF_EvalRecord,
      double, eval,   // eigenvalue of H_DWF = gamma5*R5*D_DWF(mass)
      int,    n       // mode index: pairs (+mu_0,-mu_0,+mu_1,-mu_1,...) ordered by |mu_n|
    );
  };
}

// Read a SCIDAC field and return the embedded H_DWF_EvalRecord.
// Sets eval_out=0, n_out=-1 if the file format does not carry a record.
template <class T>
void readFileRecord(T& out, std::string const fname, double& eval_out, int& n_out){
  eval_out = 0.0; n_out = -1;
  if(isHDF5file(fname) || isRawBinaryFile(fname)){
    readFile(out, fname);   // no user-record in these formats
    return;
  }
#ifdef HAVE_LIME
  Grid::H_DWF_EvalRecord rec;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out, rec);
  RD.close();
  eval_out = rec.eval;
  n_out    = rec.n;
  std::cout << Grid::GridLogMessage
            << "readFileRecord: " << fname
            << "  eval=" << eval_out << "  n=" << n_out << std::endl;
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

// ---------------------------------------------------------------------------
// gaussian_smear:  apply a 4D Gaussian low-pass filter to a complex scalar
// lattice field via FFT.   K(p) = exp(-sigma^2 |p|^2 / 2)
// with continuum lattice momentum p_mu = 2*pi*k_mu/L_mu, k_mu folded into
// [-L/2, L/2).
//
// Implements the Luchang/collaborator suggestion (sec:smear_sweep): convolving
// q_B^mgap and q_gluon with the same sigma decouples the gluonic Wilson-flow's
// non-linear gauge diffusion (which preserves instanton saddle points) from
// pure linear signal-processing bandwidth matching.   Symmetric application
// to both fields is required to test the "matched-bandwidth peak" hypothesis.
//
// Cost: ~1 forward + 1 backward 4D FFT per call + a per-site kernel build.
// The kernel is recomputed inside the function for cleanliness; if used in a
// hot loop callers should cache the kernel themselves.
// ---------------------------------------------------------------------------
LatticeComplexD gaussian_smear(const LatticeComplexD& field, double sigma){
  GridBase* gridB = field.Grid();
  GridCartesian* gridC = dynamic_cast<GridCartesian*>(gridB);
  assert(gridC && "gaussian_smear: requires GridCartesian");
  Coordinate latt = gridC->FullDimensions();
  int Nd_local = latt.size();

  // Build kernel = exp(-sigma^2/2 * sum_mu (2*pi*k_mu/L_mu)^2)
  // Iterate per local site with pokeLocalSite — runs once per sigma and
  // avoids per-lane SIMD acrobatics; cost is O(V) doubles per call.
  LatticeComplexD kernel(gridC);
  kernel = ComplexD(0.0, 0.0);
  {
    Coordinate local_dims = gridC->LocalDimensions();
    Coordinate proc_coor  = gridC->ThisProcessorCoor();
    double half_sig2 = 0.5 * sigma * sigma;
    for(int ls = 0; ls < gridC->lSites(); ls++){
      Coordinate lcoor;
      gridC->LocalIndexToLocalCoor(ls, lcoor);
      double psq = 0.0;
      for(int mu = 0; mu < Nd_local; mu++){
        int km = proc_coor[mu] * local_dims[mu] + lcoor[mu];
        if(km > latt[mu]/2) km -= latt[mu];
        double pm = 2.0 * M_PI * km / RealD(latt[mu]);
        psq += pm * pm;
      }
      ComplexD val(std::exp(-half_sig2 * psq), 0.0);
      pokeLocalSite(val, kernel, lcoor);
    }
  }

  // Apply via FFT
  LatticeComplexD field_FT(gridC), out_field(gridC);
  FFT theFFT(gridC);
  theFFT.FFT_all_dim(field_FT, field, FFT::forward);
  field_FT = field_FT * kernel;
  theFFT.FFT_all_dim(out_field, field_FT, FFT::backward);
  return out_field;
}

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
  // WeightSpec: one weight track to accumulate simultaneously.
  //   label    : output file suffix and table key (e.g. "sign", "mgap", "mext")
  //   is_sign  : true  => w = sign(mu_n)  (exact ±1; label keyword "sign")
  //   mgap_val : used when !is_sign:
  //              = 0   => w = m_gap_auto/mu_n  (auto = min_n|mu_n| from --evals)
  //              ≠ 0   => w = mgap_val/mu_n    (user-supplied literal)
  // --weights <token>[,<token>,...] controls which tracks are active (default: sign,mgap).
  //   token = "sign"       => sign track
  //   token = "mgap"       => auto-mgap track
  //   token = "label=val"  => named literal track (e.g. mext=0.011)
  // --m_gap <val>          => legacy alias; adds/updates a track named "mext".
  //
  // evals: eigenvalues of H_DWF = Gamma5R5 * D_DWF(mass).
  //   These come from Compute_DWF_G5R5.cc (eMe[i]) and INCLUDE m_f.
  //   Near-zero modes have |mu_n| ~ m_f; bulk modes have |mu_n| ~ O(1).
  //   --evals accepts either:
  //     (a) a comma-separated list:  --evals 0.012,0.015,-0.011
  //     (b) a path to a text file:   --evals /path/to/evals.txt
  //         File format: one eigenvalue per line (blank lines / '#' comments ignored).
  //   The number of eigenvalues should match the number of --f2 density files;
  //   if fewer are provided the remaining modes are treated as mu_n=0 (sign=0).
  // mass_f: input quark mass m_f; fallback for m_gap_auto when --evals is absent.
  // topo_out: output file prefix for the q_top fields.
  struct WeightSpec {
    std::string label;
    bool        is_sign;   // true  => w = sign(mu_n)
    double      mgap_val;  // false => w = mgap_val/mu_n  (0 = auto: min|mu_n|, >0 = literal)
  };
  // Two distinct fermion masses are needed (see comment above on q_naive vs Sigma_low):
  //   mass_f  : MUST match the m_f used in Compute_DWF_G5R5 to generate the
  //             eigenvectors.  Used to recover the kinetic eigenvalue
  //             lambda^(0) = sqrt((lambda^H)^2 - m_f^2).  For zero-mass evecs
  //             (lambda^H already kinetic), set --mass 0.
  //   bc_mass : independent Banks–Casher Lorentzian regulator for Sigma_low(x).
  //             Defaults to mass_f if --bc_mass is absent (legacy behaviour).
  //             Must be non-zero to get a non-trivial Sigma_low.
  double mass_f  = 0.0;
  double bc_mass = -1.0;  // sentinel: -1 -> fall back to mass_f below
  // --n_topo <k>: number of near-zero (topological) eigenmodes to EXCLUDE from
  // the direct bulk remainder B_bulk(x) = sum_{n in bulk} chi_n^B(x).  The k
  // excluded modes are chosen by SMALLEST |mu_n| (computed below), NOT by index
  // position, so this is robust to an unsorted --evals file (e.g. conf 795
  // lists the zero mode at index 1, not 0).  Set k = |Q| (or 2|Q| if the
  // solver lists each zero mode as a degenerate +-pair); k=0 (default) makes
  // B_bulk the full unweighted chiral sum sum_n chi_n^B.
  // NOTE: this assumes evals[c] corresponds to eigenvector c.  If the --evals
  // file is ordered differently from the evec-density files, q_B^{mgap,sign}
  // and B_bulk are all mismatched -- re-sort the eval file to match the evecs
  // (or fix Compute_DWF_G5R5) before trusting the Gamma_5-weighted fields.
  int    n_topo  = 0;
  // mgap-track weight: DEFAULT uses the UNSIGNED eigenvalue, w = m_gap/|mu_n|,
  // so the zero-mode weight is +m_gap/|mu_0| and the chirality chi_0^B (which
  // carries the index, int chi_0 = Q) sets the sign of q_B -> int q_B = Q for
  // either sign of Q.  The legacy SIGNED weight w = m_gap/mu_n is sign-blind:
  // for the near-zero mode sgn(mu_0) = -Q (spectral flow), so m_gap/mu_0 = -Q
  // multiplies chi_0 = +Q to give -Q^2 = -1 regardless of Q (verified: conf 702
  // Q=-1 and conf 795 Q=+1 both give int q_B^mgap = -1).  --signed_mgap
  // restores the legacy signed weight for cross-checks only.
  bool   signed_mgap = false;
  std::vector<double> evals;
  std::string topo_out = "topo_evec";
  // Default: both sign and mgap tracks always active
  std::vector<WeightSpec> weight_specs = {{"sign",true,0.0},{"mgap",false,0.0}};

  if( GridCmdOptionExists(argv,argv+argc,"--mass") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mass");
    GridCmdOptionFloat(arg, mass_f);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--bc_mass") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--bc_mass");
    GridCmdOptionFloat(arg, bc_mass);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--signed_mgap") ){
    signed_mgap = true;
    std::cout << GridLogMessage << "--signed_mgap: using LEGACY signed weight "
              << "m_gap/mu_n (sign-blind to Q; cross-check only)" << std::endl;
  } else {
    std::cout << GridLogMessage << "mgap weight: m_gap/|mu_n| (default; "
              << "signed-index-faithful)" << std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--n_topo") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--n_topo");
    GridCmdOptionInt(arg, n_topo);
    std::cout << GridLogMessage << "--n_topo: excluding the " << n_topo
              << " leading near-zero mode(s) from B_bulk" << std::endl;
  }
  if(bc_mass < 0.0) bc_mass = mass_f;  // legacy fallback
  if(bc_mass == 0.0){
    std::cerr << "WARNING: bc_mass=0 (no --bc_mass and --mass=0).  Sigma_low(x) "
              << "will be identically zero and all PCFs involving it will be "
              << "NaN.  Pass --bc_mass <m_f> for a non-trivial Banks–Casher "
              << "Lorentzian regulator." << std::endl;
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
  bool do_topo_out = GridCmdOptionExists(argv,argv+argc,"--topo_out");
  if( do_topo_out ){
    topo_out = GridCmdOptionPayload(argv,argv+argc,"--topo_out");
  }
  // --weights: override the set of active weight tracks
  if( GridCmdOptionExists(argv,argv+argc,"--weights") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--weights");
    std::vector<std::string> tokens;
    GridCmdOptionCSL(arg, tokens);
    weight_specs.clear();
    for(auto& tok : tokens){
      auto eq = tok.find('=');
      if(eq != std::string::npos){
        std::string base = tok.substr(0, eq);
        double val = std::stod(tok.substr(eq+1));
        // Auto-encode value into label: mext=0.011 → label "mext0.011"
        // so the output filename is self-documenting (q_A_mext0.011.702).
        // Strip trailing zeros: format as shortest decimal that round-trips.
        std::ostringstream ss; ss << val;
        std::string lbl = base + ss.str();
        weight_specs.push_back({lbl, false, val});
        std::cout << GridLogMessage << "--weights: track '" << lbl
                  << "' mgap_val=" << val << std::endl;
      } else if(tok == "sign"){
        weight_specs.push_back({"sign", true, 0.0});
        std::cout << GridLogMessage << "--weights: track 'sign' (w=sign(mu_n))" << std::endl;
      } else if(tok == "mgap"){
        weight_specs.push_back({"mgap", false, 0.0});
        std::cout << GridLogMessage << "--weights: track 'mgap' (w=m_gap_auto/mu_n)" << std::endl;
      } else {
        std::cerr << "WARNING: unknown --weights token '" << tok << "' — ignored" << std::endl;
      }
    }
  }
  // --m_gap <val>: legacy alias; adds/updates track named "mext" with mgap_val=val
  if( GridCmdOptionExists(argv,argv+argc,"--m_gap") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--m_gap");
    double m_gap_ext; GridCmdOptionFloat(arg, m_gap_ext);
    bool found = false;
    for(auto& ws : weight_specs)
      if(ws.label == "mext"){ ws.is_sign = false; ws.mgap_val = m_gap_ext; found = true; break; }
    if(!found) weight_specs.push_back({"mext", false, m_gap_ext});
    std::cout << GridLogMessage << "--m_gap: added/updated 'mext' track with mgap_val="
              << m_gap_ext << std::endl;
  }
  // --band_pass <lc1,lc2,...>: add Gaussian even-window weight tracks
  //   w_{lc}(mu) = exp(-mu^2 / (2 lc^2))
  // an integer-exact member of the bulk-estimator family (topo_charge.tex
  // sec:qBfamily; sp_sum App. window_tests): the window is EVEN in mu and ->1 at
  // mu=0, so sum_x q_B^{w} = Q_ov + O(e^{-a Ls}) is preserved (the +-lambda^H
  // pair integrals cancel for any common weight), while the LOCAL content is
  // band-passed -- bulk-mode I-Ibar chirality kept, deep-UV jitter rolled off at
  // the mobility edge lc.  Encoded WITHOUT a WeightSpec struct change by a
  // NEGATIVE mgap_val sentinel: mgap_val = -lc flags a window track with
  // lambda_c = -mgap_val (handled in the track_w computation below).  These add
  // ordinary tracks, so they flow automatically into the q_{A,B,C}_<lbl> output,
  // the FermFerm/topo_compare qdefs, and the corr_ip/comp_ref files.  Default: none.
  if( GridCmdOptionExists(argv,argv+argc,"--band_pass") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--band_pass");
    std::vector<std::string> bptok;
    GridCmdOptionCSL(arg, bptok);
    for(auto& tk : bptok){
      double lc = std::stod(tk);
      if(lc <= 0.0){
        std::cerr << "WARNING: --band_pass lambda_c=" << lc
                  << " must be > 0 -- ignored" << std::endl;
        continue;
      }
      std::ostringstream ss; ss << lc;
      std::string lbl = "bp" + ss.str();          // self-documenting, e.g. bp0.1
      weight_specs.push_back({lbl, false, -lc});  // NEGATIVE mgap_val = window flag
      std::cout << GridLogMessage << "--band_pass: track '" << lbl
                << "' Gaussian window lambda_c=" << lc << std::endl;
    }
  }
  // m_gap_auto: min_n|mu_n| computed from the --evals text file (mass_f fallback).
  // Must be known before the topo accumulation loop — embedded records only
  // supply per-mode mu_n, not a pre-loop global min.
  double m_gap_auto;
  if(!evals.empty()){
    m_gap_auto = std::fabs(evals[0]);
    for(auto& v : evals) if(std::fabs(v) < m_gap_auto) m_gap_auto = std::fabs(v);
    std::cout << GridLogMessage << "m_gap_auto = min|mu_n| = " << m_gap_auto
              << " (from --evals)" << std::endl;
  } else {
    m_gap_auto = mass_f;
    std::cout << GridLogMessage << "m_gap_auto = mass_f = " << m_gap_auto
              << " (--evals not provided; supply for accurate mgap weight)" << std::endl;
  }

  // Identify the n_topo near-zero (topological) modes by SMALLEST |mu_n|,
  // independent of the order in the --evals file.  The eigenvalue file is NOT
  // guaranteed |mu_n|-sorted (e.g. conf 795 lists the zero mode at index 1,
  // not 0), so an index-position cut would exclude the wrong mode.  is_bulk[c]
  // = 0 for the n_topo smallest-|mu_n| modes (topological), 1 otherwise (bulk).
  // Used to build the direct bulk remainder B_bulk = sum_{bulk} chi_n^B.
  std::vector<char> is_bulk(evals.size(), 1);
  if(n_topo > 0 && !evals.empty()){
    std::vector<int> idx(evals.size());
    for(size_t i=0;i<idx.size();i++) idx[i]=(int)i;
    std::sort(idx.begin(), idx.end(),
              [&](int a,int b){ return std::fabs(evals[a]) < std::fabs(evals[b]); });
    int k = std::min(n_topo, (int)idx.size());
    std::cout << GridLogMessage << "--n_topo: excluding " << k
              << " smallest-|mu_n| mode(s) from B_bulk:";
    for(int j=0;j<k;j++){ is_bulk[idx[j]] = 0;
      std::cout << " c=" << idx[j] << "(mu=" << evals[idx[j]] << ")"; }
    std::cout << std::endl;
  }
  // Print active weight tracks
  {
    std::cout << GridLogMessage << "Active weight tracks (" << weight_specs.size() << "):";
    for(auto& ws : weight_specs){
      if(ws.is_sign)             std::cout << "  " << ws.label << "(sign)";
      else if(ws.mgap_val == 0)  std::cout << "  " << ws.label << "(mgap_auto)";
      else                       std::cout << "  " << ws.label << "(" << ws.mgap_val << ")";
    }
    std::cout << std::endl;
  }

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
    std::cout << GridLogMessage << "Reading file1: "<<file_list1[c]<<std::endl;
    readFile(data1[c],file_list1[c]);
    std::cout << GridLogMessage << "Sum file1["<<c<<"] Q=" <<real(TensorRemove(sum(data1[c])))<<std::endl;
  }
  if(take_diff){
    for(int c=0;c<data1.size()-1;c++)
      data1[c] = data1[c+1] - data1[c];
    LatticeComplexD tmp(data1[0].Grid());
    std::cout << GridLogMessage << "Reading file1: "<<file_list1.back()<<std::endl;
    readFile(tmp,file_list1.back());
    std::cout << GridLogMessage << "Sum file1[last] Q="<<real(TensorRemove(sum(tmp)))<<std::endl;
    data1.back() = tmp - data1.back();
  }

  // 5D fields: data2
  //   data2[c]: plain sum over s (existing behaviour, used for visualisation)
  //   q_eps, q_bdy, q_mid: three q_top formulas accumulated over eigenvectors
  bool compute_topo = (Ls > 0) && !file_list2.empty();
  // Dynamic topo accumulator vectors — one entry per active weight track.
  // q_eps = formula B (eps_code chirality), q_bdy = formula C (boundary), q_mid = formula A (midpoint).
  int ntracks = (int)weight_specs.size();
  std::vector<LatticeComplexD> q_eps_all, q_bdy_all, q_mid_all;
  // Sign-correct legacy estimator and Banks-Casher locality density.
  // q_naive   : reconstruction of legacy sp_sum (Compute_DWF_SpectralFlow.cc archive),
  //             with cleanly-named per-mode pieces — see §1.3 of
  //             sp_sum_vs_qB_bulk_analysis.tex for the derivation.
  // Sigma_low : low-mode contribution to local scalar condensate density,
  //             sum_n [m_f / (lambda^(0)_n)^2 + m_f^2] rho_n(x). Banks–Casher
  //             collapses the Lorentzian to pi*delta(lambda^(0)) as m_f -> 0,
  //             yielding the "where does Sigma live" map (cf. §4.1 of the doc).
  // Both use the same eps_s = (s>=Ls/2)?+1:-1 convention as q_eps_all/q_bdy_all
  // for sign-consistency with the master-table PCFs of q_B/q_C.
  LatticeComplexD q_naive(grid);   q_naive   = Zero();
  LatticeComplexD Sigma_low(grid); Sigma_low = Zero();
  // B_bulk(x) = sum_{n in bulk} chi_n^B(x): the DIRECT (exact) bulk remainder,
  // excluding the n_topo leading near-zero modes.  This is the formal B of
  // eq:Btrunc computed from the eigenvectors without the q_naive - q_B^mgap
  // approximation, so it carries NO topological pedestal (int B_bulk = 0 to
  // the +-lambda pairing accuracy).  Requires --n_topo = |Q|.
  LatticeComplexD B_bulk(grid);    B_bulk    = Zero();
  // q_B^unit(x) = sum_{n=1}^{N_low} chi_n^B(x): the unit-weight bulk-grading
  // sum over ALL retained modes (zero modes + bulk), no m_gap/|mu_n|
  // suppression and no Banks-Casher term.  This is the "unit-weight bulk
  // estimator" discussed in topo_charge.tex §sec:qBnaive: same integral -> Q_ov as the
  // gap-weighted q_B (zero mode contributes +-1, bulk pairs cancel via the
  // sum rule for any common weight), but the bulk-mode topological structure
  // (the I-Ibar pattern of the companion numerics) is preserved at full
  // amplitude rather than suppressed by m_gap/|mu_n|.  Mass-independent.
  LatticeComplexD q_B_unit(grid); q_B_unit = Zero();
  if(compute_topo){
    q_eps_all.reserve(ntracks); q_bdy_all.reserve(ntracks); q_mid_all.reserve(ntracks);
    for(int t = 0; t < ntracks; t++){
      q_eps_all.emplace_back(grid); q_eps_all.back() = Zero();
      q_bdy_all.emplace_back(grid); q_bdy_all.back() = Zero();
      q_mid_all.emplace_back(grid); q_mid_all.back() = Zero();
    }
  }

  // Cache mgap and sign weight-track indices.  Used by the in-loop B_k sweep
  // (cumulative spectral sum diagnostic of sec:bk_sweep) and the post-loop
  // alpha-sweep / smear-sweep blocks.
  int mgap_idx = -1, sign_idx = -1;
  for(int t = 0; t < (int)weight_specs.size(); t++){
    if(weight_specs[t].is_sign)              sign_idx = t;
    else if(weight_specs[t].label == "mgap") mgap_idx = t;
  }
  // Bk-sweep enable flag.  Default OFF (extra memory + IO); enable with
  // --bk_sweep on the command line.  When ON, the per-mode loop checkpoints the
  // DIRECT cumulative bulk remainder B_bulk_k = sum_{bulk n<=k} chi_n^B (built
  // straight from the eigenvectors) after each mode k, so a post-loop block can
  // emit PCF/IP of B_k vs each gluonic data1[i] as a function of the mode cutoff.
  // (Earlier this used the proxy B_k = q_naive_k - q_B^mgap_k; B is now the same
  // direct field used everywhere else -- see B_bulk and sp_sum sec:Baxis_bk.)
  // Memory cost: N_conv * V * 16B ~ 16 * 16 MB = 256 MB on 32^4 with N_conv=16.
  bool do_bk_sweep =  GridCmdOptionExists(argv,argv+argc,"--bk_sweep");
  // Per-mode spectral diagnostics.  Default OFF; enable with --per_mode_out.
  // When ON, the per-mode loop records, for each retained mode n:
  //   mu_n, int chi_n^B d^4x, int rho_n, int rho_n^2
  // and a post-loop block writes them (with IPR_n = V4 * int rho^2 / (int rho)^2)
  // to <topo_out:{def}->permode>.dat.  Tests 2 & 6 of sp_sum App. window_tests:
  // |int chi_n^B| vs |mu_n| confirms the R5-symmetry suppression; IPR_n vs |mu_n|
  // locates the localization (mobility) edge for the --band_pass lambda_c.
  bool do_per_mode = GridCmdOptionExists(argv,argv+argc,"--per_mode_out");
  std::vector<double> pm_mu, pm_chiI, pm_rhoI, pm_rho2;
  std::vector<LatticeComplexD> B_bulk_partial;   // cumulative direct B_bulk(k) for bk_sweep

  // Eigenvalues embedded inside each SCIDAC evec_density file (H_DWF_EvalRecord).
  // Populated during data2 loading; takes priority over --evals when available.
  std::vector<double> evals_embedded;

  std::vector<LatticeComplexD> data2(file_list2.size()-take_diff,grid);
  for(int c=0;c<data2.size();c++) {
    std::cout << GridLogMessage << "Reading file2: "<<file_list2[c]<<std::endl;
    LatticeComplexD tmp(gridF);
    double emb_eval = 0.0; int emb_n = -1;
    readFileRecord(tmp, file_list2[c], emb_eval, emb_n);
    if(emb_n >= 0)   // valid embedded record
      evals_embedded.push_back(emb_eval);
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
    std::cout << GridLogMessage << "Sum file2["<<c<<"] Q="<<real(TensorRemove(sum(data2[c])))<<std::endl;

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
    //   Near-zero (topological) modes: |mu_n| ~ m_gap = m_f + m_res.
    //   Bulk (+/-mu) pairs: contributions cancel in the sum.
    //
    // Both weight modes (sign and mgap) are accumulated simultaneously — see loop body.
    // q_A midpoint is computed with BOTH weight tracks (see Formula A note below:
    // q_A is empirically a bulk-mode midpoint residue, not an index estimator).
    //
    // -----------------------------------------------------------------------
    // FORMULA B:  eps_code(s)-chirality density   (Bulk form)
    // -----------------------------------------------------------------------
    //   q_B(x) = -sum_n (m_gap/mu_n) * sum_s eps_code(s) * rho_n(x,s)
    //
    //   Derivation:
    //     q^exact = -sum_n (m_f/mu_n) chi_n^B,   chi_n^B = sum_s Gamma5_Blum(s)|psi_n|^2
    //     eps_code(s) = -Gamma5_Blum(s),  so chi_n^B = -sum_s eps_code(s)|psi_n|^2
    //     replacing m_f/mu_n -> m_gap/mu_n:
    //       q_B = -sum_n (m_gap/mu_n)*(-sum_s eps_code(s)|psi_n|^2)
    //           = -sum_n (m_gap/mu_n) * sum_s eps_code(s)|psi_n|^2  (with overall minus)
    //   Near-zero modes localized at one wall contribute ~+-1; bulk symmetric modes
    //   cancel because eps_code sums to zero over a uniformly distributed mode.
    //
    // -----------------------------------------------------------------------
    // FORMULA C:  boundary projection density      (Boundary form)
    // -----------------------------------------------------------------------
    //   q_C(x) = -sum_n (m_gap/mu_n) * [rho_n(x,Ls-1) - rho_n(x,0)]
    //
    //   Derivation:
    //     q_top(x) = -m_f * tr[gamma5 * S^{4D}(x,x)]
    //              = -m_f * sum_n (1/mu_n) * [psi_n†(x,Ls-1) P_R psi_n(x,Ls-1)
    //                                        + psi_n†(x,0)   P_L psi_n(x,0)]
    //   where P_R = (1+gamma5)/2, P_L = (1-gamma5)/2  with the STANDARD 4D gamma5.
    //   Replacing m_f/mu_n -> m_gap/mu_n and spinor projections by scalar proxy:
    //     psi†(x,Ls-1) P_R psi(x,Ls-1) ~ rho(x,Ls-1)  [right-wall mode: P_R ~ 1]
    //     psi†(x,0)    P_L psi(x,0)    ~ rho(x,0)      [left-wall mode:  P_L ~ 1]
    //   The relative minus sign in (rho(Ls-1) - rho(0)) encodes the eps_code
    //   convention (eps(Ls-1)=+1, eps(0)=-1).
    //
    // -----------------------------------------------------------------------
    // FORMULA A:  midpoint density  (Midpoint form; BOTH weight tracks computed)
    // -----------------------------------------------------------------------
    //   q_A(x) = -sum_n w(mu_n) * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
    //
    //   Empirical behaviour (conf 702, Ls=12, low-mode truncation N_conv=16):
    //     - Q_evec = integral of q_A over all x is ~ 1e-8 for BOTH weight tracks.
    //       I.e. q_A carries NO index content: topological near-zero modes have
    //       amplitude ~ e^{-alpha*Ls/2} at the midplane and are invisible to
    //       this estimator under truncation.
    //     - The UV bulk modes that would normally dominate J_5q at the midplane
    //       (anomaly inflow carriers) are NOT in the retained low-mode set, so
    //       the midplane current is "severed" under deflation.
    //     - What remains is a bulk-pair midpoint RESIDUE: local-chirality-driven
    //       asymmetry between s=Ls/2 and s=Ls/2-1 from retained bulk modes.  This
    //       residue has nontrivial spatial structure (can correlate with
    //       Wilson-flowed q_gluon when the gauge is smoothed) but integrates to
    //       zero. See topo_charge.tex App. B.3 (bulk modes vs near-zero modes
    //       in J_5q vs j_5).
    //
    //   Historical note: an earlier version of this comment claimed
    //     m_f/mu_n ~ e^{+alpha*Ls} compensates the midpoint suppression.
    //   That claim is WRONG at finite Ls with physical m_f: near-zero
    //   |mu_n| ~= m_gap = m_f + m_res (not m_res), so m_gap/mu_n ~= +-1,
    //   not e^{+alpha*Ls}.  Numerical check: Q_evec(q_A_sign) = Q_evec(q_A_mgap)
    //   ~ 1e-8 (identical up to rounding), confirming no compensation occurs.
    // ===========================================================================
    if(compute_topo) {
      // Eigenvalue of H_DWF = gamma5*R5*D_DWF (includes m_f).
      // Priority: (1) embedded H_DWF_EvalRecord, (2) --evals, (3) 0.0 (warning).
      double mu_n;
      if(c < (int)evals_embedded.size())
        mu_n = evals_embedded[c];
      else if(c < (int)evals.size())
        mu_n = evals[c];
      else {
        mu_n = 0.0;
        std::cout << GridLogMessage << "WARNING: no eigenvalue for mode " << c
                  << " — set mu_n=0, mode will not contribute" << std::endl;
      }

      // Compute weight for each active track:
      //   is_sign = true  => w = sign(mu_n)                            (exact ±1 or 0)
      //   is_sign = false, mgap_val = 0 => w = m_gap_auto / |mu_n|   (auto = min|mu_n|)
      //   is_sign = false, mgap_val ≠ 0 => w = mgap_val / |mu_n|      (user-supplied literal)
      //   (|mu_n| is the default; --signed_mgap uses mu_n, the legacy sign-blind weight)
      // q_A midpoint note: both weight tracks give integral Q_evec(q_A) ~ 0
      //   because topological modes have e^{-alpha*Ls/2}-suppressed midplane
      //   amplitude AND the UV bulk modes that dominate J_5q at the midplane
      //   are excluded by truncation.  q_A carries only a bulk-pair residue.
      //   (See FORMULA A header comment above.)
      std::vector<double> track_w(ntracks, 0.0);
      for(int t = 0; t < ntracks; t++){
        const WeightSpec& ws = weight_specs[t];
        if(ws.is_sign)
          track_w[t] = (mu_n > 0.0) ? 1.0 : (mu_n < 0.0) ? -1.0 : 0.0;
        else if(ws.mgap_val < 0.0) {
          // Band-pass Gaussian window (--band_pass): lambda_c = -mgap_val.
          // w(mu)=exp(-mu^2/(2 lc^2)) is even and ->1 at mu=0, so the integer
          // Q_ov is preserved while heavy modes roll off at the mobility edge lc.
          double lc = -ws.mgap_val;
          track_w[t] = std::exp(-(mu_n*mu_n) / (2.0*lc*lc));
        }
        else {
          double mg = (ws.mgap_val == 0.0) ? m_gap_auto : ws.mgap_val;
          // Default: m_gap/|mu_n| (signed-index-faithful; the chirality chi_n^B
          // carries the sign).  --signed_mgap restores the legacy m_gap/mu_n,
          // which is sign-blind to Q (see header comment at signed_mgap decl).
          double denom = signed_mgap ? mu_n : std::fabs(mu_n);
          track_w[t] = (mu_n != 0.0) ? (mg / denom) : 0.0;
        }
      }

      // --- Formula B: eps_code chirality sum ---
      // q_B(x) += -w * sum_s eps_code(s) * rho_n(x,s)
      {
        LatticeComplexD eps_slice(grid);
        for(int s = 0; s < Ls; s++){
          ExtractSlice(eps_slice, tmp, s, 0);
          double eps_s = (s >= Ls/2) ? 1.0 : -1.0;
          for(int t = 0; t < ntracks; t++)
            q_eps_all[t] = q_eps_all[t] - (track_w[t] * eps_s) * eps_slice;
        }
      }

      // --- Formula C: boundary projection ---
      // q_C(x) += -w * [rho_n(x,Ls-1) - rho_n(x,0)]
      {
        LatticeComplexD bdy_s0(grid), bdy_sLs(grid);
        ExtractSlice(bdy_s0,  tmp, 0,    0);
        ExtractSlice(bdy_sLs, tmp, Ls-1, 0);
        LatticeComplexD bdy_diff = bdy_sLs - bdy_s0;
        for(int t = 0; t < ntracks; t++)
          q_bdy_all[t] = q_bdy_all[t] - track_w[t] * bdy_diff;
      }

      // --- Formula A: midpoint density ---
      // q_A(x) += -w * 0.5 * [rho_n(x,Ls/2) - rho_n(x,Ls/2-1)]
      if(Ls >= 2){
        LatticeComplexD mid_lo(grid), mid_hi(grid);
        ExtractSlice(mid_lo, tmp, Ls/2-1, 0);
        ExtractSlice(mid_hi, tmp, Ls/2,   0);
        LatticeComplexD mid_diff = mid_hi - mid_lo;
        for(int t = 0; t < ntracks; t++)
          q_mid_all[t] = q_mid_all[t] - track_w[t] * 0.5 * mid_diff;
      }

      // --- q_naive (legacy sp_sum) and Sigma_low (Banks-Casher locality) ---
      //   q_naive(x)   = +sum_n chi_n^B(x)
      //                  + (1/2) sum_n sgn(lambda^H_n) * lambda^(0)_n * rho_n(x)
      // where chi_n^B(x) = sum_s Gamma_5(s) |u_n(x,s)|^2  (Gamma_5 convention,
      // matching eq:qnaive_split in sp_sum_vs_qB_bulk_analysis.tex).
      // The code achieves the + sign via "-eps_s * |u|^2" with eps_s = eps_code(s)
      // = -Gamma_5(s), i.e. the same double-negative as the archived sp_sum block
      // in Compute_DWF_SpectralFlow.cc.
      //   Sigma_low(x) = sum_n [m_f / ((lambda^(0)_n)^2 + m_f^2)] * rho_n(x)
      // with lambda^(0)_n = sqrt((lambda^H_n)^2 - m_f^2) the kinetic eigenvalue
      // (App. B.4). No --weights track loop: q_naive and Sigma_low are
      // single fields, accumulated unconditionally for every mode.
      {
        double lambda_H     = mu_n;
        double lambda_0_sq  = lambda_H*lambda_H - mass_f*mass_f;
        double lambda_0     = (lambda_0_sq > 0.0) ? std::sqrt(lambda_0_sq) : 0.0;
        double sgn_lambda   = (lambda_H > 0.0) ? 1.0
                            : (lambda_H < 0.0) ? -1.0 : 0.0;
        // Banks–Casher Lorentzian uses bc_mass (independent of mass_f) so the
        // §4.1 scalar-foil diagnostic remains well-defined for zero-mass evecs.
        double bc_weight    = bc_mass / (lambda_0_sq + bc_mass*bc_mass);

        // rho_n(x) = sum_s |u_n(x,s)|^2  (scalar 4D density)
        LatticeComplexD rho_n(grid); rho_n = Zero();
        {
          LatticeComplexD slice(grid);
          for(int s = 0; s < Ls; s++){
            ExtractSlice(slice, tmp, s, 0);
            rho_n = rho_n + slice;
          }
        }

        // q_naive first sum: same eps_s convention as q_B/q_C above.
        // q_eps_all uses "-(track_w * eps_s)" pattern; for q_naive the
        // per-mode coefficient is unity (no track weight), so it picks up
        // "-eps_s * |u(x,s)|^2" summed over s = +chi^B_doc per mode (eps_s
        // is -Gamma_5_doc). Net: q_naive's first sum has the same sign
        // convention as q_B^sign in this codebase.
        {
          LatticeComplexD chi_n(grid); chi_n = Zero();
          LatticeComplexD eps_slice(grid);
          for(int s = 0; s < Ls; s++){
            ExtractSlice(eps_slice, tmp, s, 0);
            double eps_s = (s >= Ls/2) ? 1.0 : -1.0;
            chi_n = chi_n - eps_s * eps_slice;   // = +chi^B_doc for mode c
          }
          q_naive = q_naive + chi_n + (0.5 * sgn_lambda * lambda_0) * rho_n;
          // Naive bulk estimator: q_B^naive = sum_n chi_n^B summed over ALL
          // retained modes (zero modes + bulk), unweighted.  Same chi_n as
          // q_naive's first sum.  See topo_charge.tex sec:qBnaive.
          q_B_unit = q_B_unit + chi_n;
          // Direct bulk remainder: the same per-mode chirality density, but
          // accumulated only for the bulk modes (skip the n_topo smallest-|mu|
          // near-zero modes, identified above by |mu_n| -- robust to an
          // unsorted --evals file).  If no --evals (is_bulk empty), fall back
          // to the index cut c>=n_topo.  Gives B_bulk = sum_{n in bulk} chi_n^B
          // with no pedestal, the exact form of eq:Btrunc.
          bool mode_is_bulk = (c < (int)is_bulk.size()) ? (bool)is_bulk[c]
                                                        : (c >= n_topo);
          if(mode_is_bulk) B_bulk = B_bulk + chi_n;

          // Per-mode spectral diagnostics (--per_mode_out): record this mode's
          // eigenvalue and the integrated densities for tests 2 & 6.  chi_n is
          // +chi_n^B (doc convention); rho_n is the scalar 4D density.
          if(do_per_mode){
            pm_mu.push_back(mu_n);
            pm_chiI.push_back(real(TensorRemove(sum(chi_n))));
            pm_rhoI.push_back(real(TensorRemove(sum(rho_n))));
            pm_rho2.push_back(real(TensorRemove(sum(rho_n*rho_n))));
          }
        }

        // Sigma_low: Banks-Casher Lorentzian-weighted scalar density.
        Sigma_low = Sigma_low + bc_weight * rho_n;
      }

      /****** Bk sweep checkpoint (sec:bk_sweep) ****************************/
      // If --bk_sweep, snapshot the running DIRECT bulk remainder B_bulk AFTER
      // mode c has been added (B_bulk already accumulates sum_{bulk n<=c} chi_n^B
      // straight from the eigenvectors).  A post-loop block emits PCF/IP of
      // B_k = B_bulk(k) vs each gluonic data1[i].  We don't emit here because
      // tau_wf/conf_id/get_td_tau are parsed later in main.
      if(do_bk_sweep){
        B_bulk_partial.emplace_back(grid);
        B_bulk_partial.back() = B_bulk;
      }
      /**********************************************************************/

      std::cout << "TopoContrib evec=" << c << " mu_n=" << mu_n;
      for(int t = 0; t < ntracks; t++)
        std::cout << " w_" << weight_specs[t].label << "=" << track_w[t];
      for(int t = 0; t < ntracks; t++)
        std::cout << " Q_B_" << weight_specs[t].label
                  << "=" << real(TensorRemove(sum(q_eps_all[t])));
      std::cout << std::endl;
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
  // q_A, q_B, q_C all use m_gap/mu_n weight (requires --evals for m_gap).
  bool have_evals = !evals.empty() || !evals_embedded.empty();
  if(compute_topo && have_evals && do_topo_out){
    // Lambda: substitute {def} placeholder in topo_out template string.
    auto fill_def = [](std::string tmpl, const std::string& d) -> std::string {
      auto pos = tmpl.find("{def}");
      if(pos != std::string::npos) tmpl.replace(pos, 5, d);
      return tmpl;
    };

    for(int t = 0; t < ntracks; t++){
      const std::string& lbl = weight_specs[t].label;
      writeFile(q_eps_all[t], fill_def(topo_out,"B_"+lbl));
      std::cout << "Wrote q_B_" << lbl << " -> " << fill_def(topo_out,"B_"+lbl)
                << "  Q=" << real(TensorRemove(sum(q_eps_all[t])));
      if(!weight_specs[t].is_sign){
        if(weight_specs[t].mgap_val == 0.0) std::cout << "  m_gap_auto=" << m_gap_auto;
        else                                std::cout << "  mgap=" << weight_specs[t].mgap_val;
      }
      std::cout << std::endl;

      writeFile(q_bdy_all[t], fill_def(topo_out,"C_"+lbl));
      std::cout << "Wrote q_C_" << lbl << " -> " << fill_def(topo_out,"C_"+lbl)
                << "  Q=" << real(TensorRemove(sum(q_bdy_all[t]))) << std::endl;

      writeFile(q_mid_all[t], fill_def(topo_out,"A_"+lbl));
      std::cout << "Wrote q_A_" << lbl << " -> " << fill_def(topo_out,"A_"+lbl)
                << "  Q=" << real(TensorRemove(sum(q_mid_all[t])));
      if(weight_specs[t].is_sign)          std::cout << "  (note: q_A_sign ~ 0 for Ls=48)";
      std::cout << std::endl;
    }

    // q_naive (legacy sp_sum, sign-correct) and Sigma_low (Banks-Casher locality).
    // Single fields, no weight-track loop.
    writeFile(q_naive,   fill_def(topo_out,"naive"));
    std::cout << "Wrote q_naive   -> " << fill_def(topo_out,"naive")
              << "  Q=" << real(TensorRemove(sum(q_naive))) << std::endl;
    writeFile(Sigma_low, fill_def(topo_out,"Sigma"));
    std::cout << "Wrote Sigma_low -> " << fill_def(topo_out,"Sigma")
              << "  S=" << real(TensorRemove(sum(Sigma_low))) << std::endl;
    writeFile(B_bulk,    fill_def(topo_out,"Bbulk"));
    std::cout << "Wrote B_bulk    -> " << fill_def(topo_out,"Bbulk")
              << "  intB=" << real(TensorRemove(sum(B_bulk)))
              << "  (n_topo=" << n_topo << ")" << std::endl;
    writeFile(q_B_unit, fill_def(topo_out,"Bunit"));
    std::cout << "Wrote q_B_unit -> " << fill_def(topo_out,"Bunit")
              << "  Q=" << real(TensorRemove(sum(q_B_unit))) << std::endl;

    // Per-mode spectral diagnostics (--per_mode_out): write mu_n, int chi_n^B,
    // int rho_n, int rho_n^2, and IPR_n = V4 * int rho^2 / (int rho)^2 to a text
    // file <topo_out:{def}->permode>.dat (one row per retained mode).
    if(do_per_mode && !pm_mu.empty()){
      double V4 = (double)grid->gSites();
      std::string pmpath = fill_def(topo_out, "permode") + ".dat";
      std::ofstream pmf(pmpath);
      pmf << "# n  mu_n  int_chiB  int_rho  int_rho2  IPR\n";
      for(size_t n = 0; n < pm_mu.size(); n++){
        double ipr = (pm_rhoI[n] != 0.0)
                   ? V4 * pm_rho2[n] / (pm_rhoI[n]*pm_rhoI[n]) : 0.0;
        pmf << n << " " << pm_mu[n] << " " << pm_chiI[n] << " "
            << pm_rhoI[n] << " " << pm_rho2[n] << " " << ipr << "\n";
        std::cout << "PerMode: n=" << n << " mu=" << pm_mu[n]
                  << " intChiB=" << pm_chiI[n] << " IPR=" << ipr << std::endl;
      }
      std::cout << "Wrote per-mode spectral data -> " << pmpath << std::endl;
    }
  }
  /****** IP & Corr of each fermion TCD definition vs gluonic TCD (--topo_compare) *****/
  // Requires: --files2 + Ls (so q_eps/q_bdy/q_mid accumulate), have_evals (for weights),
  //           --files1   (one file per gluonic flow time, e.g. TD_tau=0,4,16),
  //           --conf_id  (integer config number, written to output for notebook parsing).
  // Output per line (pure numeric): TD_tau  tau  conf  Corr  IP
  //   where TD_tau is the index i of files1[i], tau comes from --conf_id context,
  //   and the label "Topo_PCF: q_X" is printed to stderr for diagnostics only.
  // q_A, q_B, q_C all use m_gap/mu_n weight.
  // The 3 output blocks (one per q_def) are labelled by a comment line "# q_X"
  // so the shell can split them into 4 per-def files via grep.
  int conf_id = -1;
  if( GridCmdOptionExists(argv,argv+argc,"--conf_id") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--conf_id");
    GridCmdOptionInt(arg, conf_id);
  }
  // --topo_log <file>: legacy scratch file for machine-readable lines (kept for
  // backward compatibility).  Superseded by --data_dir which writes named output
  // files directly.  When both are absent, falls back to stdout.
  std::ofstream topo_log_f;
  std::ostream* topo_log = &std::cout;
  if( GridCmdOptionExists(argv,argv+argc,"--topo_log") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--topo_log");
    topo_log_f.open(arg, std::ios::trunc);
    if(!topo_log_f)
      std::cerr << "WARNING: cannot open --topo_log file: " << arg
                << " — falling back to stdout" << std::endl;
    else
      topo_log = &topo_log_f;
  }
  // --data_dir <dir>: write corr_ip_q_*.dat, comp_ref_q_*.dat, corr_ip_stoch.dat
  // directly in append mode — no scratch file, no awk post-processing needed.
  // --tau_wf  <val>: Wilson flow time written as the 'tau' column in output rows.
  // --td_taus <v0,v1,...>: actual TD_tau values for data1[0], data1[1], ...
  //   Default when absent: use the file index directly (0, 1, 2, ...).
  std::string data_dir;
  int tau_wf = -1;
  std::vector<int> td_taus_vec;
  if( GridCmdOptionExists(argv,argv+argc,"--data_dir") ){
    data_dir = GridCmdOptionPayload(argv,argv+argc,"--data_dir");
    std::cout << GridLogMessage << "--data_dir: " << data_dir << std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--tau_wf") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--tau_wf");
    GridCmdOptionInt(arg, tau_wf);
    std::cout << GridLogMessage << "--tau_wf: " << tau_wf << std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--td_taus") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--td_taus");
    std::vector<std::string> td_strs;
    GridCmdOptionCSL(arg, td_strs);
    for(auto& s: td_strs) td_taus_vec.push_back(std::stoi(s));
    std::cout << GridLogMessage << "--td_taus: " << td_taus_vec << std::endl;
  }
  // Map data1 index i -> actual TD_tau value
  auto get_td_tau = [&](int i) -> int {
    return (i < (int)td_taus_vec.size()) ? td_taus_vec[i] : i;
  };
  int topo_compare = GridCmdOptionExists(argv,argv+argc,"--topo_compare");
  if(compute_topo && have_evals && topo_compare && !data1.empty()){
    typedef typename PeriodicGimplR::ComplexField ComplexField;
    // data1[i] = gluonic TCD at flow time index i  (one file per TD_tau)
    struct QDef { std::string name; LatticeComplexD* field; };
    std::vector<QDef> qdefs;
    for(int t = 0; t < ntracks; t++){
      const std::string& lbl = weight_specs[t].label;
      qdefs.push_back({"q_A_"+lbl, &q_mid_all[t]});
      qdefs.push_back({"q_B_"+lbl, &q_eps_all[t]});
      qdefs.push_back({"q_C_"+lbl, &q_bdy_all[t]});
    }
    qdefs.push_back({"q_naive",   &q_naive});
    qdefs.push_back({"Sigma_low", &Sigma_low});
    qdefs.push_back({"B_bulk",    &B_bulk});
    qdefs.push_back({"q_B_unit", &q_B_unit});
    for(auto& qd : qdefs){
      if(data_dir.empty()) *topo_log << "# " << qd.name << std::endl;
      for(int i=0; i<(int)data1.size(); i++){
        LatticeComplex X(grid), Y(grid), one(grid); one = ComplexField::scalar_type(1.0,0.0);
        ComplexD avg1 = TensorRemove(sum(data1[i]))/RealD(grid->gSites());
        ComplexD avg2 = TensorRemove(sum(*qd.field))/RealD(grid->gSites());
        X = data1[i] - avg1*one;
        Y = *qd.field  - avg2*one;
        double corr = TensorRemove(sum(X*Y)).real()/sqrt(norm2(X)*norm2(Y));
        X = data1[i]; Y = *qd.field;
        double ip   = TensorRemove(innerProduct(X,Y)).real()/sqrt(norm2(X))/sqrt(norm2(Y));
        // Always print to stdout for monitoring (tau_wf identifies which eigenvector smearing)
        std::cout << "Topo PCF Corr: " << qd.name << " tau_wf=" << tau_wf
                  << " TD_tau=" << get_td_tau(i) << " conf=" << conf_id << " " << corr << std::endl;
        std::cout << "Topo PCF IP:   " << qd.name << " tau_wf=" << tau_wf
                  << " TD_tau=" << get_td_tau(i) << " conf=" << conf_id << " " << ip   << std::endl;
        if(!data_dir.empty()){
          // Direct append: TD_tau  tau_wf  conf  value
          std::ofstream of(data_dir+"/corr_ip_"+qd.name+".dat", std::ios::app);
          of << get_td_tau(i) << " " << tau_wf << " " << conf_id << " " << corr << "\n";
          of << get_td_tau(i) << " " << tau_wf << " " << conf_id << " " << ip   << "\n";
        } else {
          *topo_log << "Topo PCF Corr: " << i << " " << conf_id << " " << corr << std::endl;
          *topo_log << "Topo PCF IP:   " << i << " " << conf_id << " " << ip   << std::endl;
        }
      }
    }

    /****** Alpha-sweep: PCF/IP of q_B^mgap + alpha*B vs gluonic q_top ************/
    // Implements the §sec:alpha_sweep diagnostic of sp_sum_vs_qB_bulk_analysis.tex:
    // sweep alpha to identify the alpha^* that best matches each gluonic
    // reference (one alpha^* per TD_tau).  Since ∫B = 0, the integrated charge
    // ∫q_mix = Q_top is preserved exactly for every alpha.
    // Computes once per (alpha, TD_tau) and appends to alpha_sweep.dat.
    int sign_idx = -1, mgap_idx = -1;
    for(int t = 0; t < ntracks; t++){
      if(weight_specs[t].is_sign)            sign_idx = t;
      else if(weight_specs[t].label == "mgap") mgap_idx = t;
    }
    if(sign_idx >= 0 && mgap_idx >= 0){
      // B(x) = B_bulk(x) = sum_{n in bulk} chi_n^B, the DIRECT bulk remainder
      // (eq:Btrunc) built straight from the eigenvectors -- NOT the old proxy
      // q_naive - q_B^sign.  Direct B_bulk has int B = 0 to machine precision,
      // so int q_mix = int q_B^mgap = Q_top holds exactly for every alpha.
      LatticeComplexD B_field(grid);
      B_field = B_bulk;
      // Alpha grid (default; --alpha_sweep "0,0.1,..." overrides)
      std::vector<double> alphas = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0};
      if(GridCmdOptionExists(argv,argv+argc,"--alpha_sweep")){
        std::string a_arg = GridCmdOptionPayload(argv,argv+argc,"--alpha_sweep");
        std::vector<std::string> a_strs;
        GridCmdOptionCSL(a_arg, a_strs);
        alphas.clear();
        for(auto& s : a_strs) alphas.push_back(std::stod(s));
      }
      std::cout << GridLogMessage << "Alpha sweep: " << alphas.size() << " values" << std::endl;
      LatticeComplexD q_mix(grid), one_alpha(grid);
      one_alpha = ComplexField::scalar_type(1.0, 0.0);
      for(double alpha : alphas){
        q_mix = q_eps_all[mgap_idx] + alpha * B_field;
        for(int i = 0; i < (int)data1.size(); i++){
          ComplexD avg1 = TensorRemove(sum(data1[i]))/RealD(grid->gSites());
          ComplexD avg2 = TensorRemove(sum(q_mix))   /RealD(grid->gSites());
          LatticeComplexD X(grid), Y(grid);
          X = data1[i] - avg1*one_alpha;
          Y = q_mix    - avg2*one_alpha;
          double corr = real(TensorRemove(sum(X*Y))) / std::sqrt(norm2(X)*norm2(Y));
          double ip   = real(TensorRemove(innerProduct(data1[i], q_mix)))
                      / std::sqrt(norm2(data1[i])) / std::sqrt(norm2(q_mix));
          std::cout << "AlphaSweep: alpha=" << alpha
                    << " tau_wf=" << tau_wf
                    << " TD_tau=" << get_td_tau(i)
                    << " conf=" << conf_id
                    << " Corr=" << corr << " IP=" << ip << std::endl;
          if(!data_dir.empty()){
            std::ofstream of(data_dir+"/alpha_sweep.dat", std::ios::app);
            of << alpha << " " << tau_wf << " " << get_td_tau(i) << " "
               << conf_id << " " << corr << " " << ip << "\n";
          }
        }
      }
    } else {
      std::cout << GridLogMessage
                << "Alpha sweep: skipped (need both sign and mgap weight tracks)"
                << std::endl;
    }

    /****** Bk sweep emission (sec:bk_sweep) ********************************/
    // Cumulative spectral sum: how concentrated in the lowest modes is the
    // chirality structure of B?  Use the direct cumulative bulk remainder
    // B_k = B_bulk(k) = sum_{bulk n<=k} chi_n^B saved during the per-mode loop,
    // and emit PCF/IP vs q_gluon at each TD_tau, for each k = 1..N_conv.
    // Output schema (bk_sweep.dat):  k tau_wf TD_tau conf Corr IP
    if(do_bk_sweep && !B_bulk_partial.empty() && !data1.empty()){
      LatticeComplexD oneL(grid); oneL = ComplexD(1.0, 0.0);
      int nk = (int)B_bulk_partial.size();
      std::cout << GridLogMessage << "Bk sweep: " << nk << " modes" << std::endl;
      for(int k = 0; k < nk; k++){
        LatticeComplexD B_k(grid);
        B_k = B_bulk_partial[k];
        for(int i = 0; i < (int)data1.size(); i++){
          ComplexD avg1 = TensorRemove(sum(data1[i])) / RealD(grid->gSites());
          ComplexD avgB = TensorRemove(sum(B_k))      / RealD(grid->gSites());
          LatticeComplexD X(grid), Y(grid);
          X = data1[i] - avg1 * oneL;
          Y = B_k       - avgB * oneL;
          double n2X = norm2(X), n2Y = norm2(Y);
          double corr = (n2X > 0 && n2Y > 0)
                      ? real(TensorRemove(sum(X*Y))) / std::sqrt(n2X * n2Y)
                      : 0.0;
          double n21 = norm2(data1[i]), n2B = norm2(B_k);
          double ip   = (n21 > 0 && n2B > 0)
                      ? real(TensorRemove(innerProduct(data1[i], B_k)))
                          / std::sqrt(n21 * n2B)
                      : 0.0;
          std::cout << "BkSweep: k=" << (k+1)
                    << " tau_wf=" << tau_wf
                    << " TD_tau=" << get_td_tau(i)
                    << " conf=" << conf_id
                    << " Corr=" << corr << " IP=" << ip << std::endl;
          if(!data_dir.empty()){
            std::ofstream of(data_dir+"/bk_sweep.dat", std::ios::app);
            of << (k+1) << " " << tau_wf << " " << get_td_tau(i) << " "
               << conf_id << " " << corr << " " << ip << "\n";
          }
        }
      }
      // Free partial-sum memory now that PCFs are written
      B_bulk_partial.clear();   B_bulk_partial.shrink_to_fit();
    }

    /****** Sigma sweep: density-level Gaussian smearing (sec:smear_sweep) ***/
    // Luchang's diagnostic: apply a 4D Gaussian density smear of width sigma
    // to BOTH a fermionic estimator AND each gluonic q_g(tau_WG) (and Luchang's
    // q_L if --comp_file is supplied), then measure PCF as a function of
    // sigma.  Decouples physical topology from signal-processing bandwidth:
    // Wilson flow preserves instantons; Gaussian density smearing is a pure
    // linear low-pass.
    //
    // Default sigma grid: {0.5, 1.0, 1.5, 2.0, 3.0, 5.0} lattice units;
    // override via --smear_sweep "0.3,0.6,..."
    // Output schema (smear_sweep.dat): kind name sigma tau_wf TD_tau conf Corr IP
    //   kind in {fermion, comp};  name is the operator label (e.g. q_B_mgap, comp_0).
    if(GridCmdOptionExists(argv,argv+argc,"--smear_sweep")){
      std::vector<double> sigmas = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
      std::string s_arg = GridCmdOptionPayload(argv,argv+argc,"--smear_sweep");
      if(!s_arg.empty()){
        std::vector<std::string> s_strs;
        GridCmdOptionCSL(s_arg, s_strs);
        if(!s_strs.empty()){
          sigmas.clear();
          for(auto& s : s_strs) sigmas.push_back(std::stod(s));
        }
      }
      std::cout << GridLogMessage << "Smear sweep: " << sigmas.size()
                << " sigma values" << std::endl;

      // Fermionic operators to include in the sweep.  Limit to a small set
      // to keep IO manageable; the headline test is q_B^mgap, with q_naive
      // as a sanity-check (its sigma-sweep should plateau quickly since it
      // already includes the bulk-mode chirality content).
      struct SsField { std::string name; LatticeComplexD* field; };
      std::vector<SsField> ferm_for_smear;
      if(mgap_idx >= 0)
        ferm_for_smear.push_back({"q_B_mgap", &q_eps_all[mgap_idx]});
      ferm_for_smear.push_back({"q_naive", &q_naive});

      LatticeComplexD oneL(grid); oneL = ComplexD(1.0, 0.0);

      for(double sigma : sigmas){
        // Smear each gluonic data1[i] once per sigma (cache for reuse below)
        std::vector<LatticeComplexD> data1_smr;
        data1_smr.reserve(data1.size());
        for(int i = 0; i < (int)data1.size(); i++)
          data1_smr.push_back(gaussian_smear(data1[i], sigma));

        // Fermion ops vs gluonic
        for(auto& ff : ferm_for_smear){
          LatticeComplexD ferm_smr = gaussian_smear(*ff.field, sigma);
          for(int i = 0; i < (int)data1.size(); i++){
            ComplexD avg1 = TensorRemove(sum(data1_smr[i])) / RealD(grid->gSites());
            ComplexD avgF = TensorRemove(sum(ferm_smr))    / RealD(grid->gSites());
            LatticeComplexD X(grid), Y(grid);
            X = data1_smr[i] - avg1 * oneL;
            Y = ferm_smr     - avgF * oneL;
            double n2X = norm2(X), n2Y = norm2(Y);
            double corr = (n2X > 0 && n2Y > 0)
                        ? real(TensorRemove(sum(X*Y))) / std::sqrt(n2X * n2Y)
                        : 0.0;
            double n21 = norm2(data1_smr[i]), n2F = norm2(ferm_smr);
            double ip   = (n21 > 0 && n2F > 0)
                        ? real(TensorRemove(innerProduct(data1_smr[i], ferm_smr)))
                            / std::sqrt(n21 * n2F)
                        : 0.0;
            std::cout << "SmearSweep: kind=fermion name=" << ff.name
                      << " sigma=" << sigma
                      << " tau_wf=" << tau_wf
                      << " TD_tau=" << get_td_tau(i)
                      << " conf=" << conf_id
                      << " Corr=" << corr << " IP=" << ip << std::endl;
            if(!data_dir.empty()){
              std::ofstream of(data_dir+"/smear_sweep.dat", std::ios::app);
              of << "fermion " << ff.name << " " << sigma << " "
                 << tau_wf << " " << get_td_tau(i) << " " << conf_id << " "
                 << corr << " " << ip << "\n";
            }
          }
        }
        // (Luchang/comp-file vs gluonic at this sigma is handled in the
        // separate Luchang block below if comp_files are present; that block
        // detects --smear_sweep and emits its own SmearSweep rows with
        // kind=comp.)
      }
    }
    /************************************************************************/

  }
  /****** Compare fermion TCD definitions against each other (FermFerm) ************/
  // Computes Pearson Corr and normalised IP for every (i,j) pair with i<j.
  // Gated on compute_topo && have_evals (same as topo-charge accumulation).
  //
  // stdout: one greppable line per pair:
  //   "FermFerm: q_X vs q_Y tau_wf=N conf=M Corr=Z IP=W"
  // stdout: compact N×N Corr matrix for visual inspection
  // file:   data_dir/fermferm.dat  (append) — columns: def1 def2 tau_wf conf Corr IP
  if(compute_topo && have_evals){
    typedef typename PeriodicGimplR::ComplexField ComplexField;
    struct FFQDef { std::string name; LatticeComplexD* field; };
    std::vector<FFQDef> fqd;
    for(int t = 0; t < ntracks; t++){
      const std::string& lbl = weight_specs[t].label;
      fqd.push_back({"q_A_"+lbl, &q_mid_all[t]});
      fqd.push_back({"q_B_"+lbl, &q_eps_all[t]});
      fqd.push_back({"q_C_"+lbl, &q_bdy_all[t]});
    }
    fqd.push_back({"q_naive",   &q_naive});
    fqd.push_back({"Sigma_low", &Sigma_low});
    int nd = (int)fqd.size();
    LatticeComplexD one(grid); one = ComplexField::scalar_type(1.0, 0.0);

    // Pre-compute mean-subtracted fields (one pass through each field)
    std::vector<LatticeComplexD> qcen;
    qcen.reserve(nd);
    for(int i = 0; i < nd; i++){
      ComplexD avg = TensorRemove(sum(*fqd[i].field)) / RealD(grid->gSites());
      qcen.push_back(LatticeComplexD(grid));
      qcen.back() = *fqd[i].field - avg * one;
    }

    // Symmetric Corr/IP matrices (upper-triangle computed, rest filled by symmetry)
    std::vector<std::vector<double>> Cmat(nd, std::vector<double>(nd, 1.0));
    std::vector<std::vector<double>> Imat(nd, std::vector<double>(nd, 1.0));

    std::ofstream ff_of;
    if(!data_dir.empty()) ff_of.open(data_dir+"/fermferm.dat", std::ios::app);

    for(int i = 0; i < nd; i++){
      for(int j = i+1; j < nd; j++){
        double corr = real(TensorRemove(sum(qcen[i] * qcen[j])))
                      / std::sqrt(norm2(qcen[i]) * norm2(qcen[j]));
        double ip   = real(TensorRemove(innerProduct(*fqd[i].field, *fqd[j].field)))
                      / std::sqrt(norm2(*fqd[i].field)) / std::sqrt(norm2(*fqd[j].field));
        Cmat[i][j] = Cmat[j][i] = corr;
        Imat[i][j] = Imat[j][i] = ip;
        std::cout << "FermFerm: " << fqd[i].name << " vs " << fqd[j].name
                  << " tau_wf=" << tau_wf << " conf=" << conf_id
                  << " Corr=" << corr << " IP=" << ip << std::endl;
        if(ff_of.is_open())
          ff_of << fqd[i].name << " " << fqd[j].name << " "
                << tau_wf << " " << conf_id << " " << corr << " " << ip << "\n";
      }
    }

    // Compact Corr matrix (column-width adapts to longest name)
    int cw = 9;
    for(auto& q : fqd) cw = std::max(cw, (int)q.name.size() + 1);
    std::cout << GridLogMessage << "FermFerm Corr (tau_wf=" << tau_wf
              << " conf=" << conf_id << "):" << std::endl;
    std::cout << GridLogMessage << std::left << std::setw(cw) << "";
    for(int j = 0; j < nd; j++) std::cout << std::right << std::setw(cw) << fqd[j].name;
    std::cout << std::endl;
    for(int i = 0; i < nd; i++){
      std::cout << GridLogMessage << std::left << std::setw(cw) << fqd[i].name;
      for(int j = 0; j < nd; j++)
        std::cout << std::right << std::fixed << std::setprecision(4) << std::setw(cw) << Cmat[i][j];
      std::cout << std::defaultfloat << std::endl;
    }
  }
  /****** Compare fermion TCD definitions vs qlat reference fields (--comp_file) ***/
  // --comp_file <file>[,<file2>,...]
  //   SCIDAC files from qlat: each is a separate stochastic DWF TCD estimate
  //   (e.g. topo_field_0.scidac, topo_field_1.scidac from pickle_to_scidac.ipynb).
  //
  //   Two independent comparisons are attempted (each gated separately):
  //   (a) Fermion TCD (q_A/B/C) vs qlat reference — requires compute_topo && have_evals.
  //       Output: "CompRef: def comp_idx conf Q_evec Q_ref Corr IP rms_diff"
  //       Written to comp_ref_<def>.dat.
  //   (b) Qlat reference vs gluonic TCD — requires --topo_compare && --files1.
  //       Output: "TopoCompRef: comp TD_tau conf Corr IP" → corr_ip_stoch.dat.
  //       Runs even without eigenvector data (τ_FW=0 diagnostic).
  //
  //   Comp fields are appended to data1 for the visualisation pipeline only when (a) runs.
  bool do_fermion_comp = compute_topo && have_evals;
  bool do_gluon_comp   = (bool)topo_compare && !data1.empty();
  std::vector<LatticeComplexD> comp_fields;   // populated below; appended to data1 at end
  if(GridCmdOptionExists(argv,argv+argc,"--comp_file") && (do_fermion_comp || do_gluon_comp)){
    arg = GridCmdOptionPayload(argv,argv+argc,"--comp_file");
    std::vector<std::string> comp_fnames;
    GridCmdOptionCSL(arg, comp_fnames);

    typedef typename PeriodicGimplR::ComplexField ComplexField;
    LatticeComplexD one(grid); one = ComplexField::scalar_type(1.0, 0.0);

    // Build qdef list — only used when do_fermion_comp is true, but defined here
    // so it is in scope for the table header block below.
    struct QDef { std::string name; LatticeComplexD* field; };
    std::vector<QDef> qdefs;
    if(do_fermion_comp){
      for(int t = 0; t < ntracks; t++){
        const std::string& lbl = weight_specs[t].label;
        qdefs.push_back({"q_A_"+lbl, &q_mid_all[t]});
        qdefs.push_back({"q_B_"+lbl, &q_eps_all[t]});
        qdefs.push_back({"q_C_"+lbl, &q_bdy_all[t]});
      }
      qdefs.push_back({"q_naive",   &q_naive});
      qdefs.push_back({"Sigma_low", &Sigma_low});
      qdefs.push_back({"B_bulk",    &B_bulk});
      qdefs.push_back({"q_B_unit", &q_B_unit});
    }

    // Pre-load all comp files BEFORE printing the header so that LIME/IOobject
    // messages from readFile do not interleave with the table rows.
    std::vector<LatticeComplexD> refs;
    std::vector<double>          Q_refs;
    std::vector<ComplexD>        avg_refs;
    for(int ci = 0; ci < (int)comp_fnames.size(); ci++){
      refs.emplace_back(grid);
      readFile(refs.back(), comp_fnames[ci]);
      Q_refs .push_back(real(TensorRemove(sum(refs.back()))));
      avg_refs.push_back(TensorRemove(sum(refs.back())) / RealD(grid->gSites()));
      std::cout << GridLogMessage << "comp_file[" << ci << "] " << comp_fnames[ci]
                << "  Q=" << Q_refs.back() << std::endl;
    }

    // Header for human-readable fermion-q-vs-qlat table (only when fermion comparison runs)
    if(do_fermion_comp){
      std::cout << GridLogMessage
                << std::left  << std::setw(5)  << "ci"
                << std::setw(10) << "def"
                << std::setw(8)  << "tau_wf"
                << std::right << std::setw(16) << "Q_evec"
                              << std::setw(16) << "Q_ref"
                              << std::setw(16) << "Corr"
                              << std::setw(16) << "IP"
                              << std::setw(16) << "rms_diff" << std::endl;
    }

    for(int ci = 0; ci < (int)comp_fnames.size(); ci++){
      LatticeComplexD& ref   = refs[ci];
      double           Q_ref = Q_refs[ci];
      ComplexD       avg_ref = avg_refs[ci];

      // (a) Compare fermion TCD definitions (q_A/B/C) vs qlat reference
      if(do_fermion_comp){
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
                    << std::setw(10) << qd.name
                    << std::setw(8)  << tau_wf
                    << std::right << std::setw(16) << Q_evec
                                  << std::setw(16) << Q_ref
                                  << std::setw(16) << corr
                                  << std::setw(16) << ip
                                  << std::setw(16) << rms_diff << std::endl;

          // Write comparison result: comp_idx  tau_wf  conf  Q_evec  Q_ref  Corr  IP  rms_diff
          if(!data_dir.empty()){
            std::ofstream of(data_dir+"/comp_ref_"+qd.name+".dat", std::ios::app);
            of << ci << " " << tau_wf << " " << conf_id << " "
               << Q_evec << " " << Q_ref << " " << corr << " " << ip << " " << rms_diff << "\n";
          } else {
            *topo_log << "CompRef: " << qd.name << " " << ci << " " << conf_id << " "
                      << Q_evec << " " << Q_ref << " " << corr << " " << ip << " " << rms_diff << std::endl;
          }
        }
      }

      // (b) Compare qlat stochastic field vs gluonic TCD.
      // Gated on do_gluon_comp (--topo_compare + --files1) rather than have_evals:
      // this comparison is purely field-vs-field and does not require eigenvectors.
      // Runs even at τ_FW=0 when eigenvalue data is unavailable.
      if(do_gluon_comp){
        for(int i = 0; i < (int)data1.size(); i++){
          double   Q_gluon   = real(TensorRemove(sum(data1[i])));
          ComplexD avg_gluon = TensorRemove(sum(data1[i])) / RealD(grid->gSites());
          LatticeComplexD X(grid), Y(grid);
          X = ref      - avg_ref   * one;
          Y = data1[i] - avg_gluon * one;
          double corr_sg = real(TensorRemove(sum(X*Y))) / std::sqrt(norm2(X) * norm2(Y));
          X = ref;  Y = data1[i];
          double ip_sg = real(TensorRemove(innerProduct(X,Y))) / std::sqrt(norm2(X)) / std::sqrt(norm2(Y));
          // Always print to stdout for monitoring (tau_wf identifies eigenvector smearing)
          std::cout << "TopoCompRef: comp=" << ci << " tau_wf=" << tau_wf
                    << " TD_tau=" << get_td_tau(i) << " conf=" << conf_id
                    << " Corr=" << corr_sg << " IP=" << ip_sg << std::endl;
          // Write stoch-vs-gluon result: comp_idx  TD_tau  tau_wf  conf  Q_ref  Q_gluon  Corr  IP
          if(!data_dir.empty()){
            std::ofstream of(data_dir+"/corr_ip_stoch.dat", std::ios::app);
            of << ci << " " << get_td_tau(i) << " " << tau_wf << " " << conf_id << " "
               << Q_ref << " " << Q_gluon << " " << corr_sg << " " << ip_sg << "\n";
          } else {
            *topo_log << "TopoCompRef: " << ci << " " << i << " " << conf_id << " "
                      << Q_ref << " " << Q_gluon << " " << corr_sg << " " << ip_sg << std::endl;
          }
        }
      }

      // Stash for appending to data1 (visualisation frames) below.
      // Only needed when fermion fields were computed (do_fermion_comp); otherwise
      // no q_A/B/C data exists to display alongside the reference.
      if(do_fermion_comp)
        comp_fields.push_back(refs[ci]);
    }

    /****** Sigma sweep — Luchang vs gluonic (sec:smear_sweep, comp arm) ****/
    // Apply the same Gaussian density smear to both q_L and each gluonic
    // q_g(tau_WG), then measure PCF vs sigma.  q_L is UV-sharp (single
    // m=1 Dirac trace) so PCF(q_L, q_g(tau_WG=0)) is large but
    // PCF(q_L, q_g(tau_WG=16)) is small.  Smearing both should equalise:
    // the bandwidth-mismatch story makes a sharp prediction.
    // Output schema appended to smear_sweep.dat: kind name sigma tau_wf TD_tau conf Corr IP
    if(do_gluon_comp && GridCmdOptionExists(argv,argv+argc,"--smear_sweep")
       && !data1.empty() && !refs.empty()){
      std::vector<double> sigmas = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
      std::string s_arg = GridCmdOptionPayload(argv,argv+argc,"--smear_sweep");
      if(!s_arg.empty()){
        std::vector<std::string> s_strs;
        GridCmdOptionCSL(s_arg, s_strs);
        if(!s_strs.empty()){
          sigmas.clear();
          for(auto& s : s_strs) sigmas.push_back(std::stod(s));
        }
      }
      LatticeComplexD oneL(grid); oneL = ComplexD(1.0, 0.0);
      for(double sigma : sigmas){
        std::vector<LatticeComplexD> data1_smr;
        data1_smr.reserve(data1.size());
        for(int i = 0; i < (int)data1.size(); i++)
          data1_smr.push_back(gaussian_smear(data1[i], sigma));
        for(int ci = 0; ci < (int)refs.size(); ci++){
          LatticeComplexD ref_smr = gaussian_smear(refs[ci], sigma);
          for(int i = 0; i < (int)data1.size(); i++){
            ComplexD avg1 = TensorRemove(sum(data1_smr[i])) / RealD(grid->gSites());
            ComplexD avgR = TensorRemove(sum(ref_smr))      / RealD(grid->gSites());
            LatticeComplexD X(grid), Y(grid);
            X = data1_smr[i] - avg1 * oneL;
            Y = ref_smr      - avgR * oneL;
            double n2X = norm2(X), n2Y = norm2(Y);
            double corr = (n2X > 0 && n2Y > 0)
                        ? real(TensorRemove(sum(X*Y))) / std::sqrt(n2X * n2Y)
                        : 0.0;
            double n21 = norm2(data1_smr[i]), n2R = norm2(ref_smr);
            double ip   = (n21 > 0 && n2R > 0)
                        ? real(TensorRemove(innerProduct(data1_smr[i], ref_smr)))
                            / std::sqrt(n21 * n2R)
                        : 0.0;
            std::cout << "SmearSweep: kind=comp name=comp_" << ci
                      << " sigma=" << sigma
                      << " tau_wf=" << tau_wf
                      << " TD_tau=" << get_td_tau(i)
                      << " conf=" << conf_id
                      << " Corr=" << corr << " IP=" << ip << std::endl;
            if(!data_dir.empty()){
              std::ofstream of(data_dir+"/smear_sweep.dat", std::ios::app);
              of << "comp comp_" << ci << " " << sigma << " "
                 << tau_wf << " " << get_td_tau(i) << " " << conf_id << " "
                 << corr << " " << ip << "\n";
            }
          }
        }
      }
    }
    /************************************************************************/
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
