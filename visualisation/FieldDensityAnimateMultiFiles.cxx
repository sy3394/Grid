// Derived from VTK/Examples/Cxx/Medical2.cxx
// The example reads a volume dataset, extracts two isosurfaces that
// represent the skin and bone, and then displays them.
//
// Modified heavily by Peter Boyle to display lattice field theory data as movies and compare multiple files

#define AXES

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

#include <vtkAxesActor.h>
#include <vtkAxesActor.h>
#ifdef AXES
#include <vtkOrientationMarkerWidget.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>
#endif

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

#include <Grid/Grid.h>
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

#define USE_FLYING_EDGES
#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
typedef vtkFlyingEdges3D isosurface;
#else
#include <vtkMarchingCubes.h>
typedef vtkMarchingCubes isosurface;
#endif

int mpeg = 0;
int Ls = -1;
int take_diff = 0;
int dynm_dir = 3;           // set by --animate; default: T
std::vector<int> omit_dirs(1,4);
std::vector<int> omit_intcpts(1,0);
std::vector<int> xlate_omit_dirs;
bool save_file=0;

std::vector<std::string> dynm_labels = {"X", "Y", "Z", "T", "configs"};
std::ofstream data_f;

// HDF5 reader: dataset "field" shape (Nsites,2) float64, Grid lex order (x fastest)
// Uses the HDF5 C API (hdf5.h) — works even when C++ bindings (H5Cpp.h) are absent.
template <class T> void readFileHDF5(T& out, std::string const fname){
#if defined(HAVE_HDF5)
  typedef typename T::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  Grid::GridBase *grid   = out.Grid();
  int64_t         Nsites = grid->_gsites;

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
using namespace Grid;

int n_dims = Nd; //can remove this if define it in class FrameUpdater based on ext_latt_size

class FrameUpdater : public vtkCallbackCommand
{
public:

  FrameUpdater() {
    TimerCount = 0;
    x3         = 0;
    imageData  = nullptr;
    grid_data.clear();
    ext_latt_size = nullptr;
    coor_map   = nullptr;
    timerId    = 0;
    maxCount   = -1;
    dynmIndexF.clear();
    use_dynmIndexF = 0;
  }
  
  static FrameUpdater* New()
  {
    FrameUpdater* cb = new FrameUpdater;
    cb->TimerCount = 0;
    return cb;
  }

  virtual void Execute(vtkObject* caller, unsigned long eventId,void* vtkNotUsed(callData))
  {
    const int max=256;
    char text_string[max];

    if (this->TimerCount < this->maxCount) {

      if (vtkCommand::TimerEvent == eventId)
	{
	  ++this->TimerCount;


	  /*****  Make a new frame  **********/
	  for(int x0=0;x0<ext_latt_size[coor_map[0]];x0++){
	    for(int x1=0;x1<ext_latt_size[coor_map[1]];x1++){
	      for(int x2=0;x2<ext_latt_size[coor_map[2]];x2++){
		
		RealD value = 0.0;
		Coordinate site(std::vector(n_dims,0));
		
		// the first two frame dim's are always latt dim
		site[coor_map[0]] = x0; site[coor_map[1]] = x1;
		if(dynm_dir<n_dims) site[dynm_dir] = x3;    // dynm_dir != file index => coor_map[2] != latt index 
		else                site[coor_map[2]] = x2; // dynm_dir == file index => coor_map[2] == latt index

		for(int i=0; i<(int)omit_dirs.size(); i++) site[omit_dirs[i]] = omit_intcpts[i];
		/***  The last elem of omit_dirs can be the file index  ***/
		
		// The last omit dir != file index => one frame index can be a file index
		//   Recall: omit dir is the dir not one of the frame axes)
		if(omit_dirs.back()<n_dims){ 
		  site[omit_dirs.back()] = omit_intcpts.back();
		  value = coor_map[2] == n_dims? real(peekSite(*grid_data[x2],site)) : real(peekSite(*grid_data[x3],site));
		}
		// The last omit dir == file index => all frame dims are latt dims
		//   Recall: omit_dir != dynm_dir => each file content is shown in separate frames && x3 is used above
		else { 
		  site[coor_map[2]] = x2;
		  value = real(peekSite(*grid_data[0],site));
		}

		imageData->SetScalarComponentFromDouble(x0,x1,x2,0,value);
		if(save_file) data_f<<value<<std::endl;

	  }}}

	  
	  /*****   Put frame counter on the upper left corner   ***********/
	  if(use_dynmIndexF)
	    snprintf(text_string,max,"%s=%s",dynm_labels[dynm_dir].c_str(),dynmIndexF[x3].c_str());
	  else
	    snprintf(text_string,max,"%s=%d",dynm_labels[dynm_dir].c_str(),x3);
	  for(int i_omit=0; i_omit<(int)omit_dirs.size(); i_omit++){
	    if(std::find(xlate_omit_dirs.begin(), xlate_omit_dirs.end(), i_omit) != xlate_omit_dirs.end()){ // corresp. omit_dir == xlated_dir
	      char tmp[max];
	      strncpy(tmp,text_string, max);
	      snprintf(text_string,max,"%s %s=%d", tmp,
		       dynm_labels[omit_dirs[i_omit]].c_str(), omit_intcpts[i_omit]);
	    }
	  }

	  text->SetInput(text_string);
      

	  /*****  Advance animate index; on wrap, step any cycle dirs  *****/
	  x3 = (x3+1)%ext_latt_size[dynm_dir];
	  if( x3==0 ) {
	    for(int i_ind=0 ; i_ind<(int)xlate_omit_dirs.size(); i_ind++){
	      if(i_ind == 0 )
		omit_intcpts[xlate_omit_dirs[i_ind]] = (omit_intcpts[xlate_omit_dirs[i_ind]]+1)%ext_latt_size[omit_dirs[xlate_omit_dirs[i_ind]]];
	      else if(omit_intcpts[xlate_omit_dirs[i_ind-1]] == 0)
		omit_intcpts[xlate_omit_dirs[i_ind]] = (omit_intcpts[xlate_omit_dirs[i_ind]]+1)%ext_latt_size[omit_dirs[xlate_omit_dirs[i_ind]]];
	    }
	  }

	  /*****   Print the log to stdout   ***********/
	  std::cout << this->TimerCount<<"/"<<maxCount
		    << " " << dynm_labels[dynm_dir] <<"="<< x3 <<" ";
	  for(int ind : xlate_omit_dirs) std::cout << dynm_labels[omit_dirs[ind]] <<"="<<omit_intcpts[ind] << " ";
	  if(use_dynmIndexF) std::cout<<dynmIndexF[x3];
	  std::cout<<std::endl;
	  imageData->Modified();

	  vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
	  iren->GetRenderWindow()->Render();
	  
	}
    }
    
    if (this->TimerCount >= this->maxCount) {
      vtkRenderWindowInteractor* iren = dynamic_cast<vtkRenderWindowInteractor*>(caller);
      if (this->timerId > -1)
      {
        iren->DestroyTimer(this->timerId);
      }
    }
  }

private:
  int TimerCount;
  int x3;
public:
  std::vector<Grid::LatticeComplexD *> grid_data;
  int* ext_latt_size;
  int* coor_map;
  vtkImageData* imageData = nullptr;
  vtkTextActor* text = nullptr;
  int timerId ;
  int maxCount ;
  double rms;
  isosurface * posExtractor;
  isosurface * negExtractor;
  std::vector<std::string> dynmIndexF;
  bool use_dynmIndexF;
};


class SliderCallback : public vtkCommand
{
public:
    static SliderCallback* New()
    {
        return new SliderCallback;
    }
    virtual void Execute(vtkObject* caller, unsigned long eventId, void* callData)
    {
        vtkSliderWidget *sliderWidget = vtkSliderWidget::SafeDownCast(caller);
        if (sliderWidget)
        {
	  contour = ((vtkSliderRepresentation *)sliderWidget->GetRepresentation())->GetValue();
        }
	for(int i=0;i<fu_list.size();i++){
	  fu_list[i]->posExtractor->SetValue(0,  SliderCallback::contour*fu_list[i]->rms);
	  fu_list[i]->negExtractor->SetValue(0, -SliderCallback::contour*fu_list[i]->rms);
	  fu_list[i]->posExtractor->Modified();
	  fu_list[i]->negExtractor->Modified();
	}
    }
public:
  static double contour;
  std::vector<FrameUpdater *> fu_list;
};


double SliderCallback::contour;

int main(int argc, char* argv[])
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian   UGrid(latt_size, simd_layout, mpi_layout);
  GridCartesian * grid  = &UGrid;
    
 
  /***********   READ INPUT           *******************************************
  CLI layout (direction names: X Y Z T [Ls] configs):
    --animate dir       fast animation axis (default: T)
    --fix     dir=N     hold 'dir' at slice N (repeatable)
    --cycle   dir[=N]   outer loop: step 'dir' after each animate wrap (repeatable)
    --sum     dir       sum over all slices of 'dir' before display (repeatable)
  Display axes are the 3 dirs left unassigned.
  ************************************************************************************/
  std::string separator = "smr.";
  std::vector<std::string> file_list;
  double default_contour = 1.0;
  bool use_fname_as_frame_counter = false;
  std::ifstream index_list;
  
  std::string arg;
  
  std::cout << argc << " command Line arguments "<<std::endl;
  for(int c=0;c<argc;c++) {
    std::cout << " - "<<argv[c]<<std::endl;
  }

#ifdef MPEG
  std::string mpeg_fname = "movie.avi";
  if( GridCmdOptionExists(argv,argv+argc,"--mpeg") ){
    mpeg = 1;
    arg = GridCmdOptionPayload(argv,argv+argc,"--mpeg");
    if(!arg.empty()) mpeg_fname = arg;
  }
#endif
  if( GridCmdOptionExists(argv,argv+argc,"--Ls") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--Ls");
    GridCmdOptionInt(arg,Ls);
    assert(Ls !=0 && "Ls needs to be greater than 0 when specified");
    grid = SpaceTimeGrid::makeFiveDimGrid(Ls, &UGrid);
    latt_size = grid->GlobalDimensions();
    dynm_labels.insert(dynm_labels.begin(),"Ls");
    n_dims++;
    std::cout<<grid->GlobalDimensions()<<" "<<latt_size<<std::endl;
  }

  // =========================================================================
  // Animation layout: --animate / --fix / --cycle / --sum
  // =========================================================================
  // Build direction name map after --Ls (which may shift indices).
  // latt_size.size() == n_dims == 4 (no Ls) or 5 (with Ls).
  // configs is always at index n_dims (one past the last lattice dim).
  std::map<std::string,int> dir_map;
  if(Ls > 0)
    dir_map = {{"Ls",0},{"X",1},{"Y",2},{"Z",3},{"T",4},{"configs",n_dims}};
  else
    dir_map = {{"X",0},{"Y",1},{"Z",2},{"T",3},{"configs",n_dims}};
  const int n_dirs_total = n_dims + 1; // lattice dims + configs

  auto dir_name = [&](int d) -> std::string {
    for(auto& kv : dir_map) if(kv.second == d) return kv.first;
    return std::to_string(d);
  };
  auto parse_dir_val = [&](const char* tok, int& d, int& v) -> bool {
    std::string s(tok);
    auto eq = s.find('=');
    std::string name = (eq != std::string::npos) ? s.substr(0,eq) : s;
    v = (eq != std::string::npos) ? std::stoi(s.substr(eq+1)) : 0;
    auto it = dir_map.find(name);
    if(it == dir_map.end()) return false;
    d = it->second; return true;
  };

  // Default animate direction: T (index 3 in 4D, index 4 in 5D)
  std::string animate_name = (Ls > 0) ? "T" : "T";
  if( GridCmdOptionExists(argv,argv+argc,"--animate") )
    animate_name = GridCmdOptionPayload(argv,argv+argc,"--animate");

  std::map<int,int> fix_init;   // dir → initial/fixed slice
  std::set<int>     cycle_dirs; // dirs that cycle as outer loop
  std::set<int>     sum_dirs;   // dirs to sum over before display
  for(int i = 1; i < argc-1; i++){
    int d, v;
    if(std::string(argv[i]) == "--fix"){
      if(!parse_dir_val(argv[i+1], d, v)){
        std::cerr << "ERROR: unknown direction in --fix " << argv[i+1] << std::endl; exit(1);
      }
      fix_init[d] = v; i++;
    } else if(std::string(argv[i]) == "--cycle"){
      if(!parse_dir_val(argv[i+1], d, v)){
        std::cerr << "ERROR: unknown direction in --cycle " << argv[i+1] << std::endl; exit(1);
      }
      fix_init[d] = v; cycle_dirs.insert(d); i++;
    } else if(std::string(argv[i]) == "--sum"){
      std::string s(argv[i+1]);
      auto it = dir_map.find(s);
      if(it == dir_map.end()){
        std::cerr << "ERROR: unknown direction in --sum " << argv[i+1] << std::endl; exit(1);
      }
      sum_dirs.insert(it->second); i++;
    }
  }

  // Translate --animate → dynm_dir
  {
    auto it = dir_map.find(animate_name);
    if(it == dir_map.end()){
      std::cerr << "ERROR: unknown --animate direction '" << animate_name << "'" << std::endl; exit(1);
    }
    dynm_dir = it->second;
  }

  // --- Consistency check ---
  {
    std::map<int,std::string> role;
    role[dynm_dir] = "animate";
    std::vector<std::string> errs;
    auto assign = [&](int d, const std::string& r){
      if(role.count(d))
        errs.push_back("'" + dir_name(d) + "' assigned to both '" + role[d] + "' and '" + r + "'");
      else role[d] = r;
    };
    for(auto& [d,v] : fix_init)
      assign(d, cycle_dirs.count(d) ? "cycle" : "fix");
    for(int d : sum_dirs) assign(d, "sum");

    if(!errs.empty()){
      for(auto& e : errs) std::cerr << "ERROR: " << e << std::endl; exit(1);
    }
    std::vector<int> disp;
    for(int d = 0; d < n_dirs_total; d++) if(!role.count(d)) disp.push_back(d);
    if((int)disp.size() != 3){
      std::cerr << "ERROR: " << disp.size() << " display dims inferred (need exactly 3).\n"
                << "  Unassigned:";
      for(int d : disp) std::cerr << " " << dir_name(d);
      std::cerr << "\n  Assign extras via --fix, --cycle, --sum, or --animate." << std::endl;
      exit(1);
    }
    std::cout << "Display axes:";
    for(int d : disp) std::cout << " " << dir_name(d);
    std::cout << "  animate=" << dir_name(dynm_dir) << std::endl;
  }

  // --- Translate to internal representation ---
  omit_dirs.clear(); omit_intcpts.clear(); xlate_omit_dirs.clear();
  for(auto& [d,v] : fix_init)
    if(!cycle_dirs.count(d)){ omit_dirs.push_back(d); omit_intcpts.push_back(v); }
  for(int d : cycle_dirs){
    xlate_omit_dirs.push_back((int)omit_dirs.size());
    omit_dirs.push_back(d);
    omit_intcpts.push_back(fix_init.count(d) ? fix_init.at(d) : 0);
  }
  for(int d : sum_dirs){ omit_dirs.push_back(d); omit_intcpts.push_back(-1); }
  std::cout << "omit_dirs:";
  for(int d : omit_dirs) std::cout << " " << dir_name(d);
  std::cout << "  cycle:";
  for(int i : xlate_omit_dirs) std::cout << " " << dir_name(omit_dirs[i]);
  std::cout << "  sum:";
  for(int d : sum_dirs) std::cout << " " << dir_name(d);
  std::cout << std::endl;

  if( GridCmdOptionExists(argv,argv+argc,"--take_adj_diff") ){
    take_diff = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--isosurface") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--isosurface");
    GridCmdOptionFloat(arg,default_contour);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--files") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--files");
    std::cout<<"files "<<arg<<std::endl;
    GridCmdOptionCSL(arg,file_list);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--dynm_label") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--dynm_label");
    if(!arg.empty()){
      dynm_labels.pop_back();
      dynm_labels.push_back(arg);
    }
  }
  if( GridCmdOptionExists(argv,argv+argc,"--use_fname_as_frame_counter") ){
    use_fname_as_frame_counter = true;
    arg = GridCmdOptionPayload(argv,argv+argc,"--use_fname_as_frame_counter");
    if(!arg.empty()) separator = arg;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--index_file") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--index_file");
    std::cout<<"index file "<<arg<<std::endl;
    if(!arg.empty()) index_list.open(arg,std::ifstream::in);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--save_data_to") ){
    save_file = 1;
    arg = GridCmdOptionPayload(argv,argv+argc,"--save_data_to");
    if(!arg.empty()) data_f.open(arg,std::ios::trunc);
  }

  /***************************************************************/
  /********   Preprocess Data   **********************************/
  /***************************************************************/
  // Read in data
  std::vector<LatticeComplexD> data(file_list.size(),grid);
  for(int c=0;c<(int)data.size();c++) {
    std::cout << "Reading file: "<<file_list[c]<<std::endl;
    readFile(data[c],file_list[c]);
  }

  // --comp_file <file>[,<file2>,...]: append extra fields (HDF5 or SCIDAC) for side-by-side visualisation
  if( GridCmdOptionExists(argv,argv+argc,"--comp_file") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--comp_file");
    std::vector<std::string> comp_fnames;
    GridCmdOptionCSL(arg, comp_fnames);
    for(auto const &cf : comp_fnames){
      std::cout << "Reading comp_file: " << cf << std::endl;
      LatticeComplexD ref(grid);
      readFile(ref, cf);
      data.push_back(ref);
    }
  }

  int flag = 0; for(auto dir:omit_intcpts) if(dir<0) flag++;
  std::string display_info = flag ? "summed over ":"";
  for( int i_od=0; i_od<(int)omit_dirs.size(); i_od++){
    int odir    = omit_dirs[i_od];
    int intcpts = omit_intcpts[i_od];

    // Sum over odir if intcpts<0
    if(intcpts<0){
      typedef typename PeriodicGimplR::ComplexField ComplexField;
      ComplexField filter(grid), ones(grid), zeros(grid); ones = ComplexField::scalar_type(1.0,0.0); zeros = Zero();
      LatticeComplexD Fsum(grid);
      std::cout<<"before sum "<<odir<<" "<<intcpts<<std::endl;
      // odir != file_index
      if(odir < latt_size.size()){
	Lattice<iScalar<vInteger> > x_odir(grid); LatticeCoordinate(x_odir,odir);
	for(int c=0;c<(int)data.size();c++) {
	  Fsum = Zero();
	  for(int i=0; i<latt_size[odir];i++){
	    filter = where( x_odir==i, ones, zeros);
	    //ExtractSlice(tmp_F,data[c],i,0);
	    Fsum = Fsum + Cshift(filter*data[c],odir,i);
	  }
	  data[c] = Fsum;
	}
	display_info += dynm_labels[odir]+",";
	std::cout<<"Sum Done " <<latt_size[odir]<<std::endl;

      }
      // Sum over odir where odir == file_index
      else if(odir == latt_size.size()){
	Fsum = Zero();
	for(int c=0;c<(int)data.size();c++)
	  Fsum = Fsum + data[c];
	data.clear();
	data.push_back(Fsum);
	display_info += "all files,";
      }
      omit_intcpts[i_od] = 0;
    }
  }
  if(flag) display_info.pop_back();
  
  // take diff when demanded
  if(take_diff){
    assert(data.size() > 1 && "take_diff requires at least 2 files");
    for(int c=0;c<(int)data.size()-1;c++)
      data[c] = data[c+1] - data[c];
    data.pop_back();
  }
  
  /****************************************************************/
  /**************      Determine Frame Dimensions    **************/
  /****************************************************************/
  int coor_map[3] = {0,1,2};
  std::vector<int> ext_latt_size(latt_size.size()+1);
  for(int d=0,d_c=0; d<3; d++, d_c++){
    if( d_c == dynm_dir ) d_c++;
    for(auto omit_dir:omit_dirs) if( omit_dir == d_c ) d_c++;
    coor_map[d] = d_c; 
  }
  for(int d=0; d<(int)latt_size.size();d++) ext_latt_size[d] = latt_size[d];
  ext_latt_size[latt_size.size()] = data.size();

  // DEBUG
  for(auto F: data) std::cout<<"Max: "<<std::sqrt(maxLocalNorm2(F))<<" "<<real(TensorRemove(sum(F)))<<" "<<std::sqrt(maxLocalNorm2(F)/norm2(F)*grid->gSites() )<<std::endl;
  for(int d=0;d<3;d++) std::cout<<d<<" "<<coor_map[d]<<std::endl;
  for(int s : ext_latt_size) std::cout<<s<<std::endl;
  for(int s : omit_dirs) std::cout<<s<<std::endl;
  for(int s : omit_intcpts) std::cout<<s<<std::endl;

  /****************************************************************/
  /****************    Setup Frames    ****************************/
  /****************************************************************/
  // Common things:
  vtkNew<vtkNamedColors> colors;
  std::array<unsigned char, 4> posColor{{240, 184, 160, 255}};  colors->SetColor("posColor", posColor.data());
  std::array<unsigned char, 4> bkg{{51, 77, 102, 255}};         colors->SetColor("BkgColor", bkg.data());

  // Create the renderer, the render window, and the interactor. The renderer
  // draws into the render window, the interactor enables mouse- and
  // keyboard-based interaction with the data within the render window.
  //
  vtkNew<vtkRenderWindow> renWin;
  renWin->SetOffScreenRendering(1);   // headless: no X display needed (Frontier, batch nodes)
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  // Total frame count: animate range × outer-loop (cycle) ranges
  int frameCount = ext_latt_size[dynm_dir];
  for(auto ind : xlate_omit_dirs) frameCount *= ext_latt_size[omit_dirs[ind]];

  // fc = number of panels shown simultaneously.
  // If configs is the animate axis → one panel cycling through all files.
  // If configs is fixed/omitted    → each file gets its own panel.
  int configs_idx = dir_map.at("configs");
  int fc = (dynm_dir == configs_idx) ? 1 : (int)data.size();
  std::vector<FrameUpdater *> fu_list;
  double max_cntr = 0; // max slider value; used in interactive slider widget
  for (int f=0;f<fc;f++){

    // It is convenient to create an initial view of the data. The FocalPoint
    // and Position form a vector direction. Later on (ResetCamera() method)
    // this vector is used to position the camera to look at the data in
    // this direction.
    vtkNew<vtkCamera> aCamera;
    aCamera->SetViewUp(0, 0, -1);
    aCamera->SetPosition(0, -1000, 0);
    aCamera->SetFocalPoint(0, 0, 0);
    aCamera->ComputeViewPlaneNormal();
    aCamera->Azimuth(30.0);
    aCamera->Elevation(30.0);

    
    vtkNew<vtkRenderer> aRenderer;
    renWin->AddRenderer(aRenderer);

    //////// Set contour scale
    double vol = data[f].Grid()->gSites();

    auto nrm              = norm2(data[f]);
    auto rms              = sqrt(nrm/vol);
    double max_local_norm = std::sqrt(maxLocalNorm2(data[f]));
    if(max_cntr<max_local_norm/rms) max_cntr = max_local_norm/rms;

    double contour = default_contour < 0 ? -default_contour * max_local_norm : default_contour * rms; // default to 1 x RMS

    std::cout<<"ratio(|data|_inf/|data|_l2: "<<max_local_norm/rms<<" "<<max_local_norm<<" "<<rms<<std::endl;
    std::cout<<"defalt contour: "<<default_contour<<" countour "<<contour<<" rms: "<<rms<<std::endl;

    
    // The following reader is used to read a series of 2D slices (images)
    // that compose the volume. The slice dimensions are set, and the
    // pixel spacing. The data Endianness must also be specified. The reader
    // uses the FilePrefix in combination with the slice number to construct
    // filenames using the format FilePrefix.%d. (In this case the FilePrefix
    // is the root name of the file: quarter.)
    vtkNew<vtkImageData> imageData;
    imageData->SetDimensions(ext_latt_size[coor_map[0]],ext_latt_size[coor_map[1]],ext_latt_size[coor_map[2]]);
    imageData->AllocateScalars(VTK_DOUBLE, 1);
    
    vtkNew<isosurface> posExtractor;
    posExtractor->SetInputData(imageData);
    posExtractor->SetValue(0, contour);
  
    vtkNew<vtkStripper> posStripper;
    posStripper->SetInputConnection(posExtractor->GetOutputPort());

    vtkNew<vtkPolyDataMapper> posMapper;
    posMapper->SetInputConnection(posStripper->GetOutputPort());
    posMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> pos;
    pos->SetMapper(posMapper);
    pos->GetProperty()->SetDiffuseColor(colors->GetColor3d("posColor").GetData());
    pos->GetProperty()->SetSpecular(0.3);
    pos->GetProperty()->SetSpecularPower(20);
    pos->GetProperty()->SetOpacity(0.5);

    // An isosurface, or contour value is set
    // The triangle stripper is used to create triangle strips from the
    // isosurface; these render much faster on may systems.
    vtkNew<isosurface> negExtractor;
    negExtractor->SetInputData(imageData);
    negExtractor->SetValue(0, -contour);

    vtkNew<vtkStripper> negStripper;
    negStripper->SetInputConnection(negExtractor->GetOutputPort());

    vtkNew<vtkPolyDataMapper> negMapper;
    negMapper->SetInputConnection(negStripper->GetOutputPort());
    negMapper->ScalarVisibilityOff();

    vtkNew<vtkActor> neg;
    neg->SetMapper(negMapper);
    neg->GetProperty()->SetDiffuseColor(colors->GetColor3d("Ivory").GetData());

    // An outline provides context around the data.
    vtkNew<vtkOutlineFilter> outlineData;
    outlineData->SetInputData(imageData);

    vtkNew<vtkPolyDataMapper> mapOutline;
    mapOutline->SetInputConnection(outlineData->GetOutputPort());

    vtkNew<vtkActor> outline;
    outline->SetMapper(mapOutline);
    outline->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    ////////// create a label of the frame
    std::string txt = omit_dirs.back()<latt_size.size()?"All Files: "+display_info : file_list[f];
    if(take_diff) txt = "Diff: Next File - "+txt;
    if(mpeg) txt += " cntr="+std::to_string(contour);
    vtkNew<vtkTextActor> Text;
    Text->SetInput(txt.c_str());
    Text->SetPosition2(0,0);
    Text->GetTextProperty()->SetFontSize(48);
    Text->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());

    vtkNew<vtkTextActor> TextT;
    TextT->SetInput((dynm_labels[dynm_dir]+"=0").c_str());
    TextT->SetPosition(0,.9*1025);
    TextT->GetTextProperty()->SetFontSize(48);
    TextT->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());
    
#ifdef AXES
    // https://examples.vtk.org/site/Cxx/GeometricObjects/Axes/
    // Add axis labels
    vtkNew<vtkAxesActor> axes;

    axes->SetXAxisLabelText(dynm_labels[coor_map[0]].c_str());
    axes->SetYAxisLabelText(dynm_labels[coor_map[1]].c_str());
    axes->SetZAxisLabelText(dynm_labels[coor_map[2]].c_str());

    vtkNew<vtkTransform> transform;
    transform->Translate(-5.0, 0.0, 0.0);
    transform->Scale(2, 2, 2);
    // The axes are positioned with a user transform
    axes->SetUserTransform(transform);
    
    // Example of customizing label properties (e.g., color)
    //axes->GetXAxisLabelProperty()->SetColor(1.0, 0.0, 0.0); // Red X-axis label
    /*
    vtkNew<vtkOrientationMarkerWidget> orientationWidget;
    orientationWidget->SetOrientationMarker(axes);
    orientationWidget->SetInteractor(iren);
    orientationWidget->SetViewport(0.0, 0.0, 0.2, 0.2); // Example: Bottom-left corner, 20% of viewport size
    orientationWidget->EnabledOn();
    */
    aRenderer->AddActor(axes);
#endif
    
    // Actors are added to the renderer. An initial camera view is created.
    // The Dolly() method moves the camera towards the FocalPoint,
    // thereby enlarging the image.
    aRenderer->AddActor(Text);
    aRenderer->AddActor(TextT);
    aRenderer->AddActor(outline);
    aRenderer->AddActor(pos);
    aRenderer->AddActor(neg);

    // Sign up to receive TimerEvent
    std::vector<LatticeComplexD*> tmp;
    if(omit_dirs.back()<latt_size.size()) for(int i=0;i<data.size();i++) tmp.push_back(&data[i]);
    else tmp.push_back(&data[f]);
    
    vtkNew<FrameUpdater> fu;
    fu->imageData = imageData;
    fu->grid_data = tmp;
    fu->ext_latt_size = ext_latt_size.data();
    fu->coor_map  = coor_map;
    fu->text      = TextT;
    fu->maxCount = frameCount;
    fu->posExtractor = posExtractor;
    fu->negExtractor = negExtractor;
    fu->rms = rms;
    // Set frame index
    std::vector<std::string> ind_list;
    if( use_fname_as_frame_counter ) {

      for(const auto& fname : file_list) {
	int str_length = separator.length();
	ind_list.push_back(fname.substr(fname.rfind((separator))+str_length));
	std::cout<<ind_list.back()<<std::endl;
      }
      fu->dynmIndexF = ind_list;
      fu->use_dynmIndexF = true;
    }
    else if(index_list.is_open()){
      assert(!use_fname_as_frame_counter && "file name is not to be used as a frame counter");

      std::string line;
      while(getline(index_list, line)){
	ind_list.push_back(line);
      }
      fu->dynmIndexF = ind_list;
      fu->use_dynmIndexF = true;
    }
    
    iren->AddObserver(vtkCommand::TimerEvent, fu);

    aRenderer->SetActiveCamera(aCamera);
    aRenderer->ResetCamera();
    aRenderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
    aCamera->Dolly(1.0);

    double nf = fc;//file_list.size();
    std::cout << " Adding renderer " <<f<<" of "<<nf<<std::endl;
    aRenderer->SetViewport((1.0/nf)*f, 0.0,(1.0/nf)*(f+1) , 1.0);

    // Note that when camera movement occurs (as it does in the Dolly()
    // method), the clipping planes often need adjusting. Clipping planes
    // consist of two planes: near and far along the view direction. The
    // near plane clips out objects in front of the plane; the far plane
    // clips out objects behind the plane. This way only what is drawn
    // between the planes is actually rendered.
    aRenderer->ResetCameraClippingRange();
    
    fu_list.push_back(fu);
  }

  // Set a background color for the renderer and set the size of the
  // render window (expressed in pixels).
  // Initialize the event loop and then start it.
  renWin->SetSize(1024*fc, 1024);
  renWin->SetWindowName("FieldDensity");
  renWin->Render();

  iren->Initialize();

  double nf = fc;
  if ( mpeg ) {
#ifdef MPEG
    vtkWindowToImageFilter *imageFilter = vtkWindowToImageFilter::New();
    imageFilter->SetInput( renWin );
    imageFilter->SetInputBufferTypeToRGB();
    
    vtkFFMPEGWriter *writer = vtkFFMPEGWriter::New();
    writer->SetFileName(mpeg_fname.c_str());
    writer->SetRate(1);
    writer->SetInputConnection(imageFilter->GetOutputPort());
    writer->Start();

    for(int i=0;i<fu_list[0]->maxCount;i++){
      for(int f=0;f<fu_list.size();f++){
	fu_list[f]->Execute(iren,vtkCommand::TimerEvent,nullptr);
      }
      imageFilter->Modified();
      writer->Write();
    }
    writer->End();
    writer->Delete();
#else
    assert(-1 && "MPEG support not compiled");
#endif
  } else { 

    std::cout << "Max Slide Val: " << max_cntr << std::endl;
    // Add control of contour threshold
    // Create a slider widget
    vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
    sliderRep->SetMinimumValue(0.1);
    sliderRep->SetMaximumValue(max_cntr);
    sliderRep->SetValue(1.0);
    sliderRep->SetTitleText("Fraction RMS");
    // Set color properties:

    // Change the color of the knob that slides
    //  sliderRep->GetSliderProperty()->SetColor(colors->GetColor3d("Green").GetData());
    sliderRep->GetTitleProperty()->SetColor(colors->GetColor3d("AliceBlue").GetData());
    sliderRep->GetLabelProperty()->SetColor(colors->GetColor3d("AliceBlue").GetData());
    sliderRep->GetSelectedProperty()->SetColor(colors->GetColor3d("DeepPink").GetData());

    // Change the color of the bar
    sliderRep->GetTubeProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
    sliderRep->GetCapProperty()->SetColor(colors->GetColor3d("Yellow").GetData());
    sliderRep->SetSliderLength(0.05);
    sliderRep->SetSliderWidth(0.025);
    sliderRep->SetEndCapLength(0.02);

    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(0.1, 0.1);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(0.9/nf, 0.1);
  
    vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();
  
    // Create the slider callback
    vtkSmartPointer<SliderCallback> slidercallback = vtkSmartPointer<SliderCallback>::New();
    slidercallback->fu_list = fu_list;
    sliderWidget->AddObserver(vtkCommand::InteractionEvent, slidercallback);

    int timerId = iren->CreateRepeatingTimer(300);
    std::cout << "timerId: " << timerId << std::endl;

    // Start the interaction and timer
    iren->Start();
  }
  data_f.close();
  
  Grid_finalize();

  return EXIT_SUCCESS;
}
