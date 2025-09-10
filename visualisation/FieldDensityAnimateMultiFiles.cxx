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
int xlate = 0 ;
int take_diff = 0;
int dynm_dir = 3;
std::vector<int> omit_dirs(1,4);
std::vector<int> omit_intcpts(1,0);
int xlate_omit_dir=-1;
bool sum_omit_dir=0;
bool save_file=0;

std::vector<std::string> dynm_labels = {"X", "Y", "Z", "T", "tau"};
std::ofstream data_f;

template <class T> void readFile(T& out, std::string const fname){
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
    xoff       = 0;
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
		site[coor_map[0]] = (x0+xoff)%ext_latt_size[coor_map[0]]; site[coor_map[1]] = x1;
		if(dynm_dir<n_dims) site[dynm_dir] = x3;    // dynm_dir != file index => coor_map[2] != latt index 
		else                site[coor_map[2]] = x2; // dynm_dir == file index => coor_map[2] == latt index

		if(sum_omit_dir){
		  // ASSUME: omit_dir != file_index_dir if summed over the dir
		  int x_max = 1; for(int &dir : omit_dirs) x_max *= ext_latt_size[dir];
		  for(int xi=0; xi<x_max; xi++){
		    int x_tmp = xi;
		    for(int dir: omit_dirs){
		      site[dir] = x_tmp%ext_latt_size[dir];
		      x_tmp /= ext_latt_size[dir];
		    }
		    value += (coor_map[2] == n_dims) ? real(peekSite(*grid_data[x2],site)) : real(peekSite(*grid_data[x3],site));
		  }
		}
		else {
		  for(int i=0; omit_dirs.size()-1; i++) site[omit_dirs[i]] = omit_intcpts[i];
		  /***  The last elem of omit_dirs can be the file index  ***/

		  // The last omit dir !=t file index => one frame index can be a file index
		  if(omit_dirs.back()<n_dims){ 
		    site[omit_dirs.back()] = omit_intcpts.back();
		    value = coor_map[2] == n_dims? real(peekSite(*grid_data[x2],site)) : real(peekSite(*grid_data[x3],site));
		  }
		  // The last omit dir == file index => all frame dims are latt dims
		  else { 
		    site[coor_map[2]] = x2;
		    value = real(peekSite(*grid_data[0],site));
		  }
		}
		imageData->SetScalarComponentFromDouble(x0,x1,x2,0,value);
		if(save_file) data_f<<value<<std::endl;

	  }}}

	  
	  /*****   Put frame counter on the upper left corner   ***********/
	  if(use_dynmIndexF)
	    snprintf(text_string,max,"%s=%s",dynm_labels[dynm_dir].c_str(),dynmIndexF[x3].c_str());
	  else
	    snprintf(text_string,max,"%s=%d",dynm_labels[dynm_dir].c_str(),x3);
	  if(xlate_omit_dir>=0){
	    char tmp[max];
	    strncpy(tmp,text_string, max);
	    snprintf(text_string,max,"%s %s=%d", tmp,
		     dynm_labels[omit_dirs[xlate_omit_dir]].c_str(), omit_intcpts[xlate_omit_dir]);
	  }

	  text->SetInput(text_string);
      

	  /*****  Update frame dims   *************************/
	  if ( xlate ) {
	    // When translating the first frame index, dynamic index is updated only after a complete translation
	    xoff = (xoff + 1)%ext_latt_size[coor_map[0]];
	    if ( xoff== 0 ) x3 = (x3+1)%ext_latt_size[dynm_dir];
	  } else {
	    x3 = (x3+1)%ext_latt_size[dynm_dir];
	    if( xlate_omit_dir>=0 && x3==0 ) omit_intcpts[xlate_omit_dir] = (omit_intcpts[xlate_omit_dir]+1)%ext_latt_size[omit_dirs[xlate_omit_dir]];
	    //if ( x3 == 0 ) 	xoff = (xoff + 1)%ext_latt_size[coor_map[0]];
	  }

	  
	  /*****   Print the log to stdout   ***********/
	  std::cout << this->TimerCount<<"/"<<maxCount<< " xoff "<<xoff<<" t_updated "<< x3 <<" "<< omit_intcpts[xlate_omit_dir];
	  if(use_dynmIndexF) std::cout<<" "<<dynmIndexF[x3]<<std::endl;
	  else std::cout<<std::endl;
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
  int xoff;
  int x3;
public:
  std::vector<Grid::LatticeComplexD *> grid_data;
  int* ext_latt_size;
  int* coor_map;
  vtkImageData* imageData = nullptr;
  vtkTextActor* text = nullptr;
  vtkFFMPEGWriter *writer = nullptr;
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


    
  For now,
    - xlate: the first non-omitted coor
    - n_dim: #dimensions of the lattice
    - dynm_dir: dir along which 3D plot is updated
    - omit_dir: != dynm_dir & dir to be omitted completely from the 3D plot
        ASSUME: it is sorted when mulltiple vals are given
    - default: omit_dir = n_dim = file_index && dynm_dir = n_dim-1 = the last lattice index
                i.e., data in diff files are put to diff frames & time change is reflected in the frame updates
    - if omit_dir < n_dim-1;
        - index on input files is treated as a dimension of the frame
	- if dynm_dir = n_dim; n_dim-1 out of n_dim latt dim's are displayed
	- if dynm_dir < n_dim; only 2 out of 4 latt dim is displayed + sim. time
  ************************************************************************************/
  std::string separator = "smr.";
  std::vector<std::string> file_list, data_fname;
  double default_contour = 1.0;
  int use_fname_as_frame_counter = 0, pre_sum_Ls = 0;
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
    omit_dirs[0] = 5;
    n_dims++;
    std::cout<<grid->GlobalDimensions()<<" "<<latt_size<<std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--pre_sum_Ls") ){
    pre_sum_Ls = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--xlate") ){
    xlate = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--dynm_dir") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--dynm_dir");
    GridCmdOptionInt(arg,dynm_dir);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--omit_dirs") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--omit_dirs");
    GridCmdOptionIntVector(arg,omit_dirs);
    assert( 3 <= latt_size.size() + 1 - omit_dirs.size() && "Too much dirs are omitted to make a 3D plot");
    for(int &omit_dir : omit_dirs){
      assert( omit_dir != dynm_dir && "The omitted dir cannot be the same as updated dimension" );
      assert( !(omit_dir == 4 && sum_omit_dir) && "We do not sum over file index for now" );
    }
  }
  if( GridCmdOptionExists(argv,argv+argc,"--omit_intcpts") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--omit_intcpts");
    GridCmdOptionIntVector(arg,omit_intcpts);
    assert( omit_intcpts.size() == omit_dirs.size() && "omit_dirs and omit_intcpts should have the same size" );
    for(int &dir : omit_intcpts) if( dir==latt_size.size() ) std::cout<<"--omit_intcpt in "<<dir<<" dir has no effect"<<std::endl;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--sum_omit_dir") ){
    sum_omit_dir = 1;
  }
  if( GridCmdOptionExists(argv,argv+argc,"--xlate_omit_dir") ){
    arg=GridCmdOptionPayload(argv,argv+argc,"--xlate_omit_dir");
    GridCmdOptionInt(arg,xlate_omit_dir);
    assert( !sum_omit_dir && xlate_omit_dir<omit_dirs.size() &&
    //assert(!sum_omit_dir && std::find(omit_dirs.begin(), omit_dirs.end(), xlate_omit_dir) != omit_dirs.end() &&
	   "xlate_omit_dir is an index of the array of omit_dirs");
  }

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
  if( GridCmdOptionExists(argv,argv+argc,"--use_fname_as_frame_counter") ){
    use_fname_as_frame_counter = 1;
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
  // Read in data; take diff when demanded
  FieldMetaData header;
  std::vector<LatticeComplexD> data(file_list.size()-take_diff,grid);
  for(int c=0;c<data.size();c++) {
    std::cout << "Reading file: "<<file_list[c]<<std::endl;
    readFile(data[c],file_list[c]); data[c] = data[c];
  }
  if(pre_sum_Ls && Ls>0){
    std::vector<LatticeComplexD> tmp(data.size(), &UGrid);
    LatticeComplexD Fsum(&UGrid), tmp_F(&UGrid); 
    for(int c=0;c<data.size();c++) {
      Fsum = Zero();
      for(int i=0; i<Ls;i++){
	ExtractSlice(tmp_F,data[c],i,0);
        Fsum = Fsum + tmp_F;
      }
      tmp[c] = Fsum;
    }
    data.clear();
    for(int c=0;c<tmp.size();c++)
      data.push_back(tmp[c]);
    grid = &UGrid;
    latt_size = grid->GlobalDimensions();
    dynm_labels.erase(dynm_labels.begin());
    n_dims--;
    std::cout<<"Ls pre summed: New Dimensions = " << grid->GlobalDimensions()<<" "<<latt_size<<" "<<n_dims<<std::endl;
  }
  
  if(take_diff){
    for(int c=0;c<data.size()-1;c++)
      data[c] = data[c+1] - data[c];
    LatticeComplexD tmp(data[0].Grid());
    std::cout << "Reading file: "<<file_list.back()<<std::endl;
    readFile(tmp,file_list.back());
    data.back() = tmp - data.back();
  }
  
  /****************************************************************/
  /**************      Determine Frame Dimensions    **************/
  /****************************************************************/
  int coor_map[3] = {0,1,2}, ext_latt_size[latt_size.size()+1];
  for(int d=0,d_c=0,i_o=0; d<3; d++, d_c++){
    if( d_c == dynm_dir ) d_c++;
    for(auto omit_dir:omit_dirs) if( omit_dir == d_c ) d_c++;
    coor_map[d] = d_c; 
  }
  for(int d=0; d<latt_size.size();d++) ext_latt_size[d] = latt_size[d]; ext_latt_size[latt_size.size()] = data.size();

  // DEBUG
  for(auto F: data) std::cout<<"Max: "<<std::sqrt(maxLocalNorm2(F))<<" "<<real(TensorRemove(sum(F)))<<" "<<std::sqrt(maxLocalNorm2(F)/norm2(F)*grid->gSites() )<<std::endl;
  for(int d=0;d<3;d++) std::cout<<d<<" "<<coor_map[d]<<std::endl;
  for(int s : ext_latt_size) std::cout<<s<<std::endl;


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
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  // Set #Total Frames
  int frameCount = ext_latt_size[dynm_dir];// TODO: can take from input which dir is translated
  if ( !mpeg ) frameCount *= ext_latt_size[coor_map[0]];
  if(xlate_omit_dir>=0) frameCount *= ext_latt_size[omit_dirs[xlate_omit_dir]];

  double max_cntr = 0; // max slider value; used in interactive slider widget
  // If file index is not omitted but presented in the frame, all files are displayed in a single frame
  int fc = omit_dirs.back()<latt_size.size()? 1: data.size(); // TODO: can increase #frames corresp. to diff (omit_dir)-intercepts if omit_dir<4
  std::vector<FrameUpdater *> fu_list;
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

    auto nrm    = norm2(data[f]);
    //auto nrmbar = nrm/vol;
    auto rms    = sqrt(nrm/vol);
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
    int skip_one_omit_dir = xlate_omit_dir>=0?1:0;
    std::string info="";
    if(sum_omit_dir) for(auto dir : omit_dirs) info+=dynm_labels[dir];
    else for(int i_d = skip_one_omit_dir; i_d<omit_dirs.size(); i_d++) info+= dynm_labels[omit_dirs[i_d]]+"="+std::to_string(omit_intcpts[i_d])+" ";
    std::string ext = sum_omit_dir? "summed over "+info+"-dir":info;
    std::string txt = omit_dirs.back()<latt_size.size()?"All Files: "+ext : file_list[f];
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
    fu->ext_latt_size = ext_latt_size;
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
	ind_list.push_back(fname.substr(fname.rfind((separator))+4));
	std::cout<<ind_list.back()<<std::endl;
      }
      fu->dynmIndexF = ind_list;
      fu->use_dynmIndexF = true;
    }
    else if(index_list.is_open()){
      assert(!use_fname_as_frame_counter && "file name is not to be used as a frame counter");

      std::string line;
      while(!index_list.eof()){
	getline(index_list,line);
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

    double nf = fc;//file_list.size();
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
