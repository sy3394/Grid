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

#include <Grid/Grid.h>

#define USE_FLYING_EDGES
#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
typedef vtkFlyingEdges3D isosurface;
#else
#include <vtkMarchingCubes.h>
typedef vtkMarchingCubes isosurface;
#endif

int mpeg = 0 ;
int Ls = -1;
int xlate = 0 ;
int take_diff = 0;
int omit_dir = 4;
int dynm_dir = 3;
std::vector<std::string> dynm_labels = {"X", "Y", "Z", "T", "tau"};

template <class T> void readFile(T& out, std::string const fname){
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

class FrameUpdater : public vtkCallbackCommand
{
public:

  FrameUpdater() {
    TimerCount = 0;
    xoff       = 0;
    x3         = 0;
    imageData  = nullptr;
    grid_data.clear();
    frame_size = nullptr;
    coor_map   = nullptr;
    timerId    = 0;
    maxCount   = -1;
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
	  
	  // Make a new frame
	  int dims[5];
	  auto latt_size = grid_data[0]->Grid()->GlobalDimensions();
	  for(int d=0; d<latt_size.size(); d++) dims[d] = latt_size[d]; dims[4] = grid_data.size();
									  
	  for(int x0=0;x0<frame_size[0];x0++){
	    for(int x1=0;x1<frame_size[1];x1++){
	      for(int x2=0;x2<frame_size[2];x2++){
		Coordinate site({0,0,0,0}); // can take (omit_dir)-intercepts from input
		// the first two frame dim's are always latt dim
		site[coor_map[0]] = (x0+xoff)%frame_size[0]; site[coor_map[1]] = x1;
		RealD value;
		if(dynm_dir<4){ 
		  site[dynm_dir] = x3;
		  if(coor_map[2]<4) site[coor_map[2]] = x2;
		  value = coor_map[2] == 4? real(peekSite(*grid_data[x2],site)) : real(peekSite(*grid_data[0],site));
		}
		else {
		  site[coor_map[2]] = x2;
		  value = real(peekSite(*grid_data[x3],site));
		}
		imageData->SetScalarComponentFromDouble(x0,x1,x2,0,value);
	  }}}

	  if ( xlate ) { 
	    xoff = (xoff + 1)%frame_size[0];
	    if ( xoff== 0 ) x3 = (x3+1)%dims[dynm_dir];
	  } else {
	    x3 = (x3+1)%dims[dynm_dir];
	    if ( x3== 0 ) 	xoff = (xoff + 1)%frame_size[0];
	  }

	  snprintf(text_string,max,"%s=%d",dynm_labels[dynm_dir].c_str(),x3);
	  text->SetInput(text_string);
      
	  std::cout << this->TimerCount<<"/"<<maxCount<< " xoff "<<xoff<<" t_updated "<<x3  <<std::endl;
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
  int* frame_size;
  int* coor_map;
  vtkImageData* imageData = nullptr;
  vtkTextActor* text = nullptr;
  vtkFFMPEGWriter *writer = nullptr;
  int timerId ;
  int maxCount ;
  double rms;
  isosurface * posExtractor;
  isosurface * negExtractor;
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
  std::string arg, save_fname;
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
    save_fname = arg;
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
  std::vector<LatticeComplexD> data2(file_list2.size()-take_diff,grid);
  for(int c=0;c<data2.size();c++) {
    std::cout << "Reading file2: "<<file_list2[c]<<std::endl;
    LatticeComplexD tmp(gridF);
    readFile(tmp,file_list2[c]);
    LatticeComplexD tmp4D(grid); tmp4D = Zero(); data2[c] = Zero();
    for(int i=0; i<Ls;i++){
      ExtractSlice(tmp4D,tmp,i,0);
      data2[c] = data2[c] + tmp4D;
    }

    std::cout<<"Sum "<<c<<" "<<real(TensorRemove(sum(data2[c])))<<std::endl;
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
  if(sum_all){ // TODO: turn save_fname into a vector!!!!!!!!!!!!!
    if(!data1.empty()){
      LatticeComplexD tmp(data1[0].Grid()); tmp = Zero();
      for(int c=0;c<data1.size();c++) tmp = tmp + data1[c];
      writeFile(tmp,save_fname);
    }
    if(!data2.empty()){
      LatticeComplexD tmp(data2[0].Grid()); tmp = Zero();
      for(int c=0;c<data2.size();c++) tmp = tmp + data2[c];
      writeFile(tmp,save_fname);
    }
  }
  
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
	      //0.00010128755 0 0 1048576 0
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

  /****************  Compute Inner Product of 5D Eigen and Gluonic TCD   ************************************/


  /*
    // not for laptop
  for(auto F: data){
    for(int xi=0; xi<F.Grid()->gSites()*data.size(); xi++){
      Coordinate site({xi/(latt_size[1]*latt_size[2]*latt_size[3]),xi/(latt_size[2]*latt_size[3]),xi/latt_size[3],xi%latt_size[3]});
      data_f << real(TensorRemove(peekSite(F,site))) << std::endl;
    }
  }
  */
  /*
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

  int frameCount = dynm_dir<4?latt_size[dynm_dir]:data.size();//can take from input which dir is translated
  if ( !mpeg ) frameCount *= frame_size[0];

  int fc = omit_dir<4? 1: data.size(); // can increase #frames corresp. to diff (omit_dir)-intercepts if omit_dir<4
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
    
    double vol = data[f].Grid()->gSites();

    auto max_norm= std::sqrt(maxLocalNorm2(data[f]));
    auto nrm    = norm2(data[f]);
    auto nrmbar = nrm/vol;
    auto rms    = sqrt(nrmbar);
    std::cout<<"ratio: "<<max_norm/rms<<std::endl;
    double contour = default_contour * rms; // default to 1 x RMS

    // The following reader is used to read a series of 2D slices (images)
    // that compose the volume. The slice dimensions are set, and the
    // pixel spacing. The data Endianness must also be specified. The reader
    // uses the FilePrefix in combination with the slice number to construct
    // filenames using the format FilePrefix.%d. (In this case the FilePrefix
    // is the root name of the file: quarter.)
    vtkNew<vtkImageData> imageData;
    imageData->SetDimensions(latt_size[0],latt_size[1],latt_size[2]);
    imageData->AllocateScalars(VTK_DOUBLE, 1);
    for(int x0=0;x0<frame_size[0];x0++){
      for(int x1=0;x1<frame_size[1];x1++){
	for(int x2=0;x2<frame_size[2];x2++){
	  Coordinate site({0,0,0,0}); // can take (omit_dir)-intercepts from input
	  // the first two frame dim's are always latt dim
	  site[coor_map[0]] = x0; site[coor_map[1]] = x1;
	  if(coor_map[2]<4) site[coor_map[2]] = x2;
	  RealD value = (coor_map[2] == 4) ? real(peekSite(data[0],site)) : real(peekSite(data[f],site));
	  imageData->SetScalarComponentFromDouble(x0,x1,x2,0,value);
    }}}

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

    std::string txt = omit_dir<4?"All Files Displayed" : file_list[f];
    if(take_diff) txt = "Diff: Next File - "+txt+" )";
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
    if(omit_dir<4) for(auto e: data) tmp.push_back(&e);
    else tmp.push_back(&data[f]);
    vtkNew<FrameUpdater> fu;
    fu->imageData = imageData;
    fu->grid_data = tmp;
    fu->frame_size= frame_size;
    fu->coor_map  = coor_map;
    fu->text      = TextT;
    fu->maxCount = frameCount;
    fu->posExtractor = posExtractor;
    fu->negExtractor = negExtractor;
    fu->rms = rms;
      
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
  renWin->SetSize(1024*file_list.size(), 1024);
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
  
    // Add control of contour threshold
    // Create a slider widget
    vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
    sliderRep->SetMinimumValue(0.1);
    sliderRep->SetMaximumValue(15.0);
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

    double nf = file_list.size();
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
  */
  data_f.close();
  Grid_finalize();

  return EXIT_SUCCESS;
}
