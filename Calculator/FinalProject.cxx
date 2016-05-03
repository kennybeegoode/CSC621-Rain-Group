#include <vtkSmartPointer.h>
#include <vtkCurvatures.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkImageViewer2.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMetaImageReader.h>
#include <vtkObjectFactory.h>
#include <vtkResliceImageViewer.h>

#include "lib/scCalc.hh"
#include "lib/Registration.hh"
#include "lib/RegionGrowingNoThreshold.h"

using namespace std;

//seed coords
double seedX, seedY, seedZ;
//Helper funcion, ignore this
vtkSmartPointer<vtkActor> makeLine(double data[][3], unsigned length, double color[3])
{
  vtkSmartPointer<vtkPoints> points =
  vtkSmartPointer<vtkPoints>::New();

  for (int i = 0; i < length; ++i)
  {
    points->InsertPoint(i, data[i]);
  }

  vtkSmartPointer<vtkKochanekSpline> xSpline =
  vtkSmartPointer<vtkKochanekSpline>::New();
  vtkSmartPointer<vtkKochanekSpline> ySpline =
  vtkSmartPointer<vtkKochanekSpline>::New();
  vtkSmartPointer<vtkKochanekSpline> zSpline =
  vtkSmartPointer<vtkKochanekSpline>::New();

  vtkSmartPointer<vtkParametricSpline> spline =
  vtkSmartPointer<vtkParametricSpline>::New();
  spline->SetXSpline(xSpline);
  spline->SetYSpline(ySpline);
  spline->SetZSpline(zSpline);
  spline->SetPoints(points);

  vtkSmartPointer<vtkParametricFunctionSource> functionSource =
  vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource->SetParametricFunction(spline);
  functionSource->Update();

  // Setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper =
  vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(functionSource->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(color[0], color[1], color[2]); //(R,G,B)

  return actor;
}

//mouse event handler
class MouseInteractorStyle3 : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle3* New();

  virtual void OnLeftButtonDown() 
  {
    seedX = this->Interactor->GetEventPosition()[0];
    seedY = this->Interactor->GetEventPosition()[1];
    std::cout << "Seed Choosen at: " << seedX << " " << seedY << std::endl;
    std::cout << "Please close the window to continue...\n";
  }

};

//for mouse event
vtkStandardNewMacro(MouseInteractorStyle3);

//This is the main method for the entire project, add your part here
int main(int argc, char *argv[])
{
  ScCalc *calculator = new ScCalc();
  double* spacing1;
  double* spacing2;
  //SEED INPUT GUI
  //TODO: Ken, your class should go here
  //The output should be a double[3]
  if(argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " InputFile1\n" << "InputFile2\n";
    return EXIT_FAILURE;
  }

  //read input mhd file
  vtkSmartPointer<vtkMetaImageReader>reader =
  vtkSmartPointer<vtkMetaImageReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  spacing1 = reader->GetPixelSpacing();

  //display
  vtkSmartPointer<vtkResliceImageViewer> imageViewer =
  vtkSmartPointer<vtkResliceImageViewer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  imageViewer->SetInputConnection(reader->GetOutputPort());
  
  imageViewer->SetupInteractor(renderWindowInteractor2);
  imageViewer->SetColorLevel(500);
  imageViewer->SetColorWindow(2000);

  //set z coord always the most center slice 
  seedZ = (imageViewer->GetSliceMax())/2;

  //imageViewer->SetSize(512,512);
  imageViewer->SetSlice(seedZ);
  imageViewer->SetSliceOrientationToXY();
  imageViewer->Render();

  vtkSmartPointer<MouseInteractorStyle3> style =
  vtkSmartPointer<MouseInteractorStyle3>::New();
  
  renderWindowInteractor2->SetInteractorStyle(style);
  //renderWindowInteractor2->UpdateSize(100,100);
  renderWindowInteractor2->Start();

  seedZ = imageViewer->GetSlice();

  double seed1[3] = {seedX, seedY, seedZ};

  //////////////////////////////////////////////////////////////////////
  //read input mhd file
  vtkSmartPointer<vtkMetaImageReader>reader1 =
  vtkSmartPointer<vtkMetaImageReader>::New();
  reader1->SetFileName(argv[2]);
  reader1->Update();
  spacing2 = reader1->GetPixelSpacing();

  //display
  vtkSmartPointer<vtkResliceImageViewer> imageViewer1 =
  vtkSmartPointer<vtkResliceImageViewer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor3 = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

  imageViewer1->SetInputConnection(reader1->GetOutputPort());
  imageViewer1->SetupInteractor(renderWindowInteractor3);
  imageViewer1->SetColorLevel(500);
  imageViewer1->SetColorWindow(2000);

  //set z coord always the most center slice 
  seedZ = (imageViewer1->GetSliceMax())/2;

  imageViewer1->SetSlice(seedZ);
  imageViewer1->SetSliceOrientationToXY();
  imageViewer1->Render();
 
  vtkSmartPointer<MouseInteractorStyle3> style1 =
  vtkSmartPointer<MouseInteractorStyle3>::New();
 
  renderWindowInteractor3->SetInteractorStyle(style1);
  renderWindowInteractor3->Start();

  seedZ = imageViewer1->GetSlice();
  
  double seed2[3] = {seedX, seedY, seedZ};

  ////////////////////////////////////////////////////////////////////////////////////////////
  
  //debug log
  std::cout << "Seed1 Set at: " << seed1[0] <<" "<< seed1[1] <<" "<<seed1[2] <<endl;
  std::cout << "Seed2 Set at: " << seed2[0] <<" "<< seed2[1] <<" "<<seed2[2] <<endl;
  std::cout << "Spacing1 Set at: " << spacing1[0] << " " << spacing1[1] << " " << spacing1[2] <<endl;
  std::cout << "Spacing2 Set at: " << spacing2[0] << " " << spacing2[1] << " " << spacing2[2] <<endl;

  //Hardcoded seed to be used while ken is working on his GUI

  //SEGMENTATION
  //TODO: Marie, add a call to your segmentation class here
  //Output should be a double[spine length][3]
  //Feel free to modify hardcoded seed if ken's is not done yet

  RegionGrowingNoThreshold region_growing;
  std::vector<std::vector<int>> centroids1 = region_growing.GetCentroids(argv[1], seed1[0], seed1[1], seed1[2]);
  calculator->printVector(centroids1);
  std::vector<std::vector<int>> centroids2 = region_growing.GetCentroids(argv[2], seed2[0], seed2[1], seed2[2]);
  calculator->printVector(centroids2);
  
  //Hardcoded segmentation output
  double spiral[7][3] = {{0.0, 0.0, 0.0},
  {1.0, 1.0, 0.0},
  {0.5, 2.0, 1.0},
  {0.0, 3.0, 0.0},
  {1.0, 4.0, 0.0},
  {0.5, 5.0, 1.0},
  {0.0, 6.0, 0.0}};
  unsigned spLength = 7;

  double spiral2[7][3] = {{0.0, 1.0, 0.0},
  {1.0, 2.0, 0.0},
  {0.5, 3.0, 1.0},
  {0.0, 4.0, 0.0},
  {1.0, 5.0, 0.0},
  {0.5, 6.0, 1.0},
  {0.0, 7.0, 0.0}};
  unsigned spLength2 = 7;

  //REGISTRATION
  Registration *reg = new Registration();

  double trans[4][4]; // to be populated by registration algorithm
  // test with 1 iteration of optimizer
  reg->rigidAlign(argv[1], argv[2], trans, 1);

  double trans2[4][4];
  // test with 1 iteration of optimizer
  //reg->affineAlign(argv[1], argv[2], trans2, 1);

  //FINAL RESULT CALCULATION
  //TODO: Juris will need to improve this to produce reasonable results
  // calculator->spacing1 = spacing1;
  // calculator->spacing2 = spacing2;
  // calculator->loadSpine1(centroids1);
  // calculator->loadSpine2(centroids2);
  // //calculator->loadTransofrm(trans2);
  // //calculator->transformSpine1();

  // //Set colors for spine
  // double color1[3] = {1, 0, 0};
  // double color2[3] = {0, 1, 0};

  // // Setup render window, renderer, and interactor
  // vtkSmartPointer<vtkRenderer> renderer =
  // vtkSmartPointer<vtkRenderer>::New();
  // vtkSmartPointer<vtkRenderWindow> renderWindow =
  // vtkSmartPointer<vtkRenderWindow>::New();
  // renderWindow->AddRenderer(renderer);
  // vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  // vtkSmartPointer<vtkRenderWindowInteractor>::New();
  // renderWindowInteractor->SetRenderWindow(renderWindow);
  // renderer->AddActor(makeLine(calculator->spine1,calculator->spine1Length,color1));
  // renderer->AddActor(makeLine(calculator->spine2,calculator->spine2Length,color2));

  // //Ouput final view
  // renderWindow->Render();
  // renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
