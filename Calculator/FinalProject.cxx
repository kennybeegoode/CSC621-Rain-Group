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
#include "lib/scCalc.hh"

#include "lib/AffineRegistration.hh"

using namespace std;

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

//This is the main method for the entire project, add your part here
int main(int, char *[])
{
  //SEED INPUT GUI
  //TODO: Ken, your class should go here
  //The output should be a double[3]

  double seed[3] = {0.0, 0.0, 0.0}; //Hardcoded seed to be used while ken is working on his GUI

  //SEGMENTATION
  //TODO: Marie, add a call to your segmentation class here
  //Output should be a double[spine length][3]
  //Feel free to modify hardcoded seed if ken's is not done yet

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
  //TODO: Eric and Monte, change this to produce output!
  //Hardcoded registration output
  // double trans[4][4] = {{1.0, 0.0, 0.0, 5},
  //                           {0.0, 1.0, 0.0, 0.0},
  //                           {0.0, 0.0, 1.0, 0.0},
  //                           {0.0, 0.0, 0.0, 1.0}};
  double trans[4][4]; // to be populated by registration algorithm

  AffineRegistration *reg = new AffineRegistration();
  // run registration with default number of max optimizations (300)
  //reg->alignAffine("case1.mhd", "case2.mhd", trans);
  // test with 1 iteration of optimizer
  reg->alignAffine("case1.mhd", "case2.mhd", trans, 1);

  ScCalc *calculator = new ScCalc();

  //FINAL RESULT CALCULATION
  //TODO: Juris will need to improve this to produce reasonable results
  calculator->loadSpine1(spiral, spLength);
  calculator->loadSpine2(spiral2, spLength2);
  calculator->loadTransofrm(trans);
  calculator->transformSpine1();

  //Set colors for spine
  double color1[3] = {1, 0, 0};
  double color2[3] = {0, 1, 0};

  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(makeLine(calculator->spine1,calculator->spine1Length,color1));
  renderer->AddActor(makeLine(calculator->spine2,calculator->spine2Length,color2));

  //Ouput final view
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
