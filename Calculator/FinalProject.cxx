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
#include "scCalc.hh"

using namespace std;

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

int main(int, char *[])
{

  ScCalc calcualtor;

  double spiral[7][3] = {{0.0, 0.0, 0.0},
                         {1.0, 1.0, 0.0},
                         {0.5, 2.0, 1.0},
                         {0.0, 3.0, 0.0},
                         {1.0, 4.0, 0.0},
                         {0.5, 5.0, 1.0},
                         {0.0, 6.0, 0.0}};

  double spiral2[7][3] = {{0.0, 1.0, 0.0},
                         {1.0, 2.0, 0.0},
                         {0.5, 3.0, 1.0},
                         {0.0, 4.0, 0.0},
                         {1.0, 5.0, 0.0},
                         {0.5, 6.0, 1.0},
                         {0.0, 7.0, 0.0}};

  unsigned spLength = 7;
  double color1[3] = {1, 0, 0};
  double color2[3] = {0, 1, 0};

  double (*ptr)[3] = spiral;

  for (int i = 0; i < 7 - 5; ++i)
  {
    cout << "derivative = " << calcualtor.curveDerivative(ptr);
    cout << " second derivative = " << calcualtor.curveSecDerivative(ptr++);
     cout << endl;
  }
 
  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer = 
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = 
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(makeLine(spiral,spLength,color1));
  renderer->AddActor(makeLine(spiral2,spLength,color2));
 
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}