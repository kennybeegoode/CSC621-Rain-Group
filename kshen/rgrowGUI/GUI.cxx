#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkMetaImageReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkPointPicker.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkImageViewer2.h>
#include <iostream>


// Define interaction style
class MouseInteractorStyle3 : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStyle3* New();

  virtual void OnLeftButtonDown() 
  {
    std::cout << "Pressed left mouse button." << std::endl;
    int x = this->Interactor->GetEventPosition()[0];
    int y = this->Interactor->GetEventPosition()[1];
    std::cout << "(x,y) = (" << x << "," << y << ")" << std::endl;

    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }

};

vtkStandardNewMacro(MouseInteractorStyle3);

int main(int argc, char *argv[])
{

  if(argc < 2)
  {
    std::cerr << "Required arguments: image.mha" << std::endl;
    return EXIT_FAILURE;
  }

  int z;

  

  //read mhd image
  vtkSmartPointer<vtkMetaImageReader>reader =

  vtkSmartPointer<vtkMetaImageReader>::New();

  reader->SetFileName(argv[1]);

  reader->Update();


 

  vtkSmartPointer<vtkImageViewer2> imageViewer =

  vtkSmartPointer<vtkImageViewer2>::New();

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
  vtkSmartPointer<vtkRenderWindowInteractor>::New();

  imageViewer->SetInputConnection(reader->GetOutputPort());
  
  imageViewer->SetupInteractor(renderWindowInteractor);

  imageViewer->SetColorLevel(500);

  imageViewer->SetColorWindow(2000);

  imageViewer->SetSlice(50);

  imageViewer->SetSliceOrientationToXY();

  imageViewer->Render();

  renderWindowInteractor->UpdateSize(500,500);

  vtkSmartPointer<MouseInteractorStyle3> style =
  vtkSmartPointer<MouseInteractorStyle3>::New();
  renderWindowInteractor->SetInteractorStyle( style );

  renderWindowInteractor->Initialize();


  renderWindowInteractor->Start();

  return  EXIT_SUCCESS ;
}