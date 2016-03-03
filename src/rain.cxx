/**
 * Display image passed as argument along with some useful metadata.
 *
 * Requires libraries ITK, VTK, and ItkVtkGlue.
 * 
 * Example based on:
 *    http://itk.org/Wiki/ITK/Examples/SimpleOperations/WidthHeight
 *    http://itk.org/Wiki/ITK/Examples/IO/ImageFileReader
 */

#include "itkImage.h"
#include "itkImageFileReader.h"  // read image file
 
#include "QuickView.h"           // display image
 
int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image< double, 2 >         ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  std::stringstream desc; // holds caption

  // read image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update(); // update image properties

  // get image properties
  // dimensions: [0] corresponds to width, [1] to height
  desc << "width: " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0]
    << "px, height: " << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] 
    << "px";

  // process image
  // TODO

  // display image
  QuickView viewer;
  // add image with caption to viewer; the 2nd & 3rd parameters are optional
  viewer.AddImage<ImageType>( reader->GetOutput(), true, desc.str() );  
  viewer.Visualize();

  return EXIT_SUCCESS;
}
