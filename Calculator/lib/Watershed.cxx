#include <iostream>
// Software Guide : BeginCodeSnippet
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"
// Software Guide : EndCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"
int main( int argc, char *argv[] )
{
  if (argc < 8 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage conductanceTerm diffusionIterations lowerThreshold outputScaleLevel gradientMode " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::RGBPixel< unsigned char >       RGBPixelType;
  typedef itk::Image< RGBPixelType, 3 >        RGBImageType;
  typedef itk::Vector< float, 3 >              VectorPixelType;
  typedef itk::Image< VectorPixelType, 3>     VectorImageType;
  typedef itk::Image< itk::IdentifierType, 3 > LabeledImageType;
  typedef itk::Image< float, 3 >               ScalarImageType;

  typedef itk::ImageFileReader< RGBImageType >   FileReaderType;
  typedef itk::VectorCastImageFilter< RGBImageType, VectorImageType >
                                                 CastFilterType;
  typedef itk::VectorGradientAnisotropicDiffusionImageFilter<
                        VectorImageType, VectorImageType >
                                                 DiffusionFilterType;
  typedef itk::VectorGradientMagnitudeImageFilter< VectorImageType >
                                                 GradientMagnitudeFilterType;
  typedef itk::WatershedImageFilter< ScalarImageType >
                                                 WatershedFilterType;
  typedef itk::ImageFileWriter<RGBImageType> FileWriterType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(argv[1]);
  CastFilterType::Pointer caster = CastFilterType::New();

  DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
  diffusion->SetNumberOfIterations( atoi(argv[4]) );
  diffusion->SetConductanceParameter( atof(argv[3]) );
  diffusion->SetTimeStep(0.125);
  GradientMagnitudeFilterType::Pointer
    gradient = GradientMagnitudeFilterType::New();
  gradient->SetUsePrincipleComponents(atoi(argv[7]));
  WatershedFilterType::Pointer watershed = WatershedFilterType::New();
  watershed->SetLevel( atof(argv[6]) );
  watershed->SetThreshold( atof(argv[5]) );
  typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long>
    ColorMapFunctorType;
  typedef itk::UnaryFunctorImageFilter<LabeledImageType,
    RGBImageType, ColorMapFunctorType> ColorMapFilterType;
  ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName(argv[2]);
  caster->SetInput(reader->GetOutput());
  diffusion->SetInput(caster->GetOutput());
  gradient->SetInput(diffusion->GetOutput());
  watershed->SetInput(gradient->GetOutput());
  colormapper->SetInput(watershed->GetOutput());
  writer->SetInput(colormapper->GetOutput());
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
