// Mean squares 3D affine transform registration
//
// Usage: ./reg input1 input2 output [diff_img1]
//              [diff_img2] [step_length] [max_iterations]
// where input1, input2, output, diff_img1, diff_img2 are RAW image header
// files (*.mhd). Default step length is 0.1 and default max iterations is 300.
// diff_img1: the difference image between the fixed and moving images before
//            registration
// diff_img2: the difference image between the fixed and resampled moving image
//
// Based on examples:
// http://www.na-mic.org/svn/Slicer3-lib-mirrors/trunk/Insight/Testing/Code/IO/itkMetaImageStreamingWriterIOTest.cxx
// http://www.itk.org/Doxygen/html/Examples_2RegistrationITKv4_2ImageRegistration20_8cxx-example.html

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <fstream>
#include "stdint.h"

// image file I/O
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"

// 3d registration
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

//
//  The following piece of code implements an observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

class CommandIterationUpdate : public itk::Command
{
public:
    typedef CommandIterationUpdate Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);
protected:

    CommandIterationUpdate() { };
public:
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef const OptimizerType * OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
        Execute((const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = static_cast<OptimizerPointer> (object);
        if (!itk::IterationEvent().CheckEvent(&event))
        {
            return;
        }
        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << optimizer->GetValue() << "   ";
        std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

int main(int argc, char* argv[])
{
    // display help if not enough command line arguments
    if (argc < 4)
    {
        std::cerr <<
                "Mean squares 3D affine transform registration\n\n"
                "Usage: ./reg input1 input2 output [diff_img1]\n"
                "             [diff_img2] [step_length] [max_iterations]\n"
                "where input1, input2, output, diff_img1, diff_img2 are RAW image header\n"
                "files (*.mhd). Default step length is 0.1 and default max iterations is 300.\n"
                "diff_img1: the difference image between the fixed and moving images before\n"
                "           registration\n"
                "diff_img2: the difference image between the fixed and resampled moving image\n"
                ;
        return EXIT_FAILURE;
    }

    // Define image types
    const unsigned int Dimension = 3;
    typedef int16_t PixelType; // this determines the size of a pixel
    typedef itk::Image< PixelType, Dimension > FixedImageType;
    typedef itk::Image< PixelType, Dimension > MovingImageType;

    //  The transform type is instantiated using the code below. The template
    //  parameters of this class are the representation type of the space
    //  coordinates and the space dimension.
    typedef itk::AffineTransform< double, Dimension > TransformType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType > MetricType;
    typedef itk::LinearInterpolateImageFunction< MovingImageType, double > InterpolatorType;
    typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    RegistrationType::Pointer registration = RegistrationType::New();
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);

    //  The transform object is constructed below and passed to the registration
    //  method.
    TransformType::Pointer transform = TransformType::New();
    registration->SetTransform(transform);

    typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
    itk::MetaImageIO::Pointer fixedImageIO = itk::MetaImageIO::New();
    FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetImageIO(fixedImageIO);
    fixedImageReader->SetFileName(argv[1]);
    fixedImageReader->SetUseStreaming(true);
    fixedImageIO->SetUseStreamedReading(true);

    typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
    itk::MetaImageIO::Pointer movingImageIO = itk::MetaImageIO::New();
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetImageIO(movingImageIO);
    movingImageReader->SetFileName(argv[2]);
    movingImageReader->SetUseStreaming(true);
    movingImageIO->SetUseStreamedReading(true);

    registration->SetFixedImage(fixedImageReader->GetOutput());
    registration->SetMovingImage(movingImageReader->GetOutput());
    fixedImageReader->Update();
    registration->SetFixedImageRegion(fixedImageReader->GetOutput()->GetBufferedRegion());

    // Get & print useful statistics about the input images
    fixedImageReader->GenerateOutputInformation();
    movingImageReader->GenerateOutputInformation();
    FixedImageType::SizeType fixedImageSize;
    MovingImageType::SizeType movingImageSize;
    fixedImageSize = fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    movingImageSize = movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

    std::cout << "dimensions of input images:\n";
    std::cout << "   fixed  [x]:" << fixedImageSize[0] << " [y]:" << fixedImageSize[1]
            << " [z]:" << fixedImageSize[2] << "\n"; // dimensions of input 1
    std::cout << "   moving [x]:" << movingImageSize[0] << " [y]:" << movingImageSize[1]
            << " [z]:" << movingImageSize[2] << "\n"; // dimensions of input 2

    //  In this example, we again use the
    //  CenteredTransformInitializer helper class in order to compute
    //  a reasonable value for the initial center of rotation and the
    //  translation. The initializer is set to use the center of mass of each
    //  image as the initial correspondence correction.
    typedef itk::CenteredTransformInitializer< TransformType, FixedImageType,
            MovingImageType > TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixedImageReader->GetOutput());
    initializer->SetMovingImage(movingImageReader->GetOutput());
    initializer->MomentsOn();
    initializer->InitializeTransform();

    //  Now we pass the parameters of the current transform as the initial
    //  parameters to be used when the registration process starts.
    registration->SetInitialTransformParameters(transform->GetParameters());

    // display the translation parameters
    std::cout << "\nTransformation parameters\n";
    std::cout << "The first (NOutputDimension x NInputDimension) parameters "
            "define the matrix and the\n last NOutputDimension parameters the "
            "translation.\n"; // presumably on (x,y,z) ... ?
    std::cout << transform->GetParameters() << std::endl;

    //  Keeping in mind that the scale of units in scaling, rotation and
    //  translation are quite different, we take advantage of the scaling
    //  functionality provided by the optimizers. We know that the first $N
    //  \times N$ elements of the parameters array correspond to the rotation
    //  matrix factor, and the last $N$ are the components of the translation to
    //  be applied after multiplication with the matrix is performed.
    double translationScale = 1.0 / 1000.0;
    if (argc > 8)
    {
        translationScale = atof(argv[8]);
    }

    typedef OptimizerType::ScalesType OptimizerScalesType;
    OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = 1.0;
    optimizerScales[4] = 1.0;
    optimizerScales[5] = 1.0;
    optimizerScales[6] = 1.0;
    optimizerScales[7] = 1.0;
    optimizerScales[8] = 1.0;
    optimizerScales[9] = translationScale;
    optimizerScales[10] = translationScale;
    optimizerScales[11] = translationScale;
    optimizer->SetScales(optimizerScales);

    //  We also set the usual parameters of the optimization method. In this
    //  case we are using an RegularStepGradientDescentOptimizer. Below, we
    //  define the optimization parameters like initial step length, minimal
    //  step length and number of iterations. These last two act as stopping
    //  criteria for the optimization.
    double steplength = 0.1;
    if (argc > 6)
    {
        steplength = atof(argv[6]);
    }
    unsigned int maxNumberOfIterations = 300;
    if (argc > 7)
    {
        maxNumberOfIterations = atoi(argv[7]);
    }

    optimizer->SetMaximumStepLength(steplength);
    optimizer->SetMinimumStepLength(0.0001);
    optimizer->SetNumberOfIterations(maxNumberOfIterations);
    optimizer->MinimizeOn(); // set the optimizer to do minimization

    // Create the Command observer and register it with the optimizer.
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

    // Trigger execution of the registration method by calling the Update()
    // method.
    try
    {
        std::cout << "Executing registration. This will take awhile."
                << std::endl;
        registration->Update();
        std::cout << "Optimizer stop condition: "
                << registration->GetOptimizer()->GetStopConditionDescription()
                << std::endl;
    } catch (itk::ExceptionObject & err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    // After completion of optimization, recover parameters from the
    // registration method using GetLastTransformParameters().
    OptimizerType::ParametersType finalParameters =
            registration->GetLastTransformParameters();
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double bestValue = optimizer->GetValue();

    // Print results
    std::cout << "Result:" << std::endl;
    std::cout << "   Iterations   = " << numberOfIterations << std::endl;
    std::cout << "   Metric value = " << bestValue << std::endl;

    //  The following code is used to dump output images to files. They
    //  illustrate the final results of the registration. We will resample the
    //  moving image and write out the difference image before and after
    //  registration. We will also rescale the intensities of the difference
    //  images, so that they look better!
    typedef itk::ResampleImageFilter<
            MovingImageType,
            FixedImageType > ResampleFilterType;
    TransformType::Pointer finalTransform = TransformType::New();
    finalTransform->SetParameters(finalParameters);
    finalTransform->SetFixedParameters(transform->GetFixedParameters());
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform(finalTransform);
    resampler->SetInput(movingImageReader->GetOutput());
    FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(fixedImage->GetOrigin());
    resampler->SetOutputSpacing(fixedImage->GetSpacing());
    resampler->SetOutputDirection(fixedImage->GetDirection());
    resampler->SetDefaultPixelValue(100);
    typedef unsigned char OutputPixelType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::CastImageFilter<
            FixedImageType,
            OutputImageType > CastFilterType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    CastFilterType::Pointer caster = CastFilterType::New();
    writer->SetFileName(argv[3]);
    caster->SetInput(resampler->GetOutput());
    writer->SetInput(caster->GetOutput());
    writer->Update();
    typedef itk::SubtractImageFilter<
            FixedImageType,
            FixedImageType,
            FixedImageType > DifferenceFilterType;
    DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
    difference->SetInput1(fixedImageReader->GetOutput());
    difference->SetInput2(resampler->GetOutput());
    WriterType::Pointer writer2 = WriterType::New();
    typedef itk::RescaleIntensityImageFilter<
            FixedImageType,
            OutputImageType > RescalerType;
    RescalerType::Pointer intensityRescaler = RescalerType::New();
    intensityRescaler->SetInput(difference->GetOutput());
    intensityRescaler->SetOutputMinimum(0);
    intensityRescaler->SetOutputMaximum(255);
    writer2->SetInput(intensityRescaler->GetOutput());
    resampler->SetDefaultPixelValue(1);
    // Compute the difference image between the fixed and resampled moving
    // image.
    if (argc > 5)
    {
        writer2->SetFileName(argv[5]);
        writer2->Update();
    }
    typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
    IdentityTransformType::Pointer identity = IdentityTransformType::New();
    // Compute the difference image between the
    // fixed and moving image before registration.
    if (argc > 4)
    {
        resampler->SetTransform(identity);
        writer2->SetFileName(argv[4]);
        writer2->Update();
    }
    return EXIT_SUCCESS;
}
