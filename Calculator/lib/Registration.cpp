#include "Registration.hh" 

using namespace std;

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
  typedef   const OptimizerType *                             OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
    Execute( (const itk::Object *)caller, event);
    }
  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

void Registration::rigidAlign(string fixedImageInput, string movingImageInput, double transformParameters[][4], int maxNumberOfIterations) {

    std::cout << "Starting Rigid Alignment now with iterations:" << maxNumberOfIterations <<std::endl;
    
    const unsigned int                          Dimension = 3;
  typedef  float                              PixelType;
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::VersorRigid3DTransform< double > TransformType;

  typedef itk::RegularStepGradientDescentOptimizerv4<double>    OptimizerType;
  typedef itk::MeanSquaresImageToImageMetricv4<
                                            FixedImageType,
                                            MovingImageType >   MetricType;
  typedef itk::ImageRegistrationMethodv4<
                                      FixedImageType,
                                      MovingImageType,
                                      TransformType >           RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );

  TransformType::Pointer  initialTransform = TransformType::New();

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  fixedImageInput );
  movingImageReader->SetFileName( movingImageInput );

  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );

  typedef itk::CenteredTransformInitializer<
    TransformType,
    FixedImageType,
    MovingImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer =
    TransformInitializerType::New();
 
  initializer->SetTransform(   initialTransform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( movingImageReader->GetOutput() );
  
  initializer->MomentsOn();
  
  initializer->InitializeTransform();
  
  typedef TransformType::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  VersorType     rotation;
  VectorType     axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  const double angle = 0;
  rotation.Set(  axis, angle  );
  initialTransform->SetRotation( rotation );
 
  registration->SetInitialTransform( initialTransform );

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( initialTransform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetNumberOfIterations( 1 );
  optimizer->SetLearningRate( 0.2 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetReturnBestParametersAndValue(true);

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  const unsigned int numberOfLevels = 1;

  RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( 1 );
  shrinkFactorsPerLevel[0] = 1;

  RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( 1 );
  smoothingSigmasPerLevel[0] = 0;

  registration->SetNumberOfLevels( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  try
    {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;
    }

  const TransformType::ParametersType finalParameters =
                            registration->GetOutput()->Get()->GetParameters();

  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetFixedParameters( registration->GetOutput()->Get()->GetFixedParameters() );
  finalTransform->SetParameters( finalParameters );

  TransformType::MatrixType matrix = finalTransform->GetMatrix();
  TransformType::OffsetType offset = finalTransform->GetOffset();
  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );

  typedef  unsigned char                                          OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
  typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >                 WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( "result.mhd" );

  std::cout << "Set the resulting image as result.mhd" <<std::endl;

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  typedef itk::SubtractImageFilter<
                                  FixedImageType,
                                  FixedImageType,
                                  FixedImageType > DifferenceFilterType;
  DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

  typedef itk::RescaleIntensityImageFilter<
                                  FixedImageType,
                                  OutputImageType >   RescalerType;
  RescalerType::Pointer intensityRescaler = RescalerType::New();

  intensityRescaler->SetInput( difference->GetOutput() );
  intensityRescaler->SetOutputMinimum(   0 );
  intensityRescaler->SetOutputMaximum( 255 );

  difference->SetInput1( fixedImageReader->GetOutput() );
  difference->SetInput2( resampler->GetOutput() );

  resampler->SetDefaultPixelValue( 1 );



}

/**
 * Performs 3D affine registration using mean squares, with a default max
 * number of iterations of 300.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 */
void Registration::affineAlign(string fixedImageInput, string movingImageInput,
        double transformParameters[][4])
{
    return Registration::affineAlign(fixedImageInput, movingImageInput,
            transformParameters, 1);
}

/**
 * Performs 3D affine registration using mean squares.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 * @param maxNumberOfIterations  - max number of iterations; default is 300
 */
 void Registration::affineAlign(string fixedImageInput, string movingImageInput,
         double transformParameters[][4], int maxNumberOfIterations)
 {
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

	// Initialize readers for the two input files
	typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
    itk::MetaImageIO::Pointer fixedImageIO = itk::MetaImageIO::New();
    FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetImageIO(fixedImageIO);
    fixedImageReader->SetFileName(fixedImageInput);
    fixedImageReader->SetUseStreaming(true);
    fixedImageIO->SetUseStreamedReading(true);

    typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
    itk::MetaImageIO::Pointer movingImageIO = itk::MetaImageIO::New();
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetImageIO(movingImageIO);
    movingImageReader->SetFileName(movingImageInput);
    movingImageReader->SetUseStreaming(true);
    movingImageIO->SetUseStreamedReading(true);


	// configure registration algorithm
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

    //  Keeping in mind that the scale of units in scaling, rotation and
    //  translation are quite different, we take advantage of the scaling
    //  functionality provided by the optimizers. We know that the first $N
    //  \times N$ elements of the parameters array correspond to the rotation
    //  matrix factor, and the last $N$ are the components of the translation to
    //  be applied after multiplication with the matrix is performed.
    double translationScale = 1.0 / 1000.0;

	// configure optimizer
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

  optimizer->SetMaximumStepLength(steplength);
    optimizer->SetMinimumStepLength(0.0001);
    optimizer->SetNumberOfIterations(maxNumberOfIterations);
    optimizer->MinimizeOn(); // set the optimizer to do minimization

    // Create the Command observer and register it with the optimizer.
	// TODO: We might not need this.
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

	// Trigger execution of the registration method by calling the Update()
    // method.
    try
    {
        std::cout << "\nExecuting registration. This will take awhile."
                << std::endl;
        registration->Update();
        std::cout << "\nFinished registration.\nOptimizer stop condition: "
                << registration->GetOptimizer()->GetStopConditionDescription()
                << std::endl;
    } catch (itk::ExceptionObject & err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return;
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

    // interpret final transformation parameters as 4x4 matrix
    // tokenize the parameters
    stringstream ss;
    ss << finalParameters << endl;
    string paramsString = ss.str();
    paramsString = paramsString.substr(1, paramsString.length() - 3);
    vector<double> params_v;
    char delim = ',';
    size_t start = paramsString.find_first_not_of(delim), end = start;
    while (start != string::npos)
    {
        end = paramsString.find(delim, start);
        params_v.push_back(atof(paramsString.substr(start, end - start).c_str()));
        start = paramsString.find_first_not_of(delim, end);
    }

    // populate the translation matrix
    transformParameters[0][0] = params_v[0];
    transformParameters[1][0] = params_v[1];
    transformParameters[2][0] = params_v[2];

    transformParameters[0][1] = params_v[3];
    transformParameters[1][1] = params_v[4];
    transformParameters[2][1] = params_v[5];

    transformParameters[0][2] = params_v[6];
    transformParameters[1][2] = params_v[7];
    transformParameters[2][2] = params_v[8];

    transformParameters[3][0] = params_v[9];
    transformParameters[3][1] = params_v[10];
    transformParameters[3][2] = params_v[11];

    transformParameters[0][3] = transformParameters[1][3]
            = transformParameters[2][3] = 0.0;
    transformParameters[3][3] = 1.0;

    // display the matrix
    cout << "   Transform matrix:" << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            cout << setw(15) << transformParameters[j][i] << " ";
        cout << endl;
    }
}
