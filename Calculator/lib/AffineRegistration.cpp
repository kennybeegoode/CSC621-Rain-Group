#include "AffineRegistration.hh"

using namespace std;

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

/**
 * Performs 3D affine registration using mean squares, with a default max
 * number of iterations of 300.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 */
void AffineRegistration::alignAffine(string fixedImageInput, string movingImageInput,
        double transformParameters[][4])
{
    return AffineRegistration::alignAffine(fixedImageInput, movingImageInput,
            transformParameters, 300);
}

/**
 * Performs 3D affine registration using mean squares.
 * @param fixedImageInput        - a metaimage header
 * @param movingImageInput       - a metaimage header
 * @param transformParameters    - array[4][4] to write transform matrix
 * @param maxNumberOfIterations  - max number of iterations; default is 300
 */
 void AffineRegistration::alignAffine(string fixedImageInput, string movingImageInput,
         double transformParameters[][4], int maxNumberOfIterations)
 {
    //put everything that beliongs to main in here
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
