#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
//segmentation Filters
#include "itkConnectedThresholdImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkIsolatedConnectedImageFilter.h"
//Smoothing Filters
#include "itkMedianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkBilateralImageFilter.h"
//caster
#include "itkCastImageFilter.h"

//define input Image Type to be used throughout the program
const unsigned int    Dimension = 3;

typedef signed short InputPixelType;
typedef itk::Image< InputPixelType,Dimension > InputImageType;

typedef float FloatPixelType;
typedef itk::Image<FloatPixelType, Dimension> FloatImageType;

//declare dicom type
typedef itk::GDCMImageIO ImageIOType;
//Sets the value of the pixels that lie within the thresholded region. 255 --> white
const int gReplaceValue = 255;

class DcmReadWrite{
	char * inputPath; char * outputPath;

	ImageIOType::Pointer gdcmImageIO;
	
    //Read the Image in DICOM
	typedef itk::ImageFileReader< InputImageType > ReaderType;
	ReaderType::Pointer reader;
	 
	//Dicom Image Writer
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writer;     
 
	//constructor takes input and output path of image
	public: DcmReadWrite(char * iInputPath, char * iOutputPath){
		inputPath = (char *) malloc (100*sizeof(char));
		inputPath = iInputPath;
		
		outputPath = (char *) malloc (100*sizeof(char));
		outputPath = iOutputPath;
		std::cout  << inputPath << " " << outputPath << "\n"; 

		//instantiate Dicom Image type
		gdcmImageIO = ImageIOType::New();

		//instantiate reader
		reader = ReaderType::New();
		reader->SetImageIO( gdcmImageIO );

		//instantiate Writer
		writer = WriterType::New();
		writer->SetImageIO(gdcmImageIO);
	}

	//Read function should return InputImageType
	InputImageType::Pointer ReadDicom(){
		 
		reader->SetFileName( inputPath);	
		try
		{
			    reader->Update();
    	}
	    catch (itk::ExceptionObject & e)
    	{
		    std::cerr << "exception in file reader " << std::endl;
		    std::cerr << e << std::endl;
		    //return NULL;
	    }	
		return reader->GetOutput();
	}
	
	//Write Function returns status
	int WriteDicom(InputImageType::Pointer iImage){

		writer->SetFileName( outputPath );
		writer->SetInput(iImage);

  		try
    	{
		    writer->Update();
    	}
		catch (itk::ExceptionObject & e)
    	{
    		std::cerr << "exception in file writer " << std::endl;
		    std::cerr << e << std::endl;
		    return EXIT_FAILURE;
    	}
		return 1;
	}
};


class Segmenter{
	int xSeed, ySeed, zSeed;
	typedef itk::ConnectedThresholdImageFilter< InputImageType,InputImageType > ConnectedFilterType;
	typedef itk::ConfidenceConnectedImageFilter< InputImageType,InputImageType > ConfidenceFilterType;
	typedef itk::NeighborhoodConnectedImageFilter< InputImageType,InputImageType > NeighborhoodConnectedFilterType;
	typedef itk::IsolatedConnectedImageFilter< InputImageType, InputImageType > IsolatedFilterType;
	

	// Connected Threshold
	public:
	InputImageType::Pointer ConnectedThreshold(InputImageType::Pointer iImage, int iXSeed, int iYSeed, int iZSeed, float iLowThreshold, float iHighThreshold){

		//std::cout<<iXSeed<<iYSeed<<iZSeed<<iLowThreshold<<iHighThreshold;	
        
		//Set up the Connected threshold object
		ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		const InputPixelType lowerThreshold =   iLowThreshold ;
	        const InputPixelType upperThreshold =   iHighThreshold ;

		connectedThreshold->SetLower(  lowerThreshold  );
  		connectedThreshold->SetUpper(  upperThreshold  );
        connectedThreshold->SetReplaceValue( gReplaceValue );

        InputImageType::IndexType  index;

        index[0] = iXSeed;
        index[1] = iYSeed;
		index[2] = iZSeed;
 
		connectedThreshold->SetSeed( index );

		//set the input Image
		connectedThreshold->SetInput(iImage);

		//Execute the Connected threshold

		try{
			connectedThreshold->Update();

		}catch(itk::ExceptionObject & e){
			std::cerr << "exception Connected Threshold " << std::endl;
		    std::cerr << e << std::endl;
		}		
  		return connectedThreshold->GetOutput();
	}	

	//Confidence Connected
	public:
		InputImageType::Pointer ConfidenceConnected(InputImageType::Pointer iImage, int iXSeed, int iYSeed, int iZSeed, float iMultiplier, int iIterations	){
			
			//Initialize connected Threshold filter
			ConfidenceFilterType::Pointer confidenceConnected = ConfidenceFilterType::New();
			
			//initialize Connected Threshold Parameters			  
			confidenceConnected->SetMultiplier( iMultiplier );
 		    confidenceConnected->SetNumberOfIterations( iIterations );
            confidenceConnected->SetReplaceValue( gReplaceValue);

			//Set the seed
            InputImageType::IndexType  index;

            index[0] = iXSeed;
            index[1] = iYSeed;
			index[2] = iZSeed;
 
            confidenceConnected->SetSeed( index );  

			//setInput Image
			confidenceConnected->SetInput(iImage);
			
			//Executing the Filter
			try{
				confidenceConnected->Update();

			}catch(itk::ExceptionObject  & e){
				 std::cerr << "exception confidence Conneceted Threshold " << std::endl;
			         std::cerr << e << std::endl;
			}
			
			return confidenceConnected->GetOutput();
		}


	//Neighborhood Connected

	public:
		InputImageType::Pointer NeighborhoodConnected(InputImageType::Pointer iImage, int iXSeed, int iYSeed, int iZSeed, float iLowIntensity, float iHighIntensity, int iRadius){
			
			NeighborhoodConnectedFilterType::Pointer neighborhoodConnected = NeighborhoodConnectedFilterType::New();

			//set threshold value
  			const InputPixelType lowerThreshold =  iLowIntensity;
			const InputPixelType upperThreshold =  iHighIntensity;

			neighborhoodConnected->SetLower(  lowerThreshold  );
			neighborhoodConnected->SetUpper(  upperThreshold  );
			//Set Replace Value
			neighborhoodConnected->SetReplaceValue( gReplaceValue );

			//set Seed
			InputImageType::IndexType  index;

			index[0] = iXSeed;
			index[1] = iYSeed;
			index[2] = iZSeed;
 
			neighborhoodConnected->SetSeed( index );

			InputImageType::SizeType   radius;

			radius[0] = iRadius;   // pixels along X
			radius[1] = iRadius;   // pixels along Y
			radius[2] = iRadius;   // pixels along Z

			neighborhoodConnected->SetRadius( radius );

			//set the image
			neighborhoodConnected->SetInput(iImage);	
	
			try{
				neighborhoodConnected->Update();
			}catch(itk::ExceptionObject & e){
				std::cerr << "exception Neighborhood Connected Threshold " << std::endl;
			    std::cerr << e << std::endl;
			}	

		return neighborhoodConnected->GetOutput();	
		}

	//Isolated Connected 	

	public: 
		InputImageType::Pointer IsolatedConnected(InputImageType::Pointer iImage, int iXSeed1, int iYSeed1, int iZSeed1, int iXSeed2, int iYSeed2, int iZSeed2, float iLowThreshold ){
			
			IsolatedFilterType::Pointer isolatedConnected = IsolatedFilterType::New();

			//initialize the filter with parameters
			
			//Threshold
			const InputPixelType lowerThreshold =  iLowThreshold ;

			//Seed 1
			InputImageType::IndexType  indexSeed1;

  			indexSeed1[0] = iXSeed1;
			indexSeed1[1] = iYSeed1;
			indexSeed1[2] = iZSeed1;

			//seed 2
			InputImageType::IndexType  indexSeed2;

			indexSeed2[0] = iXSeed2;
			indexSeed2[1] = iYSeed2;
			indexSeed2[2] = iZSeed2;

			isolatedConnected->SetLower(  lowerThreshold  );
			isolatedConnected->SetSeed1( indexSeed1 );
			isolatedConnected->SetSeed2( indexSeed2 );
			
			//set replace value
			isolatedConnected->SetReplaceValue( gReplaceValue );

			//set Input Image
			isolatedConnected->SetInput(iImage);

			try
    		{
				isolatedConnected->Update();
    		}
			catch( itk::ExceptionObject & excep )
 			{
				std::cerr << "Isolated Connected Error" << std::endl;
				std::cerr << excep << std::endl;
			}
			return isolatedConnected->GetOutput();
		}
};


class SmoothingFilter{
	typedef itk::MedianImageFilter<InputImageType, InputImageType > MedianFilterType;
	typedef itk::DiscreteGaussianImageFilter<InputImageType, InputImageType>  DiscreteGaussianFilterType;
	typedef itk::CurvatureFlowImageFilter< FloatImageType, FloatImageType > CurvatureFlowFilterType;
        typedef itk::BilateralImageFilter<InputImageType, InputImageType > BilateralFilterType;
	
	//declare caster from short to float
	typedef itk::CastImageFilter< InputImageType, FloatImageType > CastingShortToFloatFilterType;
	typedef itk::CastImageFilter< FloatImageType,InputImageType> CastingFloatToShortFilterType;



	//Bilateral Image Smoothing
	public:
	InputImageType::Pointer Bilateral(InputImageType::Pointer iImage, double iDomainSigma, double iRangeSigma){
		
		//initialize the filter
		BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
		//initialize Parameters
		bilateralFilter->SetDomainSigma (iDomainSigma);
		bilateralFilter->SetRangeSigma(iRangeSigma);	

		//set Input Image
		bilateralFilter->SetInput(iImage);

		try{
			bilateralFilter->Update();
		}catch(itk::ExceptionObject & e){
			std::cerr<<"error in bilateral Filter" <<std::endl;
			std::cerr<<e<<std::endl;
		}

		return bilateralFilter->GetOutput();
	}

	//curvature Flow Filter

	public: 
		InputImageType::Pointer CurvatureFlow(InputImageType::Pointer iImage, int iIterations, double iTimeStep){
			
			//curvature Flow needs Input of type float.
			//since DICOM images come as type signed Short a Caster is needed
			CurvatureFlowFilterType::Pointer curvatureFlowFilter = CurvatureFlowFilterType::New();

			curvatureFlowFilter->SetNumberOfIterations(iIterations);
			curvatureFlowFilter->SetTimeStep(iTimeStep);
			
			//Initialize casting filter from 'unsgined char' to 'float'
			CastingShortToFloatFilterType::Pointer casterShortToFloat = CastingShortToFloatFilterType::New();

			//initialize casting filter from 'float' to 'unsigned char'
			CastingFloatToShortFilterType::Pointer casterFloatToShort = CastingFloatToShortFilterType::New();

			//Set pipeline
			casterShortToFloat->SetInput(iImage);
			curvatureFlowFilter->SetInput(casterShortToFloat->GetOutput());
			casterFloatToShort->SetInput(curvatureFlowFilter->GetOutput());
			
			try{
				casterFloatToShort->Update();
			}catch(itk::ExceptionObject & e){
				std::cerr<<"error in Curvature Flow Filtering"<<std::endl;
				std::cerr<<e<<std::endl;
			}			
		
		return iImage;
	}

	//discrete Gaussian Filter	

	public :
		InputImageType::Pointer DiscreteGaussian(InputImageType::Pointer iImage, double iVariance){
			DiscreteGaussianFilterType::Pointer discreteGaussianFilter = DiscreteGaussianFilterType::New();

			//initialize the filter
			discreteGaussianFilter->SetVariance(iVariance);

			//set inputImage
			discreteGaussianFilter->SetInput(iImage);

			try{
				discreteGaussianFilter->Update();
			}catch(itk::ExceptionObject & e){
				std::cerr<<"Discrete Gaussian Filter Error"<<std::endl;
				std::cerr<<e<<std::endl;
			}
							
		return discreteGaussianFilter->GetOutput();
		}

	//Median Filter

	public:
	InputImageType::Pointer Median(InputImageType::Pointer iImage, int iRadius){

		//initialize Median Filter
		MedianFilterType::Pointer medianFilter = MedianFilterType::New();

		//create radius
		MedianFilterType::InputSizeType radius;
		radius.Fill(iRadius);

		//set radius within filter
		medianFilter->SetRadius(radius);

		//set input Image
		medianFilter->SetInput(iImage);

		try{
			medianFilter->Update();
		}catch(itk::ExceptionObject & e){
			std::cerr<<"Median Filter Error"<<std::endl;
			std::cerr<<e<<std::endl;
		}
		
		return medianFilter->GetOutput();
	}
};


int main( int argc, char * argv[] )
{

	//segmentation variables
	int segmentation;
	int xseed, yseed, zseed, segmentationIteration, xseed2, yseed2, zseed2, neighborhoodRadius;
	float lowThreshold, highThreshold, multiplier;
	
	//smoothing variables
	int smoothing;
	int radius, variance, timestep, smoothingIteration, domainSigma, rangeSigma;

	//Input and output Images
	char * inputImagePath = (char *) malloc(100*sizeof(char)); 
	char * outputImagePath = (char * )malloc(100*sizeof(char));

	//Create object to read and Write Dicom
	std::cout << "enter: InputImagePath OutputImagePath\n";
	std::cin >> inputImagePath >> outputImagePath;	

	DcmReadWrite readerWriter (inputImagePath,outputImagePath);

	//Read the Image specified in the constructor into object InputImageType
	InputImageType::Pointer inputImage = readerWriter.ReadDicom();
	
	//prompt user menu
    std::cout << " 1. Bilateral Smoothing \n 2. Curvature Flow Smoothing\n 3. Discrete Gaussian Smoothing \n 4. Median Smoothing \n 5. None\n --> ";
	std::cin >> smoothing;

	//create smoothingFilter object
	SmoothingFilter smoothingFilter;

	switch(smoothing){
 
		case 1:
			std::cout<<"DomainSigma RangeSigma: \n";
			std::cin>> domainSigma >> rangeSigma;
			inputImage = smoothingFilter.Bilateral(inputImage, domainSigma, rangeSigma);
			break;
		case 2:
			std::cout<<"Iterations Timestep: \n";
			std::cin>> smoothingIteration >> timestep;
			inputImage = smoothingFilter.CurvatureFlow(inputImage,smoothingIteration,timestep);
			break;
		case 3:
			std::cout<<"Variance: \n";
			std::cin>>variance;
			inputImage = smoothingFilter.DiscreteGaussian(inputImage, variance);
			break;
		case 4:	
			std::cout<<"Radius: \n";
			std::cin>>radius;
			inputImage = smoothingFilter.Median(inputImage,radius);
			break;
		default:
			std::cout<< smoothing <<"No Smoothing Algorithm was selected\n";
			break;
	}

	//get the input using getline to collect the '\n' that is left as a result of the previous call to "cin"
	std::string str;
	getline(std::cin, str);

	//prompt user menu
	std::cout << " 1. Connected Threshold \n 2. Confidence Connected\n 3. Neighborhood Connected \n 4. Isolated Connected \n 5. None\n --> ";
	std::cin >> segmentation;

	//object that handles Segmentation
	Segmenter segmenter;

	/*switch Segmentation Filter*/
	switch (segmentation){

		//Connected Threshold
		case 1:
			std::cout<<" XSeed YSeed ZSeed LowThreshold HighThreshold\n";
			std::cin >> xseed >> yseed >> zseed >> lowThreshold >> highThreshold;
			inputImage = segmenter.ConnectedThreshold(inputImage, xseed, yseed, zseed, lowThreshold, highThreshold);		
			break;
		//Confidence Connected
		case 2:
			std::cout<<" XSeed YSeed ZSeed Multiplier Iterations\n";
			std:: cin >> xseed >> yseed >> zseed >> multiplier >> segmentationIteration;
			inputImage = segmenter.ConnectedThreshold(inputImage, xseed, yseed, zseed, multiplier, segmentationIteration);
			break;
		//Neighborhood connected
		case 3:
			std::cout<<" XSeed YSeed ZSeed LowThreshold HighThreshold Radius\n";
			std::cin >> xseed>> yseed >> zseed >> lowThreshold >> highThreshold >> neighborhoodRadius;
			inputImage = segmenter.NeighborhoodConnected(inputImage, xseed, yseed, zseed, lowThreshold, highThreshold, neighborhoodRadius);
			break;
		//Isolated Connected
		case 4:
			std::cout<<" XSeed1  YSeed1 ZSeed1 XSeed2 YSeed2 ZSeed2 LowThreshold\n";
			std:: cin >> xseed >> yseed >> zseed >> xseed2 >> yseed2 >> zseed2 >> lowThreshold;			
			inputImage = segmenter.IsolatedConnected(inputImage, xseed, yseed, zseed, xseed2, yseed2, zseed2, lowThreshold);
			break;
		default:
			std::cout << segmentation <<" No Segmentation Algorithm was selected\n";
			break;
	
	}
	
	//Write processed inputImage into the output image	
    readerWriter.WriteDicom(inputImage);

return 1;
}