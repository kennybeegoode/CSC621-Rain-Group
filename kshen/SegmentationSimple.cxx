//image I/O
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"

//segmentation Filters
#include "itkConnectedThresholdImageFilter.h"


const unsigned int    Dimension = 3;

typedef signed short InputPixelType;
typedef itk::Image< InputPixelType,Dimension > InputImageType;

typedef float FloatPixelType;
typedef itk::Image<FloatPixelType, Dimension> FloatImageType;

//declare dicom type
typedef itk::GDCMImageIO ImageIOType;
//Sets the value of the pixels that lie within the thresholded region. 255 --> white
const int gReplaceValue = 255;


int main( int argc, char * argv[] )
{


	//Input and output dir
	char * inputImagePath = (char *) malloc(100*sizeof(char)); 
	char * outputImagePath = (char * )malloc(100*sizeof(char));

	//Create object to read and Write Dicom
	std::cout << "enter: InputImagePath OutputImagePath\n";
	std::cin >> inputImagePath >> outputImagePath;	

	DcmReadWrite readerWriter (inputImagePath,outputImagePath);

	//Read the Image specified in the constructor into object InputImageType
	InputImageType::Pointer inputImage = readerWriter.ReadDicom();

	//seed coordinates
	int xseed, yseed, zseed;
	int lowThreshold, highThreshold;

	//manunally plant seed and set thresholds
	std::cout<<" XSeed YSeed ZSeed LowThreshold HighThreshold\n";
	std::cin >> xseed >> yseed >> zseed >> lowThreshold >> highThreshold;

	//calling connectedThreshold segmentation filter
	inputImage = segmenter.ConnectedThreshold(inputImage, xseed, yseed, zseed, lowThreshold, highThreshold);

	//Write processed inputImage into the output image	
    readerWriter.WriteDicom(inputImage);

    return 1;
}






class Segmenter{
	int xSeed, ySeed, zSeed;
	typedef itk::ConnectedThresholdImageFilter< InputImageType,InputImageType > ConnectedFilterType;
	typedef itk::ConfidenceConnectedImageFilter< InputImageType,InputImageType > ConfidenceFilterType;
	typedef itk::NeighborhoodConnectedImageFilter< InputImageType,InputImageType > NeighborhoodConnectedFilterType;
	typedef itk::IsolatedConnectedImageFilter< InputImageType, InputImageType > IsolatedFilterType;
	

	// Connected Threshold
	public:
	InputImageType::Pointer ConnectedThreshold(InputImageType::Pointer iImage, int iXSeed, int iYSeed, int iZSeed, float iLowThreshold, float iHighThreshold){


	//debug log
	std::cout<<iXSeed<<iYSeed<<iZSeed<<iLowThreshold<<iHighThreshold;	
        
		//Set up the Connected threshold object
		ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		const InputPixelType lowerThreshold =   iLowThreshold ;
	    const InputPixelType upperThreshold =   iHighThreshold ;

		connectedThreshold->SetLower(  lowerThreshold  );
  		connectedThreshold->SetUpper(  upperThreshold  );
        connectedThreshold->SetReplaceValue( gReplaceValue );

        InputImageType::IndexType  index;

        //SetSeed function only takes an array pointer
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
}





// image I/O class

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