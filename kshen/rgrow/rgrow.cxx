//file I/O
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//segmentation filters
#include "itkConnectedThresholdImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkIsolatedConnectedImageFilter.h"


typedef unsigned char PixelType;
typedef itk::Image<PixelType, 2> InputImageType;
typedef itk::ConnectedThresholdImageFilter< InputImageType,InputImageType > ConnectedFilterType;

InputImageType::Pointer ConnectedThreshold(InputImageType::Pointer iImage, int iXSeed, int iYSeed, float iLowThreshold, float iHighThreshold){

		//std::cout<<iXSeed<<iYSeed<<iZSeed<<iLowThreshold<<iHighThreshold;	
        
		//Set up the Connected threshold object
		ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		const PixelType lowerThreshold =  iLowThreshold ;
	    const PixelType upperThreshold =  iHighThreshold ;

		connectedThreshold->SetLower(lowerThreshold);
  		connectedThreshold->SetUpper(upperThreshold);
        connectedThreshold->SetReplaceValue(255);

        InputImageType::IndexType  index;

        index[0] = iXSeed;
        index[1] = iYSeed;
		
 
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


int main () {
	int xseed, yseed;
	char * filename = (char *)malloc(100*sizeof(char));

	std::cout << "Enter the filename: \n";
	std::cin >> filename;
	std::cout << "input: xseed, yseed\n";
	std::cin >> xseed;
	std::cin >> yseed;


	//read an image

	itk::ImageFileReader<InputImageType>::Pointer reader = itk::ImageFileReader<InputImageType>::New();

	reader->SetFileName(filename);

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr<<e.GetDescription()<<std::endl;
		return EXIT_FAILURE;
	}

	InputImageType::Pointer image = reader->GetOutput();


	//segmentation

image = ConnectedThreshold(image, xseed, yseed,  0, 50);

//write an image

itk::ImageFileWriter<InputImageType>::Pointer writer = itk::ImageFileWriter<InputImageType>::New();

writer->SetFileName("output.jpg");
writer->SetInput(image);

try
{
	writer->Update();
}
catch (itk::ExceptionObject & e)
{
	std::cerr << e.GetDescription() << std::endl;
	return EXIT_FAILURE;
}
return 1;
}