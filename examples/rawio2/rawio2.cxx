// (1) Read two *.mhd headers then load the corresponding raw image.
// (2) Perform some segmentation (maybe) then use these as the fixed and moving images
// for registration purposes.
// (3) Finally, output the resultant image in a third raw image.


// This is based on example:
// http://www.na-mic.org/svn/Slicer3-lib-mirrors/trunk/Insight/Testing/Code/IO/itkMetaImageStreamingWriterIOTest.cxx 

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <fstream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"
#include "stdint.h"

int main(int argc, char* argv[])
{
  if (argc<4)
  {
    std::cerr << "Usage: " << argv[0] << " input1 input2 output" << std::endl;
    std::cerr << "  where input1, input2, and output are RAW image header"
		<< " files (*.mhd)." << std::endl;
  }
      
  // remove the output file
  itksys::SystemTools::RemoveFile(argv[2]);
    
  typedef int16_t                    PixelType;
  typedef itk::Image<PixelType,3>   ImageType;

  itk::MetaImageIO::Pointer metaImageIO = itk::MetaImageIO::New();

  typedef itk::ImageFileReader<ImageType>         ReaderType;
  typedef itk::ImageFileWriter<ImageType >        WriterType;

  // reader for input 1
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetImageIO(metaImageIO);
  reader1->SetFileName(argv[1]);
  reader1->SetUseStreaming(true);
  metaImageIO->SetUseStreamedReading(true);

  // reader for input 2
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetImageIO(metaImageIO);
  reader2->SetFileName(argv[2]);
  reader2->SetUseStreaming(true);
  metaImageIO->SetUseStreamedReading(true);

  ImageType::RegionType region;
  ImageType::SizeType size;
  ImageType::SizeType fullsize;
  ImageType::IndexType index;
  
  // unsigned int m_NumberOfPieces = 10;

  reader1->GenerateOutputInformation();
  reader2->GenerateOutputInformation();
  size = reader1->GetOutput()->GetLargestPossibleRegion().GetSize();

  index.Fill(0);
  size[0] = fullsize[0];
  size[1] = fullsize[1];
  size[2] = fullsize[2];

  std::cout << "input dimensions\n";
  std::cout << "x: size[0] = " << size[0] << "\n";
  std::cout << "y: size[1] = " << size[1] << "\n";
  std::cout << "z: size[2] = " << size[2] << " (# layers)\n";

  // unsigned int zsize = fullsize[2]/m_NumberOfPieces;

  // Setup the writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[3]);
  
  // for(unsigned int i=0;i<m_NumberOfPieces;i++)
  //   {
  //   std::cout << "Reading piece " << i+1 << " of " << m_NumberOfPieces << std::endl;

  //   index[2] += size[2];

  //   // At the end we need to adjust the size to make sure
  //   // we are reading everything
  //   if(i == m_NumberOfPieces-1)
  //     {
  //     size[2] = fullsize[2]-index[2];
  //     }
  //   else
  //     {
  //     size[2] = zsize;
  //     }

  //   region.SetIndex(index);
  //   region.SetSize(size);
  //   reader->GetOutput()->SetRequestedRegion(region);
  //   try
  //     {
  //     reader->Update();
  //     }
  //   catch (itk::ExceptionObject &ex)
  //     {
  //     std::cout << "ERROR : " << ex << std::endl;
  //     return EXIT_FAILURE;
  //     }
   
    // Write the image     
    itk::ImageIORegion  ioregion(3);
    itk::ImageIORegion::IndexType index2;
    index2.push_back(region.GetIndex()[0]);
    index2.push_back(region.GetIndex()[1]);
    index2.push_back(region.GetIndex()[2]);
    ioregion.SetIndex(index2);
    itk::ImageIORegion::SizeType size2;
    size2.push_back(region.GetSize()[0]);
    size2.push_back(region.GetSize()[1]);
    size2.push_back(region.GetSize()[2]);
    ioregion.SetSize(size2);
    writer->SetIORegion(ioregion);
    writer->SetInput(reader->GetOutput());
    
    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  //   } // end for pieces
   
   // writer->SetInput(reader->GetOutput());
   // write->update();
     
  return EXIT_SUCCESS;
}
