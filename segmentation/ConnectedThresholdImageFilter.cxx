/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRGBPixel.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char *argv[])
{
   if(argc < 8) {
      std::cerr << "Missing Parameters " << std::endl;
      std::cerr << "Usage: " << argv[0];
      std::cerr << " inputImage  outputImage seedX seedY seedZ lowerThreshold upperThreshold" << std::endl;
      return EXIT_FAILURE;
    }
    typedef   float           InternalPixelType;
    const     unsigned int    Dimension = 3;
    typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
    typedef unsigned char                            OutputPixelType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::CastImageFilter< InternalImageType, OutputImageType >
                                                     CastingFilterType;

    CastingFilterType::Pointer caster = CastingFilterType::New();

    typedef  itk::ImageFileReader< InternalImageType > ReaderType;
    typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(argv[1]);
    writer->SetFileName(argv[2]);

    // Smoothing the image before appplying the region growing.
    typedef itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType >
       CurvatureFlowImageFilterType;
    CurvatureFlowImageFilterType::Pointer smoothing =
       CurvatureFlowImageFilterType::New();
    smoothing->SetInput(reader->GetOutput());

    typedef itk::ConnectedThresholdImageFilter< InternalImageType,
        InternalImageType > ConnectedFilterType;
    ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
    connectedThreshold->SetInput(smoothing->GetOutput());
    caster->SetInput(connectedThreshold->GetOutput());
    writer->SetInput(caster->GetOutput());
    smoothing->SetNumberOfIterations(2);
    smoothing->SetTimeStep(0.05);

    // Setting the thresholds.
    const InternalPixelType lowerThreshold = atof(argv[6]);
    const InternalPixelType upperThreshold = atof(argv[7]);
    connectedThreshold->SetLower(lowerThreshold);
    connectedThreshold->SetUpper(upperThreshold);

    // The color we aregoing to paint the region with.
    connectedThreshold->SetReplaceValue(255);

    // Setting up seeds.
    InternalImageType::IndexType  index;
    index[0] = atoi(argv[3]);
    index[1] = atoi(argv[4]);
    index[2] = atoi(argv[5]);
    connectedThreshold->SetSeed(index);

    InternalImageType::Pointer image = smoothing->GetOutput();
    reader->Update();

    std::cout << "You selected pixel value " <<  image->GetPixel(index) << std::endl;

    // Running all filters at once
    try {
       writer->Update();
    }
    catch(itk::ExceptionObject & excep) {
       std::cerr << "Exception caught !" << std::endl;
       std::cerr << excep << std::endl;
    }
    return EXIT_SUCCESS;
}
