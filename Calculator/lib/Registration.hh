#ifndef _REGISTRATION_HH__
#define _REGISTRATION_HH__

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// image file I/O
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"

// 3d Affine Registration
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

//3D Rigid Transform Registration
#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkVersorRigid3DTransform.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkExtractImageFilter.h"

#include <fstream>
#include "stdint.h"
#include <string>
#include <sstream>  // std::stringstream
#include <stdlib.h> // std::atof
#include <iomanip>  // std::setw

//
//  The following piece of code implements an observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

//Spine Curvature Calculator
class Registration
{
public:
	void rigidAlign(std::string fixedImageInput, std::string movingImageInput, double transformParameters[][4], int maxNumberOfIterations);
	void affineAlign(std::string fixedImageInput, std::string movingImageInput,
			double transformParameters[][4]);
	void affineAlign(std::string fixedImageInput,	std::string movingImageInput,
			double transformParameters[][4], int maxNumberOfIterations);
};

#endif
