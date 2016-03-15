#ifndef REGIONGROWING
#define REGIONGROWING

#include "RegionGrowingGUI.h"

#include "fltkImageViewer.h"
#include "fltkVTKImageViewer.h"

//includes for vtk rendering
#include "vtkRenderer.h"
#include "vtkRenderWindow.h" 

class  RegionGrowing : public RegionGrowingGUI {
	
	public:
	
		/** image dimensions */
		enum { ImageDimension = RegionGrowingBase::ImageDimension };
		
		/** input pixel type */
		typedef RegionGrowingBase::InputPixelType InputPixelType;
		
		/** output pixel type */
		typedef RegionGrowingBase::OutputPixelType OutputPixelType;
		
		/** internal pixel type */
		typedef RegionGrowingBase::InternalPixelType InternalPixelType;
		
		/** dicom pixel type */
		typedef RegionGrowingBase::DicomPixelType DicomPixelType;
		
		/** input/output image type */
		typedef RegionGrowingBase::InputImageType InputImageType;
		
		/** output image type */
		typedef RegionGrowingBase::OutputImageType OutputImageType;
		
		/** internal image type */
		typedef RegionGrowingBase::InternalImageType InternalImageType;
		
		/** dicom image type */
		typedef RegionGrowingBase::DicomImageType DicomImageType;
		
		/** image viewer types */
		typedef OutputPixelType OverlayPixelType;
		
		typedef fltk::ImageViewer< InputPixelType, 
                                    OverlayPixelType >          InputImageViewerType;

		typedef fltk::ImageViewer< InternalPixelType, 
                                    OverlayPixelType >          InternalImageViewerType;

		typedef fltk::ImageViewer< InputPixelType, 
                                    OverlayPixelType >          OutputImageViewerType;
									
		typedef fltk::ImageViewer< DicomPixelType, 
                                    OverlayPixelType >          DicomImageViewerType;

		typedef fltk::VTKImageViewer< OutputPixelType >               VTKImageViewerType;
		
	public:
	
		float float_xStartValueInput;
		float float_yStartValueInput;
		float float_zStartValueInput;
		float float_xEndValueInput;
	    float float_yEndValueInput;
		float float_zEndValueInput;
		
		RegionGrowing();
		virtual ~RegionGrowing();
		
		virtual void LoadInputImage();
		virtual void LoadInputImageSeries();
		virtual void WriteOutputImage();
		virtual void ShowConsole();
		virtual void ShowStatus( const char * text);
		virtual void SaveConfConSeries();
		virtual void SaveCustomSeries();
		
		/** preprocessing show functions */
		virtual void ShowGradientAnisotropicDiffusionImage();
		virtual void ShowCurvatureAnisotropicDiffusionImage();
		
		/** filter / image show functions */
		virtual void ShowInputImage();
		virtual void ShowConnectedThresholdImage();
		virtual void ShowConfidenceConnectedImage();
		virtual void ShowCustomRegionGrowingImage();
		virtual void ShowHomogeneousImage();
		virtual void ShowHomogeneousImageWithVTK();
		
		/** volume function */
		virtual void ShowVolume();
		virtual void ShowCustomVolume();
		
		/** writer definitions */
		virtual void WriteConnectedThresholdImage();
		virtual void WriteConfidenceConnectedImage();
		
		virtual void Quit();
		
		virtual void SelectSeedPoint( float x, float y, float z);
		
		static void ClickSelectCallback( float x, float y, float z, float value, void * args );
		
	private:
	
		/** filtered image viewers */
		InternalImageViewerType     m_HomogeneousImageViewer;
		
		/** input image viewer */
		InputImageViewerType 		m_InputImageViewer;
		DicomImageViewerType 		m_DicomImageViewer;
		
		/** prerprocessing image viewers */
		InternalImageViewerType		m_GradientAnisotropicDiffusionImageViewer;
		InternalImageViewerType		m_CurvatureAnisotropicDiffusionImageViewer;
		
		/** region growing filters */
		InternalImageViewerType       m_ConnectedThresholdImageViewer;
		InternalImageViewerType       m_ConfidenceConnectedImageViewer;
		InternalImageViewerType		  m_CustomRegionGrowingImageViewer;
		
		/** segmented image viewer */
		VTKImageViewerType::Pointer   m_VTKSegmentedImageViewer;
		
		
		
};

#endif