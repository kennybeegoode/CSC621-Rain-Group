#include "RegionGrowing.h"
#include <FL/fl_file_chooser.H>

/** constructor */
RegionGrowing::RegionGrowing() {

	m_InputImageViewer.SetLabel("Input Image");
	m_DicomImageViewer.SetLabel("DICOM Image");
	
	/** Image Viewer parameters */
	m_ConnectedThresholdImageViewer.SetLabel("Connected Threshold Image");
	m_ConfidenceConnectedImageViewer.SetLabel("Confidence Connected Image");
	m_CustomRegionGrowingImageViewer.SetLabel("Custom Region Growing Image");
	
	m_CurvatureAnisotropicDiffusionImageViewer.SetLabel("Curvature Anisotropic Diffusion Image");
	m_CurvatureAnisotropicDiffusionImageViewer.ClickSelectCallBack( ClickSelectCallback, (void *)this);
	m_CurvatureAnisotropicDiffusionImageViewer.SetImage( m_CurvatureAnisotropicDiffusionImageFilter->GetOutput() ); 
	
	m_GradientAnisotropicDiffusionImageViewer.SetLabel("Gradient Anisotropic Diffusion Image");
	m_GradientAnisotropicDiffusionImageViewer.ClickSelectCallBack( ClickSelectCallback, (void *)this);
	m_GradientAnisotropicDiffusionImageViewer.SetImage( m_GradientAnisotropicDiffusionImageFilter->GetOutput() );
	
	m_HomogeneousImageViewer.SetLabel("Homogeneous Image");
	m_HomogeneousImageViewer.ClickSelectCallBack( ClickSelectCallback, (void *)this);
	m_HomogeneousImageViewer.SetImage( m_NullImageFilter->GetOutput() );
	
	/** connect GUI values to preprocessing filter parameters */
	m_CurvatureAnisotropicDiffusionImageFilter->SetNumberOfIterations(
         static_cast<unsigned int>(curvatureAnisotropicDiffusionIterationsValueInput->value()) );

	m_CurvatureAnisotropicDiffusionImageFilter->SetTimeStep(
		curvatureAnisotropicDiffusionTimeStepValueInput->value() );

	m_CurvatureAnisotropicDiffusionImageFilter->SetConductanceParameter(
		curvatureAnisotropicDiffusionConductanceValueInput->value() );
		
	m_GradientAnisotropicDiffusionImageFilter->SetNumberOfIterations(
         static_cast<unsigned int>(gradientAnisotropicDiffusionIterationsValueInput->value()) );

	m_GradientAnisotropicDiffusionImageFilter->SetTimeStep(
                                   gradientAnisotropicDiffusionTimeStepValueInput->value() );

	m_GradientAnisotropicDiffusionImageFilter->SetConductanceParameter(
                                   gradientAnisotropicDiffusionConductanceValueInput->value() );
	
	// init itk filter
	//m_ConnectedThresholdImageFilter->SetLower( 
    //  static_cast<unsigned long>( lowerThresholdCounter->value() ) );

	//m_ConnectedThresholdImageFilter->SetUpper( 
    //  static_cast<unsigned long>( upperThresholdCounter->value() ) );
	  
	m_ConfidenceConnectedImageFilter->SetMultiplier( multiplierValueInput->value() );

	m_ConfidenceConnectedImageFilter->SetNumberOfIterations( 
      static_cast<InputPixelType>( iterationsConfidenceValueInput->value() ) );
	  
	/*********** HARD VALUES FOR OUR CUSTOM REGION GROWING FILTER **************/
	//m_CustomRegionGrowingImageFilter->SetMultiplier( static_cast<InputPixelType>( 2.5 ) );
	//m_CustomRegionGrowingImageFilter->SetNumberOfIterations( static_cast<InputPixelType>( 2 ) );
	/***************************************************************************/
	
	m_VTKSegmentedImageViewer = VTKImageViewerType::New();
	//m_VTKSegmentedImageViewer->SetImage( m_ConfidenceConnectedImageFilter->GetOutput() );
	m_VTKSegmentedImageViewer->SetImage( m_CastImageFilter2->GetOutput() );
	
	// GUI Observers
	//inputImageButton->Observe( m_ImageReader.GetPointer() );
	inputImageButton->Observe( m_DicomReader.GetPointer() );
	//thresholdConnectedImageButton->Observe( m_ConnectedThresholdImageFilter.GetPointer() );
	confidenceConnectedImageButton->Observe( m_ConfidenceConnectedImageFilter.GetPointer() );
	customRegionGrowingImageButton->Observe( m_CustomRegionGrowingImageFilter.GetPointer() );
	gradientAnisotropicDiffusionImageButton->Observe( m_GradientAnisotropicDiffusionImageFilter.GetPointer() );
	curvatureAnisotropicDiffusionImageButton->Observe( m_CurvatureAnisotropicDiffusionImageFilter.GetPointer() );
	homogeneousImageVTKButton->Observe( m_NullImageFilter.GetPointer() );
	
	//progressSlider->Observe( m_NullImageFilter.GetPointer() );
	//progressSlider->Observe( m_ConnectedThresholdImageFilter.GetPointer() );
	//progressSlider->Observe( m_ConfidenceConnectedImageFilter.GetPointer() );
	
}

/** Destructor */
RegionGrowing::~RegionGrowing() {}

/** show main console */
void RegionGrowing::ShowConsole( void )
{
	consoleWindow->show();
}

/** quit */
void RegionGrowing::Quit(void) 
{
	m_InputImageViewer.Hide();
	m_DicomImageViewer.Hide();
	m_ConnectedThresholdImageViewer.Hide();
	m_ConfidenceConnectedImageViewer.Hide();
	m_CustomRegionGrowingImageViewer.Hide();
	m_HomogeneousImageViewer.Hide();
	m_VTKSegmentedImageViewer->Hide();
	consoleWindow->hide();
}

/** save output image */
void RegionGrowing::WriteOutputImage( void )
{	
	const char * filename = fl_file_chooser("Output Image filename","*.*","");
	if( !filename )
	{
		return;
	}

	this->ShowStatus("Writing output image file...");
  
	try 
	{
		RegionGrowingBase::WriteOutputImage( filename );
	}
	catch( ... ) 
	{
		this->ShowStatus("Problems writing image");
		return;
	}

	this->ShowStatus("Output Image saved");
}

void RegionGrowing::SaveConfConSeries( void )
{
	const  char * dirname = fl_dir_chooser("Choose DICOM image save folder",0,0);
	
	try
	{
		RegionGrowingBase::SaveConfConSeries( dirname );
	}
	catch( ... )
	{
		this->ShowStatus("Problems reading file format");
		controlsGroup->deactivate();
		return;
	}
}

void RegionGrowing::SaveCustomSeries( void )
{
	const  char * dirname = fl_dir_chooser("Choose DICOM image save folder",0,0);
	
	try
	{
		RegionGrowingBase::SaveCustomSeries( dirname );
	}
	catch( ... )
	{
		this->ShowStatus("Problems reading file format");
		controlsGroup->deactivate();
		return;
	}
}

/** write confidence connected image */
void RegionGrowing::WriteConfidenceConnectedImage( void )
{
	m_ImageWriter->SetInput( m_ConfidenceConnectedImageFilter->GetOutput() );
	this->WriteOutputImage();
}

/** write connected threshold image */
void RegionGrowing::WriteConnectedThresholdImage( void )
{
	m_ImageWriter->SetInput( m_ConnectedThresholdImageFilter->GetOutput() );
	this->WriteOutputImage();
}

/** load input image */
void RegionGrowing::LoadInputImage( void )
{
	const char * filename = fl_file_chooser("Input Image filename","*.*","");
	if (!filename)
	{
		return;
	}
	
	this->ShowStatus("Loading input image file...");
	
	try
	{
		RegionGrowingBase::LoadInputImage( filename );
	}
	catch( ... )
	{
		this->ShowStatus("Problems reading file format");
		controlsGroup->deactivate();
		return;
	}
	
	this->ShowStatus("Input Image Loaded");
	
	controlsGroup->activate();
	
	InputImageType::RegionType region = m_ImageReader->GetOutput()->GetBufferedRegion();

	InputImageType::IndexType start = region.GetIndex();
	InputImageType::SizeType  size  = region.GetSize();
 
	/** x/y/zStartValueInput meant to come from GUI for region selection */
	/**
	xStartValueInput->value( start[0] );
	yStartValueInput->value( start[1] );
	zStartValueInput->value( start[2] );

	xEndValueInput->value( start[0]+size[0] );
	yEndValueInput->value( start[1]+size[1] );
	zEndValueInput->value( start[2]+size[2] );
	*/
	
	/** default region values, remove if region selected implemented */
	float_xStartValueInput = start[0];
    float_yStartValueInput = start[1];
	float_zStartValueInput = start[2];
	
	float_xEndValueInput =  start[0]+size[0];
	float_yEndValueInput =  start[0]+size[0];
	float_zEndValueInput =  start[0]+size[0];
	
}

/** load input image */
void RegionGrowing::LoadInputImageSeries( void )
{
	const  char * dirname = fl_dir_chooser("Choose DICOM image folder",0,0);
	
	try
	{
		RegionGrowingBase::LoadInputImageSeries( dirname );
	}
	catch( ... )
	{
		this->ShowStatus("Problems reading file format");
		controlsGroup->deactivate();
		return;
	}
	
	this->ShowStatus("Input Image Loaded");
	
	m_DicomReader->Update();
	
	controlsGroup->activate();
	
	InputImageType::RegionType region = m_DicomReader->GetOutput()->GetBufferedRegion();

	InputImageType::IndexType start = region.GetIndex();
	InputImageType::SizeType  size  = region.GetSize();
	
	float_xStartValueInput = start[0];
    float_yStartValueInput = start[1];
	float_zStartValueInput = start[2];
	
	float_xEndValueInput =  start[0]+size[0];
	float_yEndValueInput =  start[0]+size[0];
	float_zEndValueInput =  start[0]+size[0];
	
	this->UpdateExtract();
	
}

/** show status */
void RegionGrowing::ShowStatus( const char * message )
{
	//statusTextOutput->value( message );
	//Fl::check();
}

/** show input image */
void RegionGrowing::ShowInputImage( void )
{
	if( !m_InputImageIsLoaded )
    {
		return;
    }
	
	m_DicomToInternalImageTypeFilter->SetInput( m_DicomReader->GetOutput() );
	m_DicomToInternalImageTypeFilter->Update();
	m_HomogeneousImageViewer.SetImage( m_DicomToInternalImageTypeFilter->GetOutput() );		
	m_HomogeneousImageViewer.Show();
	
}

/***********************************************
	Window Call for Custom Region Growing Filter
***********************************************/

void RegionGrowing::ShowCustomRegionGrowingImage( void )
{
	m_CustomRegionGrowingImageFilter->Update();
	if (m_InputImageIsDICOM)
	{
		m_DicomToInternalImageTypeFilter->SetInput( m_DicomReader->GetOutput() );
		m_DicomToInternalImageTypeFilter->Update();
		m_CustomRegionGrowingImageViewer.SetImage( m_DicomToInternalImageTypeFilter->GetOutput() );
	} else {
		m_CustomRegionGrowingImageViewer.SetImage( m_ImageReader->GetOutput() );  
	} 
	m_CustomRegionGrowingImageViewer.SetOverlay( m_CustomRegionGrowingImageFilter->GetOutput() );
	m_CustomRegionGrowingImageViewer.Show();
}

/***********************************************/

/** show threshold connected image */
void RegionGrowing::ShowConnectedThresholdImage( void )
{
	m_ConnectedThresholdImageFilter->Update();
	if (m_InputImageIsDICOM)
	{
		m_DicomToInternalImageTypeFilter->SetInput( m_DicomReader->GetOutput() );
		m_DicomToInternalImageTypeFilter->Update();
		m_ConnectedThresholdImageViewer.SetImage( m_DicomToInternalImageTypeFilter->GetOutput() );
	} else {
		m_ConnectedThresholdImageViewer.SetImage( m_ImageReader->GetOutput() );  
	}
	
	m_ConnectedThresholdImageViewer.SetOverlay( m_ConnectedThresholdImageFilter->GetOutput() );
	m_ConnectedThresholdImageViewer.Show();
}

/** show confidence connected image */
void RegionGrowing::ShowConfidenceConnectedImage( void )
{
	m_NullImageFilter->Update();
	m_ConfidenceConnectedImageFilter->SetInput(m_NullImageFilter->GetOutput());
	m_ConfidenceConnectedImageFilter->Update();
	m_ConfidenceConnectedImageViewer.SetImage( m_NullImageFilter->GetOutput() );
	m_ConfidenceConnectedImageViewer.SetOverlay( m_ConfidenceConnectedImageFilter->GetOutput() );
	m_ConfidenceConnectedImageViewer.Show();
}

void RegionGrowing::ShowVolume( void )
{
	std::cout << "ShowVolume: Calculating volume...";
	unsigned int volume = 0;
	unsigned int tot_volume = 0;
	std::cout << "done!\n";
	std::cout << "ShowVolume: getting output...";
	OutputImageType::Pointer img = m_ConfidenceConnectedImageFilter->GetOutput();
	std::cout << "done!\n";
	std::cout << "ShowVolume: instantiating iterator...";
	ImageIterator  it ( img, img->GetRequestedRegion() );
	std::cout << "done!\n";
	std::cout << "ShowVolume: go to iterator beginning...";
	it.GoToBegin();
	std::cout << "done!\n";
	std::cout << "ShowVolume: begin iterating...";
	while( !it.IsAtEnd() )
	{
		//std::cout << "ShowVolume: While: 1\n";
		OutputPixelType val = it.Get();
		//std::cout << "ShowVolume: While: 2\n";
		if (val>0)
		{
			//std::cout << "ShowVolume: While: 3\n";
			volume++;
			//std::cout << "ShowVolume: While: 4\n";
		}
		//std::cout << "ShowVolume: While: 5\n";
		++tot_volume;
		++it;
    }
	std::cout << "done!\n";
	
	std::cout << "ShowVolume: volume = " << volume << "\n";
	
	std::cout << "ShowVolume: setting volume output in GUI...";
	volumeOutput->value( volume );
	totalVolumeOutput->value( tot_volume );
	std::cout << "done!\n";

}

void RegionGrowing::ShowCustomVolume( void )
{
	std::cout << "ShowVolume: Calculating volume...";
	unsigned int volume = 0;
	unsigned int tot_volume = 0;
	std::cout << "done!\n";
	std::cout << "ShowVolume: getting output...";
	OutputImageType::Pointer img = m_CustomRegionGrowingImageFilter->GetOutput();
	std::cout << "done!\n";
	std::cout << "ShowVolume: instantiating iterator...";
	ImageIterator  it ( img, img->GetRequestedRegion() );
	std::cout << "done!\n";
	std::cout << "ShowVolume: go to iterator beginning...";
	it.GoToBegin();
	std::cout << "done!\n";
	std::cout << "ShowVolume: begin iterating...";
	while( !it.IsAtEnd() )
	{
		//std::cout << "ShowVolume: While: 1\n";
		OutputPixelType val = it.Get();
		//std::cout << "ShowVolume: While: 2\n";
		if (val>0)
		{
			//std::cout << "ShowVolume: While: 3\n";
			volume++;
			//std::cout << "ShowVolume: While: 4\n";
		}
		//std::cout << "ShowVolume: While: 5\n";
		++tot_volume;
		++it;
    }
	std::cout << "done!\n";
	
	std::cout << "ShowVolume: volume = " << volume << "\n";
	
	std::cout << "ShowVolume: setting volume output in GUI...";
	volumeCustomVolumeOutput->value( volume );
	totalCustomVolumeOutput->value( tot_volume );
	std::cout << "done!\n";

}

/** show homogeneous image */
void RegionGrowing::ShowHomogeneousImage( void )
{
	m_NullImageFilter->Update();
	m_HomogeneousImageViewer.SetImage( m_NullImageFilter->GetOutput() );  
	m_HomogeneousImageViewer.Show();
}

/** show homogeneous image w/ vtk */
void RegionGrowing::ShowHomogeneousImageWithVTK( void )
{
	//m_VTKSegmentedImageViewer->SetColorLevel(127.5);
	//m_VTKSegmentedImageViewer->SetColorWindow(255); 
	m_VTKSegmentedImageViewer->Show();
}

/** show gradient anisotropic diffusion image */
void RegionGrowing::ShowGradientAnisotropicDiffusionImage( void )
{
  m_GradientAnisotropicDiffusionImageFilter->Update();
  m_GradientAnisotropicDiffusionImageViewer.SetImage( m_GradientAnisotropicDiffusionImageFilter->GetOutput() );  
  m_GradientAnisotropicDiffusionImageViewer.Show();

}

/** Show curvature anisotropic diffusion image */
void RegionGrowing::ShowCurvatureAnisotropicDiffusionImage( void )
{
  m_CurvatureAnisotropicDiffusionImageFilter->Update();
  m_CurvatureAnisotropicDiffusionImageViewer.SetImage( m_CurvatureAnisotropicDiffusionImageFilter->GetOutput() );  
  m_CurvatureAnisotropicDiffusionImageViewer.Show();

}

/** INSERT OTHER FILTER METHODS HERE */

/** click select seed point callback */
void RegionGrowing::ClickSelectCallback(float x, float y, float z, float itkNotUsed(value), void * args)
{
	RegionGrowing * self = 
		static_cast<RegionGrowing *>( args );

	self->SelectSeedPoint( x, y, z );
}

/** select seed point callback */
void RegionGrowing::SelectSeedPoint(float x, float y, float z)
{
	xSeedPointValueOutput->value( x );
	ySeedPointValueOutput->value( y );
	zSeedPointValueOutput->value( z );

	typedef ConnectedThresholdImageFilterType::IndexType IndexType;
	IndexType seed;
	seed[0] = static_cast<IndexType::IndexValueType>( x );
	seed[1] = static_cast<IndexType::IndexValueType>( y );
	seed[2] = static_cast<IndexType::IndexValueType>( z );

	m_ConnectedThresholdImageFilter->SetSeed( seed );
	m_ConfidenceConnectedImageFilter->SetSeed( seed );
	m_CustomRegionGrowingImageFilter->SetSeed( seed );
}

/**
RegionGrowing::ShowGPURayTracedVolume(InternalImageType * sourceImage, const char * title)
{
	try
	{
		//Prepare VTK data
		m_InternalToVTKImageTypeFilter = InternalToVTKImageTypeFilterType::New();
		m_InternalToVTKImageTypeFilter->SetOutputMinimum(   0 );
		m_InternalToVTKImageTypeFilter->SetOutputMaximum( 4128 );
		m_InternalToVTKImageTypeFilter->SetInput(sourceImage);
		
		// create connector
		ConnectorType::Pointer m_Connector = ConnectorType::New();
		m_Connector->SetInput(m_InternalToVTKImageTypeFilter->GetOutput());
		m_Connector->GetImporter()->SetDataScalarTypeToUnsignedShort();
		m_Connector->Update();	
		
		
		vtkRenderer * arenderer = vtkRenderer::New(); 
		vtkRenderWindow * renwin = vtkRenderWindow::New(); 
		vtkVolume * vol=vtkVolume::New(); 

		
		// Create and initialize the mapper
		vtkKWEVolumeMapper *mapper = vtkKWEVolumeMapper::New();
		mapper->SetInput(m_Connector->GetOutput());
		mapper->SetBlendModeToComposite();
		
		// Create and initialize the volume
		vtkVolume *volume = vtkVolume::New();
		volume->SetMapper(mapper);
		// Create transfer function - mapping scalar values [0-...] to COLOR [0-1, 0-1, 0-1]
		vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
		//colorTransferFunction->AddRGBPoint(   0.0, 0.0,0.0,0.0);
		//colorTransferFunction->AddRGBPoint( 500.0, 0.9,0.5,0.3);
		//colorTransferFunction->AddRGBPoint(1100.0, 0.8,0.8,0.6);
		//colorTransferFunction->AddRGBPoint(1200.0, 0.6,0.6,0.6);
		colorTransferFunction->AddRGBPoint(1024, 1.0,0.0,0.0);
		colorTransferFunction->AddRGBPoint(2048, 0.0,1.0,0.0);
		colorTransferFunction->AddRGBPoint(3096, 0.0,0.0,1.0);

		
		vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();
		opacityTransferFunction->AddPoint(    0, 0.0);
		opacityTransferFunction->AddPoint(    1, 0.0);
		//opacityTransferFunction->AddPoint(  980, 0.1);
		//opacityTransferFunction->AddPoint(  1055, 0.2);
		opacityTransferFunction->AddPoint(  1024, 0.05);
		opacityTransferFunction->AddPoint(  2048, 0.2);
		opacityTransferFunction->AddPoint( 3096, 1.0);
		// The property describes how the data will look
		vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
		volumeProperty->SetColor(colorTransferFunction);
		volumeProperty->SetScalarOpacity(opacityTransferFunction);
		//volumeProperty->ShadeOn(); // request 3d shading (german:> beleuchtungsberechnung anfordern) //has to be implemented in vtkCastRay_OwnWork
			
		volume->SetProperty(volumeProperty);

		// basic camera
		vtkCamera *cam=vtkCamera::New(); 
		cam->SetViewUp(0,0,-1); 
		cam->SetPosition(0,1,0); 
		cam->SetFocalPoint(0,0,0); 
		cam->ComputeViewPlaneNormal();

		arenderer->AddActor(volume); 
		arenderer->SetActiveCamera(cam); 
		arenderer->ResetCamera(); 
		arenderer->SetBackground(1,1,1);

		renwin->AddRenderer(arenderer);
		renwin->SetWindowName(title);
		vtkRenderWindowInteractor *iren=vtkRenderWindowInteractor::New(); 
		iren->SetRenderWindow(renwin); 
		renwin->Render(); 
		iren->Initialize(); 
		iren->Start();
		iren->Delete();
		renwin->Delete();
		cam->Delete();
		volumeProperty->Delete();
		volume->Delete();
	}
	catch(itk::ExceptionObject &ex)
	{
		std::cout << ex.GetDescription() << std::endl;
	}
	return true;
}
*/

/** main */
int main()
{
	try 
	{
		RegionGrowing * console = new RegionGrowing();
		console->ShowConsole();
		Fl::run();
		delete console;
    }
	catch( itk::ExceptionObject & e )
    {
		std::cerr << "ITK exception caught in main" << std::endl;
		std::cerr << e << std::endl;
    }
	catch( std::exception & e )
    {
		std::cerr << "STD exception caught in main" << std::endl;
		std::cerr << e.what() << std::endl;
    }
	catch( ... )
    {
		std::cerr << "unknown exception caught in main" << std::endl;
    }


	return 0;
}