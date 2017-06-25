# Spine Curavture Analyais Algorithm (CSC621 Term Project - Team Rain)
A CT scan image analysis utility designed to assess the severity of scoliosis or other spine illness. Workflow consists of (i) GUI seed selection (ii) spine segmentstion (iii) location of centroids (iv) curve fitting (v) curvature calculation (vi) max curavture angle detection. Analysis is performed in 3D, outputs include calculated curvature angles as well as curvature visualization.

## Team
[Maria Khvalchik](https://github.com/mkhvalchik "Maria Khvalchik's Github account") (Lead)  
[Juris Puchin](https://github.com/JurisPuchin "Juris Puchin's Github account") (Lead)  
[Monte Thakkar](https://github.com/monte9 "Monte Thakkar's Github account")  
[Kenneth Shen](https://github.com/kennybeegoode "Kenneth Shen's Github account")  
[Eric C. Black](https://github.com/bitacoustic "Eric C. Black's Github account") 

## Dependencies
[ITK](http://itk.org/)  
[VTK](http://www.vtk.org/)  
[FLTK](http://www.fltk.org/)  
[ItkVtkGlue](https://github.com/InsightSoftwareConsortium/ITKWikiExamples/raw/master/ItkVtkGlue.tar.gz) (direct link)

## Compilation
*Linux* - From within ```/Calculator``` directory, run:
```
mkdir bin
cd bin
cmake ../   // or ccmake if you want to use a GUI for this step
make
(on mac) /FinalProject.app/Contents/MacOS/FinalProject ../output1.mhd ../output2.mhd
(on linux) please add linux command here
```
