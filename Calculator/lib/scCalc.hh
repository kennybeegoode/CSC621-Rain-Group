// in myclass.h
#ifndef __SCCALC_HH_INCLUDED__   // if x.h hasn't been included yet...
#define __SCCALC_HH_INCLUDED__

#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <iostream>

//Spine Curvature Calculator
class ScCalc
{
private:
    double (*spine1)[3];
    unsigned spine1Length;

    double (*spine2)[3];
    unsigned spine2Length;

    double transMatrix[3][3];
    void loadSpineX(double spine[][3], unsigned length, unsigned spineNumber);
public:
    void loadSpine1(double spine[][3], unsigned length);
    void loadSpine2(double spine[][3], unsigned length);
    void loadTransofrm(double matrix[3][3]);
    // void runTransofrm();
    void compareSpines();
	//TODO: add transpose/rotate function
	//TODO: add mean square distance function
	//TODO: add circular calculation funcion
	//TODO: add dot product over determinant product funcion
	//TODO: add smoothing function
	static double firstDerivative(double currentData[5]);
	static double curveDerivative(double points[5][3]);
	static double curveSecDerivative(double points[9][3]);
};

#endif