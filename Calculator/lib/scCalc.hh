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
public:
 static double firstDerivative(double currentData[5]);
 static double curveDerivative(double points[5][3]);
 static double curveSecDerivative(double points[9][3]);
};

#endif