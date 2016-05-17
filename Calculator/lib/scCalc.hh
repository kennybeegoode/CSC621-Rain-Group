// in myclass.h
#ifndef __SCCALC_HH_INCLUDED__   // if x.h hasn't been included yet...
#define __SCCALC_HH_INCLUDED__

#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>

// #include "Riostream.h"
// #include "TROOT.h"
// #include "TApplication.h"
// #include "TCanvas.h"
// #include "TH1.h"
// #include "TSystem.h"
// #include "TBrowser.h"
// #include "TFile.h"
// #include "TRandom.h"
// #include "TMultiDimFit.h"
// #include "TVectorD.h"
// #include "TMath.h"

#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <math.h>

//Spine Curvature Calculator
class ScCalc
{
private:
    double transMatrix[4][4];
    void loadSpineX(double spine[][3], unsigned length, unsigned spineNumber);
public:
    std::vector<std::vector<int> > spine1vec;
    double (*spine1)[3];
    double (*fit1)[3];
    unsigned spine1Length;
    double* spacing1;
    double* angles1;

    std::vector<std::vector<int> > spine2vec;
    double (*spine2)[3];
    double (*fit2)[3];
    unsigned spine2Length;
    double* spacing2;
    double* angles2;

    // TRandom* gRandom;
    
    void printVector(std::vector<std::vector<int> > vec);
    void saveVector(std::vector<std::vector<int> > vec, char* fileName);
    std::vector<std::vector<int> > loadVector(char* fileName);
    void loadSpine1(double spine[][3], unsigned length);
    void loadSpine2(double spine[][3], unsigned length);
    void loadSpine1(std::vector<std::vector<int>> spine);
    void loadSpine2(std::vector<std::vector<int>> spine);
    void loadTransofrm(double matrix[4][4]);
    void transformSpine1();
    void compareSpines();
    void printSpine1();
    void printTransform();
    void printAngles();
	//TODO: add mean square distance function
    //TODO: curve fit function

	//TODO: add circular calculation funcion
	//TODO: add dot product over determinant product funcion
	//TODO: add smoothing function
	static double firstDerivative(double currentData[5]);
	static double curveDerivative(double points[5][3]);
	static double curveSecDerivative(double points[9][3]);

    void crateSpineFit(int spineNum, unsigned order = 7);
    bool polynomialfit(int obs, int degree, double *dx, 
        double *dy, double *store);

    void anglesX(double points[][3], double angles[], int npoints);
    void anglesY(double points[][3], double angles[], int npoints);
    void maXanglesX(double points[][3], double angles[], int npoints);
    void getMax3Dangles(double points[][3], double angles[], int npoints);
    // Int_t multidimfit(bool doFit = true);
    // int CompareResults(TMultiDimFit *fit, bool doFit);
    // void makeData(Double_t* x, Double_t& d, Double_t& e);
};

#endif
