#include "scCalc.hh" 

using namespace std;

void ScCalc::printVector(vector<vector<int> > vec)
{
  unsigned length = vec.size();
  for (unsigned i = 0; i < length; ++i)
  {
      cout << "Point " << i << ": [";
      for (unsigned j = 0; j < 3; ++j)
      {
          cout << vec[i][j];
          if (j != 2) cout << ", ";
      }
      cout << "]" << endl;
  }
}

void ScCalc::saveVector(vector<vector<int> > vec, char* fileName)
{
  unsigned length = vec.size();
  ofstream f(fileName);
  for (unsigned i = 0; i < length; ++i)
  {
      for (unsigned j = 0; j < 3; ++j)
      {
          f << vec[i][j];
          if (j != 2) f << ',';
      }
      f << '\n';
  }
}

vector<vector<int>> ScCalc::loadVector(char* fileName)
{
    FILE *fp;
    vector<vector<int> > vec;

    string sFileName = fileName;
    ifstream fileStream(sFileName);
    if (!fileStream.is_open())
    {
        cout << "Exiting unable to open file " << fileName << endl;
    }

    string line;
    while(getline(fileStream, line, '\n')) {
        stringstream ss(line);
        vector<int> numbers;
        string in_line;
        while(getline (ss, in_line, ',')) 
        {
          int i = stoi(in_line, 0);
          numbers.push_back(i);
        }
        vec.push_back(numbers);
    }

    return vec;
}

void ScCalc::loadSpine1(vector<vector<int>> spine)
{
    unsigned length = spine.size();
    double (*newSpine)[3] = new double[length][3];

    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spacing1[j] * spine[i][j];
        }
    }

    spine1 = newSpine; 
    spine1Length = length;
}

void ScCalc::loadSpine2(vector<vector<int>> spine)
{
    unsigned length = spine.size();
    double (*newSpine)[3] = new double[length][3];

    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spacing2[j] * spine[i][j];
        }
    }
    
    spine2 = newSpine; 
    spine2Length = length;
}

void ScCalc::loadSpine1(double spine[][3], unsigned length)
{
    loadSpineX(spine, length, 0);
}
void ScCalc::loadSpine2(double spine[][3], unsigned length)
{
    loadSpineX(spine, length, 1);
}

void ScCalc::loadSpineX(double spine[][3], unsigned length, unsigned spineNumber)
{
    double (*newSpine)[3] = new double[length][3];
    double* spacing;

    if (!spineNumber) spacing = spacing1;
    else spacing = spacing2;

    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spacing[j] * spine[i][j];
        }
    }

    if (!spineNumber) {spine1 = newSpine; spine1Length = length;}
    else {spine2 = newSpine; spine2Length = length;}
}

void ScCalc::loadTransofrm(double matrix[4][4])
{

    for (unsigned i = 0; i < 4; ++i)
    {
        for (unsigned j = 0; j < 4; ++j)
        {
            transMatrix[i][j] = matrix[i][j];
        }
    }
}

void ScCalc::transformSpine1()
{
    double (*newSpine)[3] = new double[spine1Length][3];

    for (unsigned i = 0; i < spine1Length; ++i)
    {
        for(int j=0; j<3; ++j) 
        {
            newSpine[i][j] = 0;

            for(int k=0; k<3; ++k)
            {
                newSpine[i][j]+=spine1[i][k]*transMatrix[j][k];
            }
            newSpine[i][j]+=transMatrix[j][3];
        }
    }

    spine1 = newSpine;
}

void ScCalc::compareSpines()
{
    double (*ptr)[3] = spine1;

    for (unsigned i = 0; i < spine1Length - 9; ++i)
    {
        cout << "derivative = " << curveDerivative(ptr);
        cout << " second derivative = " << curveSecDerivative(ptr++);
        cout << endl;
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            cout << transMatrix[i][j] << " ";
        }

        cout << endl;
    }
}

void ScCalc::printSpine1()
{
    for (unsigned i = 0; i < spine1Length; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            cout << spine1[i][j] << " ";
        }
        cout << endl;
    }
}

void ScCalc::printAngles()
{
    cout << endl;
    cout << "=======================================" << endl;
    cout << "= Spine Curvature Summary (Degrees)   =" << endl;
    cout << "=======================================" << endl;
    double max1 = 0;
    double max2 = 0;

    unsigned length;
    if (spine1Length > spine2Length)
        length = spine1Length;
    else
        length = spine2Length;

    for (unsigned i = 0; i < length - 3; ++i)
    {
        if (i < spine1Length - 3) 
        {
            cout << "Spine1: " << angles1[i] << ", ";
            if (angles1[i] > max1) max1 = angles1[i];
        } else cout << "               ";
        if (i < spine2Length - 3)
        {
            cout << "Spine2: " << angles2[i];
            if (angles2[i] > max2) max2 = angles2[i];
        }
        cout << endl;
    }
    cout << "=======================================" << endl;
    cout << "Spine1 MAX: " << max1 << " Spine2 MAX: " << max2;
}

double ScCalc::firstDerivative(double a[5])
{
    double firstDer(0.0);
    double t(0.001);

    // no loop here as 5 points is just enough to calculate the derivative
    // of one elem (the middle one). If you do use a loop, it has to start
    // at 2 and end 2 elems before the end (so i - 2 and i + 2 are always
    // valid indexes for the array.)

    size_t i(2);

    firstDer = (   4.0 / 3.0 * (a[i + 1] - a[i - 1]) / (2.0 * t)
                 - 1.0 / 3.0 * (a[i + 2] - a[i - 2]) / (4.0 * t) );
    // Rather than use (double)4 you can just use 4.0. The .0 tells the
    // compiler you mean a double (if you need a float instead, append
    // an f, like 3.0f )

    return firstDer;
}

double ScCalc::curveDerivative(double points[5][3])
{
	double *x = new double[5];
	double *y = new double[5];
	double *z = new double[5];

	for (int i = 0; i < 5; ++i)
	{
		x[i] = points[i][0];
		y[i] = points[i][1];
		z[i] = points[i][2];
	}

	double dX = firstDerivative(x);
	double dY = firstDerivative(y);
	double dZ = firstDerivative(z);

	return sqrt(dX * dX + dY * dY + dZ * dZ);
}

double ScCalc::curveSecDerivative(double points[9][3])
{
	double *dOne = new double[9];

	double (*ptr)[3] = points;

	for (int i = 0; i < 5; ++i)
	{
		dOne[i] = curveDerivative(ptr++);
	}

	return firstDerivative(dOne);
}

// std::vector<T> polyfit( const std::vector<T>& oX, 
//     const std::vector<T>& oY, int nDegree )
// {
//     using namespace boost::numeric::ublas;
 
//     if ( oX.size() != oY.size() )
//         throw std::invalid_argument( "X and Y vector sizes do not match" );
 
//     // more intuative this way
//     nDegree++;
    
//     size_t nCount =  oX.size();
//     matrix<T> oXMatrix( nCount, nDegree );
//     matrix<T> oYMatrix( nCount, 1 );
    
//     // copy y matrix
//     for ( size_t i = 0; i < nCount; i++ )
//     {
//         oYMatrix(i, 0) = oY[i];
//     }
 
//     // create the X matrix
//     for ( size_t nRow = 0; nRow < nCount; nRow++ )
//     {
//         T nVal = 1.0f;
//         for ( int nCol = 0; nCol < nDegree; nCol++ )
//         {
//             oXMatrix(nRow, nCol) = nVal;
//             nVal *= oX[nRow];
//         }
//     }
 
//     // transpose X matrix
//     matrix<T> oXtMatrix( trans(oXMatrix) );
//     // multiply transposed X matrix with X matrix
//     matrix<T> oXtXMatrix( prec_prod(oXtMatrix, oXMatrix) );
//     // multiply transposed X matrix with Y matrix
//     matrix<T> oXtYMatrix( prec_prod(oXtMatrix, oYMatrix) );
 
//     // lu decomposition
//     permutation_matrix<int> pert(oXtXMatrix.size1());
//     const std::size_t singular = lu_factorize(oXtXMatrix, pert);
//     // must be singular
//     BOOST_ASSERT( singular == 0 );
 
//     // backsubstitution
//     lu_substitute(oXtXMatrix, pert, oXtYMatrix);
 
//     // copy the result to coeff
//     return std::vector<T>( oXtYMatrix.data().begin(), oXtYMatrix.data().end() );
// }

void ScCalc::crateSpineFit(int spineNum, unsigned order)
{
    double (*spine)[3];
    double (*fit)[3];
    unsigned spLength;
    double *maXanX;

    if (spineNum == 1) {
        spine = spine1;
        spLength = spine1Length;
        fit1 = new double[spLength][3];
        fit = fit1;
        angles1 = new double[spLength - 3];
        maXanX = angles1;
    } else {
        spine = spine2;
        spLength = spine2Length;
        fit2 = new double[spLength][3];
        fit = fit2;
        angles2 = new double[spLength - 3];
        maXanX = angles2;
    }

    double *dx = new double[spLength];
    double *dy = new double[spLength];
    double *dz = new double[spLength];

    double *xStore = new double[order];
    double *yStore = new double[order];

    for (unsigned i = 0; i < spLength; ++i) 
    {
        dx[i] = spine[i][0];
        dy[i] = spine[i][1];
        dz[i] = spine[i][2];
    }

    polynomialfit(spLength, order, dz, dx, xStore);
    polynomialfit(spLength, order, dz, dy, yStore);

    cout << endl;
    cout << "Spine" << spineNum << " polynomial curve:" << endl;
    cout << endl;

    for (int i = 0; i < order; ++i)
    {
        cout << i << " order term ";
        cout << "X: " << xStore[i];
        cout << " Y: " << yStore[i] << endl;
    }

    for (unsigned i = 0; i < spLength; ++i) 
    {
        double sumX = 0;
        double sumY = 0;

        for (unsigned j = 0; j < order; ++j) 
        {
            sumX += xStore[j] * pow(dz[i], j);
            sumY += yStore[j] * pow(dz[i], j);
        }

        //cout << "Next point: " << sumX << ", " << sumY << ", " << dz[i] << endl;

        fit[i][0] = sumX;
        fit[i][1] = sumY;
        fit[i][2] = dz[i];
    }

    double *anX = new double[spLength];
    double *anY = new double[spLength];
    

    //anglesX(fit, anX, spLength);
    //anglesY(fit, anY, spLength);
    getMax3Dangles(fit, maXanX, spLength);
}


bool ScCalc::polynomialfit(int obs, int degree, 
           double *dx, double *dy, double *store) /* n, p */
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return true; /* we do not "analyse" the result (cov matrix mainly)
          to know if the fit is "good" */
}

void ScCalc::getMax3Dangles(double points[][3], double angles[], int npoints)
{
    for (int i = 1; i < npoints - 2; ++i) 
    {
        double *vec1 = new double[3];
        for (int k = 0; k < 3; ++k) vec1[k] = points[i][k] - points[i-1][k];
        angles[i - 1] = 0;

        for (int j = i + 2; j < npoints; ++j)
        {

            double *vec2 = new double[3];
            for (int k = 0; k < 3; ++k) vec2[k] = points[j][k] - points[j - 1][k];

            double dotProduct = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
            double mag1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
            double mag2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);

            double newAngle = acos(dotProduct / (mag1 * mag2)) * 180/3.1415926358979323;
            if (newAngle > angles[i - 1]) angles[i - 1] = newAngle;
            else break;
        }
    }
}

void ScCalc::maXanglesX(double points[][3], double angles[], int npoints){
    for (int j = 0; j < npoints; j++)
    {
        double previousAngle = 360;

        for(int i = j; i < npoints; i++)
        {
            int last = (j - 1 + npoints) % npoints;
            int next = (i + 1) % npoints;
            double x1 = points[i][0] - points[last][0];
            double z1 = points[i][2] - points[last][2];
            double x2 = points[next][0] - points[i][0];
            double z2 = points[next][2] - points[i][2];
            double theta1 = atan2(z1, x1)*180/3.1415926358979323;
            double theta2 = atan2(z2, x2)*180/3.1415926358979323;
            double thisAngle = (180 + theta1 - theta2 + 360);
            while(thisAngle>360)thisAngle-=360;

            if (previousAngle < thisAngle)
                break;
            else
                previousAngle = thisAngle;
        } 
        angles[j] = previousAngle;
    }
}

void ScCalc::anglesX(double points[][3], double angles[], int npoints){
for(int i = 0; i < npoints; i++){
    int last = (i - 1 + npoints) % npoints;
    int next = (i + 1) % npoints;
    double x1 = points[i][0] - points[last][0];
    double z1 = points[i][2] - points[last][2];
    double x2 = points[next][0] - points[i][0];
    double z2 = points[next][2] - points[i][2];
    double theta1 = atan2(z1, x1)*180/3.1415926358979323;
    double theta2 = atan2(z2, x2)*180/3.1415926358979323;
    angles[i] = (180 + theta1 - theta2 + 360);
    while(angles[i]>360)angles[i]-=360;
} }

void ScCalc::anglesY(double points[][3], double angles[], int npoints){
for(int i = 0; i < npoints; i++){
    int last = (i - 1 + npoints) % npoints;
    int next = (i + 1) % npoints;
    double y1 = points[i][1] - points[last][1];
    double z1 = points[i][2] - points[last][2];
    double y2 = points[next][1] - points[i][1];
    double z2 = points[next][2] - points[i][2];
    double theta1 = atan2(z1, y1)*180/3.1415926358979323;
    double theta2 = atan2(z2, y2)*180/3.1415926358979323;
    angles[i] = (180 + theta1 - theta2 + 360);
    while(angles[i]>360)angles[i]-=360;
} }

// //____________________________________________________________________
// void ScCalc::makeData(Double_t* x, Double_t& d, Double_t& e)
// {
//   // Make data points
//   Double_t upp[5] = { 10, 10, 10, 10,  1 };
//   Double_t low[5] = {  0,  0,  0,  0, .1 };
//   for (int i = 0; i < 2; i++)
//     x[i] = (upp[i] - low[i]) * gRandom->Rndm() + low[i];

//   d = x[0] * TMath::Sqrt(x[1] * x[1]);

//   e = gRandom->Gaus(upp[4],low[4]);
// }

// //____________________________________________________________________
// int ScCalc::CompareResults(TMultiDimFit *fit, bool doFit)
// {
//    //Compare results with reference run


//    // the right coefficients (before fit)
//   double GoodCoeffsNoFit[] = {
//   -4.37056,
//   43.1468,
//   13.432,
//   13.4632,
//   13.3964,
//   13.328,
//   13.3016,
//   13.3519,
//   4.49724,
//   4.63876,
//   4.89036,
//   -3.69982,
//   -3.98618,
//   -3.86195,
//   4.36054,
//   -4.02597,
//   4.57037,
//   4.69845,
//   2.83819,
//   -3.48855,
//   -3.97612
// };

//    // the right coefficients (after fit)
//   double GoodCoeffs[] = {
//      -4.399,
//      43.15,
//      13.41,
//      13.49,
//      13.4,
//      13.23,
//      13.34,
//      13.29,
//      4.523,
//      4.659,
//      4.948,
//      -4.026,
//      -4.045,
//      -3.939,
//      4.421,
//      -4.006,
//      4.626,
//      4.378,
//      3.516,
//      -4.111,
//      -3.823,
// };

// // Good Powers
//   int GoodPower[] = {
//   1,  1, 
//   2,  1, 
//   1,  1,  
//   1,  1, 
//   1,  2,  
//   2,  2,  
//   2,  1,  
//   2,  1, 
//   1,  1,  
//   1,  3,  
//   1,  1,  
//   1,  1,  
//   1,  2,  
//   1,  2,  
//   2,  1,  
//   2,  2,  
//   2,  1,  
//   2,  3,  
//   1,  2,  
//   2,  1,  
//   2,  2
// };

//   Int_t nc = fit->GetNCoefficients();
//   Int_t nv = fit->GetNVariables();
//   const Int_t *powers = fit->GetPowers();
//   const Int_t *pindex = fit->GetPowerIndex();
//   if (nc != 21) return 1;
//   const TVectorD *coeffs = fit->GetCoefficients();
//   int k = 0;
//   for (Int_t i=0;i<nc;i++) {
//      if (doFit) {
//         if (!TMath::AreEqualRel((*coeffs)[i],GoodCoeffs[i],1e-3)) return 2;
//      }
//      else {
//         if (TMath::Abs((*coeffs)[i] - GoodCoeffsNoFit[i]) > 5e-5) return 2;
//      }
//      for (Int_t j=0;j<nv;j++) {
//         if (powers[pindex[i]*nv+j] != GoodPower[k]) return 3;
//         k++;
//      }
//   }

//   // now test the result of the generated function
//   gROOT->ProcessLine(".L MDF.C");

//   Double_t refMDF = (doFit) ? 43.95 : 43.98;
//   // this does not work in CLing since the function is not defined
//   //Double_t x[]    = {5,5,5,5};
//   //Double_t rMDF   = MDF(x);
//   //LM:  need to return the address of the result since it is casted to a long (this should not be in a tutorial !)
//   Long_t iret = gROOT->ProcessLine(" Double_t x[] = {5,5,5,5}; double result=MDF(x); &result;");
//   Double_t rMDF = * ( (Double_t*)iret);
//   //printf("%f\n",rMDF);
//   if (TMath::Abs(rMDF -refMDF) > 1e-2) return 4;
//   return 0;
// }

// //____________________________________________________________________
// Int_t ScCalc::multidimfit(bool doFit)
// {

//   cout << "*************************************************" << endl;
//   cout << "*             Multidimensional Fit              *" << endl;
//   cout << "*                                               *" << endl;
//   cout << "* By Christian Holm <cholm@nbi.dk> 14/10/00     *" << endl;
//   cout << "*************************************************" << endl;
//   cout << endl;

//   // Initialize global TRannom object.
//   gRandom = new TRandom();

//   // Open output file
//   TFile* output = new TFile("mdf.root", "RECREATE");

//   // Global data parameters
//   Int_t nVars       = 2;
//   Int_t nData       = 500;
//   Double_t x[2];

//   // make fit object and set parameters on it.
//   TMultiDimFit* fit = new TMultiDimFit(nVars, TMultiDimFit::kMonomials,"v");

//   Int_t mPowers[]   = { 6 , 6};
//   fit->SetMaxPowers(mPowers);
//   fit->SetMaxFunctions(1000);
//   fit->SetMaxStudy(1000);
//   fit->SetMaxTerms(30);
//   fit->SetPowerLimit(1);
//   fit->SetMinAngle(10);
//   fit->SetMaxAngle(10);
//   fit->SetMinRelativeError(.01);

//   // variables to hold the temporary input data
//   Double_t d;
//   Double_t e;

//   // Print out the start parameters
//   fit->Print("p");

//   printf("======================================\n");

//   // Create training sample
//   Int_t i;
//   for (i = 0; i < nData ; i++) {

//     // Make some data
//     makeData(x,d,e);

//     // Add the row to the fit object
//     fit->AddRow(x,d,e);
//   }

//   // Print out the statistics
//   fit->Print("s");

//   // Book histograms
//   fit->MakeHistograms();

//   // Find the parameterization
//   fit->FindParameterization();

//   // Print coefficents
//   fit->Print("rc");

//   // Get the min and max of variables from the training sample, used
//   // for cuts in test sample.
//   Double_t *xMax = new Double_t[nVars];
//   Double_t *xMin = new Double_t[nVars];
//   for (i = 0; i < nVars; i++) {
//     xMax[i] = (*fit->GetMaxVariables())(i);
//     xMin[i] = (*fit->GetMinVariables())(i);
//   }

//   nData = fit->GetNCoefficients() * 100;
//   Int_t j;

//   // Create test sample
//   for (i = 0; i < nData ; i++) {
//     // Make some data
//     makeData(x,d,e);

//     for (j = 0; j < nVars; j++)
//       if (x[j] < xMin[j] || x[j] > xMax[j])
//     break;

//     // If we get through the loop above, all variables are in range
//     if (j == nVars)
//       // Add the row to the fit object
//       fit->AddTestRow(x,d,e);
//     else
//       i--;
//   }
//   //delete gRandom;

//   // Test the parameterizatio and coefficents using the test sample.
//   if (doFit)
//      fit->Fit("M");

//   // Print result
//   fit->Print("fc v");

//   // Write code to file
//   fit->MakeCode();

//   // Write histograms to disk, and close file
//   output->Write();
//   output->Close();
//   delete output;

//   // Compare results with reference run
//   Int_t compare = CompareResults(fit, doFit);
//   if (!compare) {
//      printf("\nmultidimfit ..............................................  OK\n");
//   } else {
//      printf("\nmultidimfit ..............................................  fails case %d\n",compare);
//   }

//   //What we will need:
//   //fit->Eval(x);

//   // We're done
//   delete fit;
//   return compare;
// }
