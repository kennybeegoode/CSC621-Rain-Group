#include "scCalc.hh" 

using namespace std;

void ScCalc::loadSpineX(double spine[][3], unsigned length, unsigned spineNumber)
{
    double (*newSpine)[3] = new double[length][3];

    for (unsigned i = 0; i < length; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            newSpine[i][j] = spine[i][j];
        }
    }

    if (!spineNumber) {spine1 = newSpine; spine1Length = length;}
    else {spine2 = newSpine; spine2Length = length;}
}

void ScCalc::loadTransofrm(double matrix[3][3])
{
    //transMatrix = {{0}};

    for (unsigned i = 0; i < 3; ++i)
    {
        for (unsigned j = 0; j < 3; ++j)
        {
            transMatrix[i][j] = matrix[i][j];
        }
    }
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

void ScCalc::loadSpine1(double spine[][3], unsigned length)
{
    loadSpineX(spine, length, 0);
}
void ScCalc::loadSpine2(double spine[][3], unsigned length)
{
    loadSpineX(spine, length, 1);
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