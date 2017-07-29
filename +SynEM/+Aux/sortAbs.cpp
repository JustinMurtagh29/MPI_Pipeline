/* sortAbs.cpp
* In placing sorting of the columns of a single/double matrix by absolute
* value.
*
* Written by
* Benedikt Staffler <benedikt.staffler@brain.mpg.de>
*/

#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>

using namespace std;

void mexFunction(
	int nlhs, mxArray * plhs[],
	int nrhs, const mxArray * prhs[])
{
    if (nrhs != 1) mexErrMsgTxt("Invalid number of input arguments");
    const mwSize * dims = mxGetDimensions(prhs[0]);

	if (mxIsSingle(prhs[0]))
	{
	    float *m = (float*)mxGetData(prhs[0]);
	    for (int i = 0; i < dims[1]; i++)
	    {
	        std::sort(m + (i * dims[0]), m + ((i + 1) * dims[0]),
	                [](float a, float b) { return abs(a) < abs(b); });
	    }
	}
	else if (mxIsDouble(prhs[0]))
	{
		double *m = (double*)mxGetData(prhs[0]);
	    for (int i = 0; i < dims[1]; i++)
	    {
	        std::sort(m + (i * dims[0]), m + ((i + 1) * dims[0]),
	                [](double a, double b) { return abs(a) < abs(b); });
	    }
	}
	else
    {
        mexErrMsgTxt("Array must be single or double");
    }
}
