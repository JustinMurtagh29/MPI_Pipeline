/* scaleRows.cpp
*
* Written by
* Benedikt Staffler <benedikt.staffler@brain.mpg.de>
*/

#include "mex.h"
#include "matrix.h"
#include <algorithm>

void mexFunction(
	int nlhs, mxArray * plhs[],
	int nrhs, const mxArray * prhs[])
{
    if (nrhs != 2) mexErrMsgTxt("Invalid number of input arguments");
    const mwSize * dims = mxGetDimensions(prhs[0]);

	if (mxIsSingle(prhs[0]))
	{
	    float *m = (float*)mxGetData(prhs[0]);
        float *fact = (float*)mxGetData(prhs[1]);
        
	    for (int i = 0; i < dims[1]; i++)
	    {
            for (int j = 0; j < dims[0]; j++)
            {
                m[i * dims[0] + j] = m[i * dims[0] + j] * fact[i];
            }
	    }
	}
	else
    {
        mexErrMsgTxt("Array must be single or double");
    }
}