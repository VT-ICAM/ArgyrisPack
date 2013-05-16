#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "physical_values.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C;
        double *refValues;
        double *values;
        mwSignedIndex numPoints;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalValuesMex",
                                  "Requires two arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalValuesMex",
                                  "The C matrix must be the first argument.");
        }
        if (21 != mxGetM(prhs[1])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalValuesMex",
                                  "There should be 21 basis values "
                                  "corresponding to 21 rows.");
        }

        C = mxGetPr(prhs[0]);
        refValues = mxGetPr(prhs[1]);
        numPoints = mxGetN(prhs[1]);

        /* check output. */
        if (1 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalValuesMex",
                                  "Requires exactly one output.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        values = mxGetPr(plhs[0]);

        ap_physical_values(C, refValues, numPoints, values);
}
