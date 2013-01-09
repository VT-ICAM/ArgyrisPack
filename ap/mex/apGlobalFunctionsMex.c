#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "global_functions.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C;
        double *refValues;
        double *values;
        mwSignedIndex quadPoints;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalFunctionsMex",
                                  "Requires two arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalFunctionsMex",
                                  "The C matrix must be the first argument.");
        }
        if (21 != mxGetM(prhs[1])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalFunctionsMex",
                                  "There should be 21 basis functions "
                                  "corresponding to 21 rows.");
        }

        C = mxGetPr(prhs[0]);
        refValues = mxGetPr(prhs[1]);
        quadPoints = mxGetN(prhs[1]);

        /* check output. */
        if (1 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalFunctionsMex",
                                  "Requires exactly one output.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        values = mxGetPr(plhs[0]);

        ap_global_functions(C, refValues, quadPoints, values);
}
