#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "local_functions.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x, *y;
        double *refValues;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalFunctionsMex",
                                  "Requires two arguments.");
        }
        for (i = 0; i < 2; i++) {
                if (mxGetN(prhs[i]) != 1 && mxGetM(prhs[i]) != 1) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apLocalFunctionsMex",
                                          "Inputs must be vectors.");
                }
        }
        if (mxGetM(prhs[0])*mxGetN(prhs[0]) != (mxGetM(prhs[1])*mxGetN(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalFunctionsMex",
                                  "Inputs must have same size.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        quadPoints = mxGetN(prhs[1])*mxGetM(prhs[1]);

        /* check output. */
        if (1 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalFunctionsMex",
                                  "Requires exactly one output.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        refValues = mxGetPr(plhs[0]);

        ap_local_functions(x, y, quadPoints, refValues);
}
