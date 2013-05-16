#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "ref_gradients.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x, *y;
        double *refdx, *refdy;
        mwSignedIndex numPoints;
        int i;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefGradientsMex",
                                  "Requires two arguments.");
        }
        for (i = 0; i < 2; i++) {
                if (mxGetN(prhs[i]) != 1 && mxGetM(prhs[i]) != 1) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apRefGradientsMex",
                                          "Inputs must be vectors.");
                }
        }
        if (mxGetM(prhs[0])*mxGetN(prhs[0]) != (mxGetM(prhs[1])*mxGetN(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefGradientsMex",
                                  "Inputs must have same size.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        numPoints = mxGetN(prhs[1])*mxGetM(prhs[1]);

        /* check output. */
        if (2 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefGradientsMex",
                                  "Requires exactly two outputs.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        refdx = mxGetPr(plhs[0]);
        refdy = mxGetPr(plhs[1]);

        ap_ref_gradients(x, y, numPoints, refdx, refdy);
}
