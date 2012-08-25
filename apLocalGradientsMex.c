#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "local_gradients.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x, *y;
        double *refdx, *refdy;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalGradientsMex",
                                  "Requires two arguments.");
        }
        for (i = 0; i < 2; i++) {
                if (mxGetN(prhs[i]) != 1 && mxGetM(prhs[i]) != 1) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apLocalGradientsMex",
                                          "Inputs must be vectors.");
                }
        }
        if (mxGetM(prhs[0])*mxGetN(prhs[0]) != (mxGetM(prhs[1])*mxGetN(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalGradientsMex",
                                  "Inputs must have same size.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        quadPoints = mxGetN(prhs[1])*mxGetM(prhs[1]);

        /* check output. */
        if (2 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apLocalGradientsMex",
                                  "Requires exactly two outputs.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        refdx = mxGetPr(plhs[0]);
        refdy = mxGetPr(plhs[1]);

        ap_local_gradients(x, y, quadPoints, refdx, refdy);
}
