#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "ref_values.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x;
        double *y;
        double *refValues;
        mwSignedIndex numPoints;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefValuesMex",
                                  "Two inputs required.");
        }
        if (((mxGetM(prhs[0])) != 1) && (mxGetN(prhs[0]) != 1)) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefValuesMex",
                                  "Inputs must be vectors.");
        }
        if (((mxGetM(prhs[1])) != 1) && (mxGetN(prhs[1]) != 1)) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefValuesMex",
                                  "Inputs must be vectors.");
        }

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefValuesMex",
                                  "Requires exactly one output.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        numPoints = mxGetN(prhs[0])*mxGetM(prhs[0]);

        plhs[0] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        refValues = mxGetPr(plhs[0]);

        ap_ref_values(x, y, numPoints, refValues);
}
