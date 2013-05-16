#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "physical_gradients.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *B, *C;
        double *ref_dx, *ref_dy;
        double *dx, *dy;
        mwSignedIndex quadPoints;

        /* check input. */
        if (nrhs != 4) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "Four inputs required.");
        }
        if (21 != mxGetN(prhs[0]) || 21 != mxGetM(prhs[0])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "The B matrix must be the second argument.");
        }
        if (mxGetN(prhs[2]) != mxGetN(prhs[3])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "Mismatch in number of quadrature points.");
        }
        if (21 != mxGetM(prhs[2]) || 21 != mxGetM(prhs[3])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "There should be 21 basis functions "
                                  "corresponding to 21 rows.");
        }

        C = mxGetPr(prhs[0]);
        B = mxGetPr(prhs[1]);
        ref_dx = mxGetPr(prhs[2]);
        ref_dy = mxGetPr(prhs[3]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalGradientsMex",
                                  "Requires two outputs.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        dx = mxGetPr(plhs[0]);
        dy = mxGetPr(plhs[1]);

        ap_physical_gradients(C, B, ref_dx, ref_dy, quadPoints, dx, dy);
}
