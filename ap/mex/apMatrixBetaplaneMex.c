#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "order_logic.h"

#include "multiply_by_diagonal.c"
#include "physical_gradients.c"
#include "physical_values.c"
#include "matrix_betaplane.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C, *B;
        double *ref_values, *ref_dx, *ref_dy;
        double* weights;
        double* betaplane;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 6) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                  "Requires six arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                  "The B matrix must be the second argument.");
        }
        for (i = 2; i < 5; i++) {
                if (mxGetN(prhs[i]) != mxGetN(prhs[2])) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                          "Mismatch in number of quadrature "
                                          "points.");
                }

                if (mxGetM(prhs[i]) != 21) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                          "There should be 21 basis functions "
                                          "corresponding to 21 rows.");
                }

        }
        if (mxGetN(prhs[2]) != mxGetN(prhs[5]) &&
            mxGetN(prhs[2]) != mxGetM(prhs[5])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                  "Mismatch in number of weights "
                                  "and number of quadrature points.");
        }
        if (mxGetM(prhs[5]) != 1 && mxGetN(prhs[5]) != 1) {
            mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                              "The weights must be in a one-dimensional array.");
        }

        C          = mxGetPr(prhs[0]);
        B          = mxGetPr(prhs[1]);
        ref_values = mxGetPr(prhs[2]);
        ref_dx     = mxGetPr(prhs[3]);
        ref_dy     = mxGetPr(prhs[4]);
        weights    = mxGetPr(prhs[5]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBetaplaneMex",
                                  "Requires one output.");
        }

        plhs[0]   = mxCreateDoubleMatrix(21, 21, mxREAL);
        betaplane = mxGetPr(plhs[0]);
        ap_matrix_betaplane(C, B, ref_values, ref_dx, ref_dy, weights,
                            quadPoints, betaplane);
}
