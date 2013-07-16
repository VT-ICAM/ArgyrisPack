#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "order_logic.h"

#include "diagonal_multiply.c"
#include "physical_gradients.c"
#include "matrix_stiffness.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C, *B;
        double *ref_dx, *ref_dy;
        double* weights;
        double* stiffness;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 5) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                  "Requires five arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                  "The B matrix must be the second argument.");
        }
        for (i = 2; i < 4; i++) {
                if (mxGetN(prhs[i]) != mxGetN(prhs[2])) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                          "Mismatch in number of quadrature "
                                          "points.");
                }

                if (mxGetM(prhs[i]) != 21) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                          "There should be 21 basis functions "
                                          "corresponding to 21 rows.");
                }

        }
        if (mxGetN(prhs[2]) != mxGetN(prhs[4]) &&
            mxGetN(prhs[2]) != mxGetM(prhs[4])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                  "Mismatch in number of weights "
                                  "and number of quadrature points.");
        }
        if (mxGetM(prhs[4]) != 1 && mxGetN(prhs[4]) != 1) {
            mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                              "The weights must be in a one-dimensional array.");
        }

        C          = mxGetPr(prhs[0]);
        B          = mxGetPr(prhs[1]);
        ref_dx     = mxGetPr(prhs[2]);
        ref_dy     = mxGetPr(prhs[3]);
        weights    = mxGetPr(prhs[4]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixStiffnessMex",
                                  "Requires one output.");
        }

        plhs[0]   = mxCreateDoubleMatrix(21, 21, mxREAL);
        stiffness = mxGetPr(plhs[0]);
        ap_matrix_stiffness(C, B, ref_dx, ref_dy, weights, quadPoints,
                            stiffness);
}
