#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "order_logic.h"

#include "multiply_by_diagonal.c"
#include "physical_hessians.c"
#include "matrix_biharmonic.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C, *B;
        double *ref_dxx, *ref_dxy, *ref_dyy;
        double* weights;
        double* biharmonic;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 6) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "Requires six arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "The B matrix must be the second argument.");
        }
        for (i = 2; i < 5; i++) {
                if (mxGetN(prhs[i]) != mxGetN(prhs[2])) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                          "Mismatch in number of quadrature "
                                          "points.");
                }

                if (mxGetM(prhs[i]) != 21) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                          "There should be 21 basis functions "
                                          "corresponding to 21 rows.");
                }

        }
        if (mxGetN(prhs[5]) != mxGetN(prhs[2]) &&
            mxGetM(prhs[5]) != mxGetN(prhs[2])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "Mismatch in number of weights "
                                  "and number of quadrature points.");
        }
        if (mxGetM(prhs[5]) != 1 && mxGetN(prhs[5]) != 1) {
            mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                              "The weights must be in a one-dimensional array.");
        }

        C          = mxGetPr(prhs[0]);
        B          = mxGetPr(prhs[1]);
        ref_dxx    = mxGetPr(prhs[2]);
        ref_dxy    = mxGetPr(prhs[3]);
        ref_dyy    = mxGetPr(prhs[4]);
        weights    = mxGetPr(prhs[5]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "Requires one output.");
        }

        plhs[0]    = mxCreateDoubleMatrix(21, 21, mxREAL);
        biharmonic = mxGetPr(plhs[0]);

        ap_matrix_biharmonic(C, B, ref_dxx, ref_dxy, ref_dyy, weights,
                             quadPoints, biharmonic);
}
