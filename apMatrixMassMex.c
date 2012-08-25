#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "order_logic.h"

#include "multiply_by_diagonal.c"
#include "global_functions.c"
#include "matrix_mass.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C, *B;
        double *ref_values;
        double* weights;
        double* mass;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 4) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                                  "Requires four arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                                  "The B matrix must be the second argument.");
        }
        if (mxGetM(prhs[2]) != 21) {
            mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                              "There should be 21 basis functions "
                              "corresponding to 21 rows.");
        }
        if (mxGetN(prhs[2]) != mxGetN(prhs[3]) &&
            mxGetN(prhs[2]) != mxGetM(prhs[3])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                                  "Mismatch in number of weights "
                                  "and number of quadrature points.");
        }
        if (mxGetM(prhs[3]) != 1 && mxGetN(prhs[3]) != 1) {
            mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                              "The weights must be in a one-dimensional array.");
        }

        C          = mxGetPr(prhs[0]);
        B          = mxGetPr(prhs[1]);
        ref_values = mxGetPr(prhs[2]);
        weights    = mxGetPr(prhs[3]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixMassMex",
                                  "Requires one output.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, 21, mxREAL);
        mass    = mxGetPr(plhs[0]);
        ap_matrix_mass(C, B, ref_values, weights, quadPoints, mass);
}
