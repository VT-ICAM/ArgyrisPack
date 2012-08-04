#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "matrix_biharmonic.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *C, *B, *Th;
        double *ref_dxx, *ref_dxy, *ref_dyy;
        double* weights;
        double* biharmonic;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 7) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "Requires seven arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "The B matrix must be the second argument.");
        }
        if ((3 != mxGetN(prhs[2])) || (3 != mxGetM(prhs[2]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "The Th matrix must be the third argument.");
        }
        for (i = 3; i < 5; i++) {
                if (mxGetN(prhs[i]) != mxGetN(prhs[6])) {
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

        C = mxGetPr(prhs[0]);
        B = mxGetPr(prhs[1]);
        Th = mxGetPr(prhs[2]);
        ref_dxx = mxGetPr(prhs[3]);
        ref_dxy = mxGetPr(prhs[4]);
        ref_dyy = mxGetPr(prhs[5]);
        weights = mxGetPr(prhs[6]);
        quadPoints = mxGetN(prhs[3]);

        /* check output. */
        if (nlhs != 1) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apMatrixBiharmonicMex",
                                  "Requires one output.");
        }

        plhs[0]   = mxCreateDoubleMatrix(21, 21, mxREAL);
        biharmonic = mxGetPr(plhs[0]);
        ap_matrix_biharmonic(C, B, Th, ref_dxx, ref_dxy, ref_dyy,
                             weights, quadPoints, biharmonic);
}
