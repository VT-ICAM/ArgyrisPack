#include "mex.h"
#include "blas.h"

#include "ArgyrisFunctions.c"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *C;
    double *refValues;
    double *values;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisFunctionsMex",
                          "Requires two arguements.");
    }
    C = mxGetPr(prhs[0]);
    rArgyrisFunctions = mxGetPr(prhs[1]);
    quadPoints = mxGetN(prhs[1]);
    rows = mxGetM(prhs[0]);

    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisFunctionsMex",
                          "The C matrix must be the first arguement.");
    }
    if (21 != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisFunctionsMex",
          "There should be 21 basis functions corresponding to 21 rows.");
    }

    // check output.
    if (1 != nlhs){
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisFunctionsMex",
          "Requires exactly one output.");
    }
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    argyrisFunctions = mxGetPr(plhs[0]);

    // generate output.
    evaluateArgyrisFunctions(C, refValues, values, quadPoints, rows);
}
