#include "mex.h"

#include "ArgyrisMaps.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x;
    double *C;
    double *B;
    double *b;
    double *Th;

    // check output.
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:fillArgyrisMaps","Four outputs required.");
    }

    // check input.
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:fillArgyrisMaps","One input required.");
    }

    if (2 != mxGetM(prhs[0]) || 3 != mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:Fill3", "Input must be a 2x3 matrix.");
    }

    // get the coordinates in.
    x = mxGetPr(prhs[0]);

    // create the output matricies.
    plhs[0] = mxCreateDoubleMatrix(21,21,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2,2,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(2,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3,3,mxREAL);

    // create pointers to the output matricies.
    C = mxGetPr(plhs[0]);
    B = mxGetPr(plhs[1]);
    b = mxGetPr(plhs[2]);
    Th = mxGetPr(plhs[3]);

    // fill the output matricies.
    fillArgyrisMaps(x, C, B, b, Th);
}
