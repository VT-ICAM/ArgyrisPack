#include "mex.h"

#include "global_maps.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x;
    double *C;
    double *B;
    double *b;
    double *Th;

    /* check input. */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalMapsMex","One input required.");
    }

    if (2 != mxGetM(prhs[0]) || 3 != mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalMapsMex", "Input must be a 2x3 matrix.");
    }

    x = mxGetPr(prhs[0]);

    /* check output. */
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGlobalMapsMex","Four outputs required.");
    }

    plhs[0] = mxCreateDoubleMatrix(21,21,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2,2,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(2,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3,3,mxREAL);
    C = mxGetPr(plhs[0]);
    B = mxGetPr(plhs[1]);
    b = mxGetPr(plhs[2]);
    Th = mxGetPr(plhs[3]);

    ap_global_maps(x, C, B, b, Th);
}
