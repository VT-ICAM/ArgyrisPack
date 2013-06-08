#include <math.h>
#include <string.h>
#include "mex.h"

#include "order_logic.h"
#include "physical_maps.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x;
        double *y;
        double *C;
        double *B;
        double *b;

        /* check input. */
        if(nrhs!=2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalMapsMex",
                                  "Two inputs required.");
        }

        if (3 != mxGetN(prhs[0])*mxGetM(prhs[0]) ||
            3 != mxGetN(prhs[1])*mxGetM(prhs[1])) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalMapsMex",
                                  "Inputs must be 1x3 matricies.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);

        /* check output. */
        if(nlhs!=3) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalMapsMex",
                                  "Three outputs required.");
        }

        plhs[0] = mxCreateDoubleMatrix(21,21,mxREAL);
        plhs[1] = mxCreateDoubleMatrix(2,2,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(2,1,mxREAL);

        C  = mxGetPr(plhs[0]);
        B  = mxGetPr(plhs[1]);
        b  = mxGetPr(plhs[2]);
        ap_physical_maps(x, y, C, B, b);
}
