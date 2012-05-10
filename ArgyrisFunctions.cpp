#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

#define columnOrder(row,col,entriesPerColumn) col * entriesPerColumn + row

void evaluateArgyrisFunctions(double *C, double *rArgyrisFunctions,
                              double *argyrisFunctions,
                              mwSignedIndex quadPoints, mwSignedIndex rows)
{
    // stuff for DGEMM
    char *chn = "N";
    double one = 1.0, zero = 0.0;

    // perform the transformation using the C matrix.
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          rArgyrisFunctions, &rows, &zero, argyrisFunctions, &rows);
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *C;
    double *rArgyrisFunctions;
    double *argyrisFunctions;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisFunctions",
                          "Requires two arguements.");
    }
    C = mxGetPr(prhs[0]);
    rArgyrisFunctions = mxGetPr(prhs[1]);
    quadPoints = mxGetN(prhs[1]);
    rows = mxGetM(prhs[0]);

    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0])))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisFunctions",
                          "The C matrix must be the first arguement.");
    }
    if (21 != mxGetM(prhs[1]))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisFunctions",
          "There should be 21 basis functions corresponding to 21 rows.");
    }

    // check output.
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    argyrisFunctions = mxGetPr(plhs[0]);

    // generate output.
    evaluateArgyrisFunctions(C, rArgyrisFunctions, argyrisFunctions, quadPoints, rows);
}
