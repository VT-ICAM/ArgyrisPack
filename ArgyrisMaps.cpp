#include <math.h>
#include "mex.h"

#define columnOrder(row,col,entriesPerColumn) col * entriesPerColumn + row

void fillArgyrisMaps(double* x, double* C, double* B, double* b, double* Th)
{
    // temporary values.
    double x0, x1, x2, y0, y1, y2;
    double B00, B01, B10, B11;
    double w00, w01, w02, w10, w11, w12, w20, w21, w22;
    double v00, v01, v02, v10, v11, v12;
    double norm0, norm1, norm2;
    double norm0squared, norm1squared, norm2squared;
    double CConstant0, CConstant1, CConstant2;
    double SQRT2;

    // extract coordinates.
    x0 = x[columnOrder(0, 0, 2)];
    y0 = x[columnOrder(1, 0, 2)];
    x1 = x[columnOrder(0, 1, 2)];
    y1 = x[columnOrder(1, 1, 2)];
    x2 = x[columnOrder(0, 2, 2)];
    y2 = x[columnOrder(1, 2, 2)];

    // fill B.
    B00 = -x0 + x1;
    B01 = -x0 + x2;
    B10 = -y0 + y1;
    B11 = -y0 + y2;

    B[columnOrder(0, 0, 2)] = B00;
    B[columnOrder(0, 1, 2)] = B01;
    B[columnOrder(1, 0, 2)] = B10;
    B[columnOrder(1, 1, 2)] = B11;

    // fill b.
    b[0] = x0;
    b[1] = y0;

    // v values
    v00 = -x0 + x1;
    v01 = -x0 + x2;
    v02 = -x1 + x2;
    v10 = -y0 + y1;
    v11 = -y0 + y2;
    v12 = -y1 + y2;

    // w values
    w00 = v00*v00;
    w10 = v01*v01;
    w20 = v02*v02;
    w01 = 2*v10*v00;
    w11 = 2*v11*v01;
    w21 = 2*v12*v02;
    w02 = v10*v10;
    w12 = v11*v11;
    w22 = v12*v12;

    // constants.
    norm0 = sqrt(w00 + w02);
    norm1 = sqrt(w10 + w12);
    norm2 = sqrt(w20 + w22);
    SQRT2 = sqrt(2.0);
    norm0squared = norm0*norm0;
    norm1squared = norm1*norm1;
    norm2squared = norm2*norm2;
    CConstant0 = (B00*v02 + B01*v02 + B10*v12 + B11*v12);
    CConstant1 = (B00*v01 + B10*v11);
    CConstant2 = (B01*v00 + B11*v10);

    // fill Th. Note that this is the transpose of the usual definition.
    Th[columnOrder(0, 0, 3)] = B00*B00;
    Th[columnOrder(1, 0, 3)] = 2*B00*B10;
    Th[columnOrder(2, 0, 3)] = B10*B10;
    Th[columnOrder(0, 1, 3)] = B00*B01;
    Th[columnOrder(1, 1, 3)] = B00*B11 + B01*B10;
    Th[columnOrder(2, 1, 3)] = B10*B11;
    Th[columnOrder(0, 2, 3)] = B01*B01;
    Th[columnOrder(1, 2, 3)] = 2*B01*B11;
    Th[columnOrder(2, 2, 3)] = B11*B11;

    // fill C. Note that this is the transpose of the usual definition.
    C[columnOrder(0, 0, 21)]   = 1;
    C[columnOrder(1, 1, 21)]   = 1;
    C[columnOrder(2, 2, 21)]   = 1;
    C[columnOrder(3, 3, 21)]   = B00;
    C[columnOrder(4, 3, 21)]   = B10;
    C[columnOrder(3, 4, 21)]   = B01;
    C[columnOrder(4, 4, 21)]   = B11;
    C[columnOrder(5, 5, 21)]   = B00;
    C[columnOrder(6, 5, 21)]   = B10;
    C[columnOrder(5, 6, 21)]   = B01;
    C[columnOrder(6, 6, 21)]   = B11;
    C[columnOrder(7, 7, 21)]   = B00;
    C[columnOrder(8, 7, 21)]   = B10;
    C[columnOrder(7, 8, 21)]   = B01;
    C[columnOrder(8, 8, 21)]   = B11;
    C[columnOrder(9, 9, 21)]   = B00*B00;
    C[columnOrder(10, 9, 21)]  = 2*B00*B10;
    C[columnOrder(11, 9, 21)]  = B10*B10;
    C[columnOrder(9, 10, 21)]  = B00*B01;
    C[columnOrder(10, 10, 21)] = B00*B11 + B01*B10;
    C[columnOrder(11, 10, 21)] = B10*B11;
    C[columnOrder(9, 11, 21)]  = B01*B01;
    C[columnOrder(10, 11, 21)] = 2*B01*B11;
    C[columnOrder(11, 11, 21)] = B11*B11;
    C[columnOrder(12, 12, 21)] = B00*B00;
    C[columnOrder(13, 12, 21)] = 2*B00*B10;
    C[columnOrder(14, 12, 21)] = B10*B10;
    C[columnOrder(12, 13, 21)] = B00*B01;
    C[columnOrder(13, 13, 21)] = B00*B11 + B01*B10;
    C[columnOrder(14, 13, 21)] = B10*B11;
    C[columnOrder(12, 14, 21)] = B01*B01;
    C[columnOrder(13, 14, 21)] = 2*B01*B11;
    C[columnOrder(14, 14, 21)] = B11*B11;
    C[columnOrder(15, 15, 21)] = B00*B00;
    C[columnOrder(16, 15, 21)] = 2*B00*B10;
    C[columnOrder(17, 15, 21)] = B10*B10;
    C[columnOrder(15, 16, 21)] = B00*B01;
    C[columnOrder(16, 16, 21)] = B00*B11 + B01*B10;
    C[columnOrder(17, 16, 21)] = B10*B11;
    C[columnOrder(15, 17, 21)] = B01*B01;
    C[columnOrder(16, 17, 21)] = 2*B01*B11;
    C[columnOrder(17, 17, 21)] = B11*B11;

    // row 18
    C[columnOrder(0, 18, 21)]  = -15.0/8.0*CConstant2/norm0squared;
    C[columnOrder(1, 18, 21)]  = 15.0/8.0*CConstant2/norm0squared;
    C[columnOrder(3, 18, 21)]  = -7.0/16.0*CConstant2*v00/norm0squared;
    C[columnOrder(4, 18, 21)]  = -7.0/16.0*CConstant2*v10/norm0squared;
    C[columnOrder(5, 18, 21)]  = -7.0/16.0*CConstant2*v00/norm0squared;
    C[columnOrder(6, 18, 21)]  = -7.0/16.0*CConstant2*v10/norm0squared;
    C[columnOrder(9, 18, 21)]  = -1.0/32.0*CConstant2*w00/norm0squared;
    C[columnOrder(10, 18, 21)] = -1.0/32.0*CConstant2*w01/norm0squared;
    C[columnOrder(11, 18, 21)] = -1.0/32.0*CConstant2*w02/norm0squared;
    C[columnOrder(12, 18, 21)] = 1.0/32.0*CConstant2*w00/norm0squared;
    C[columnOrder(13, 18, 21)] = 1.0/32.0*CConstant2*w01/norm0squared;
    C[columnOrder(14, 18, 21)] = 1.0/32.0*CConstant2*w02/norm0squared;
    C[columnOrder(18, 18, 21)] = -(B01*v10 - B11*v00)/norm0;

    // row 19
    C[columnOrder(0, 19, 21)]  = 15.0/8.0*CConstant1/norm1squared;
    C[columnOrder(2, 19, 21)]  = -15.0/8.0*CConstant1/norm1squared;
    C[columnOrder(3, 19, 21)]  = 7.0/16.0*CConstant1*v01/norm1squared;
    C[columnOrder(4, 19, 21)]  = 7.0/16.0*CConstant1*v11/norm1squared;
    C[columnOrder(7, 19, 21)]  = 7.0/16.0*CConstant1*v01/norm1squared;
    C[columnOrder(8, 19, 21)]  = 7.0/16.0*CConstant1*v11/norm1squared;
    C[columnOrder(9, 19, 21)]  = 1.0/32.0*CConstant1*w10/norm1squared;
    C[columnOrder(10, 19, 21)] = 1.0/32.0*CConstant1*w11/norm1squared;
    C[columnOrder(11, 19, 21)] = 1.0/32.0*CConstant1*w12/norm1squared;
    C[columnOrder(15, 19, 21)] = -1.0/32.0*CConstant1*w10/norm1squared;
    C[columnOrder(16, 19, 21)] = -1.0/32.0*CConstant1*w11/norm1squared;
    C[columnOrder(17, 19, 21)] = -1.0/32.0*CConstant1*w12/norm1squared;
    C[columnOrder(19, 19, 21)] = (B00*v11 - B10*v01)/norm1;

    // row 20
    C[columnOrder(1, 20, 21)]  = 15.0/16.0*CConstant0*SQRT2/norm2squared;
    C[columnOrder(2, 20, 21)]  = -15.0/16.0*CConstant0*SQRT2/norm2squared;
    C[columnOrder(5, 20, 21)]  = 7.0/32.0*CConstant0*SQRT2*v02/norm2squared;
    C[columnOrder(6, 20, 21)]  = 7.0/32.0*CConstant0*SQRT2*v12/norm2squared;
    C[columnOrder(7, 20, 21)]  = 7.0/32.0*CConstant0*SQRT2*v02/norm2squared;
    C[columnOrder(8, 20, 21)]  = 7.0/32.0*CConstant0*SQRT2*v12/norm2squared;
    C[columnOrder(12, 20, 21)] = 1.0/64.0*CConstant0*SQRT2*w20/norm2squared;
    C[columnOrder(13, 20, 21)] = 1.0/64.0*CConstant0*SQRT2*w21/norm2squared;
    C[columnOrder(14, 20, 21)] = 1.0/64.0*CConstant0*SQRT2*w22/norm2squared;
    C[columnOrder(15, 20, 21)] = -1.0/64.0*CConstant0*SQRT2*w20/norm2squared;
    C[columnOrder(16, 20, 21)] = -1.0/64.0*CConstant0*SQRT2*w21/norm2squared;
    C[columnOrder(17, 20, 21)] = -1.0/64.0*CConstant0*SQRT2*w22/norm2squared;
    C[columnOrder(20, 20, 21)] = 0.5*(B00*v12 + B01*v12 - B10*v02 - B11*v02)*SQRT2/norm2;
}

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

    if (2 != mxGetM(prhs[0]) || 3 != mxGetN(prhs[0]))
    {
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
