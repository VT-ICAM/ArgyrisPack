void ap_global_maps(double* restrict x, double* restrict y,
                    double* restrict C, double* restrict B,
                    double* restrict b, double* restrict Th)
{
        /* temporary values. */
        double x0, x1, x2, y0, y1, y2;
        double B00, B01, B10, B11;
        double w00, w01, w02, w10, w11, w12, w20, w21, w22;
        double v00, v01, v02, v10, v11, v12;
        double norm0, norm1, norm2;
        double norm0squared, norm1squared, norm2squared;
        double C_constant0, C_constant1, C_constant2;

        /* extract coordinates. */
        x0 = x[0];
        x1 = x[1];
        x2 = x[2];
        y0 = y[0];
        y1 = y[1];
        y2 = y[2];

        /* fill B. */
        B00 = -x0 + x1;
        B01 = -x0 + x2;
        B10 = -y0 + y1;
        B11 = -y0 + y2;

        B[ORDER(0, 0, 2, 2)] = B00;
        B[ORDER(0, 1, 2, 2)] = B01;
        B[ORDER(1, 0, 2, 2)] = B10;
        B[ORDER(1, 1, 2, 2)] = B11;

        /* fill b. */
        b[ORDER(0,0,2,1)] = x0;
        b[ORDER(1,0,2,1)] = y0;

        /* v values */
        v00 = -x0 + x1;
        v01 = -x0 + x2;
        v02 = -x1 + x2;
        v10 = -y0 + y1;
        v11 = -y0 + y2;
        v12 = -y1 + y2;

        /* w values */
        w00 = v00*v00;
        w10 = v01*v01;
        w20 = v02*v02;
        w01 = 2*v10*v00;
        w11 = 2*v11*v01;
        w21 = 2*v12*v02;
        w02 = v10*v10;
        w12 = v11*v11;
        w22 = v12*v12;

        /* Manually folded constants. */
        norm0 = sqrt(w00 + w02);
        norm1 = sqrt(w10 + w12);
        norm2 = sqrt(w20 + w22);
        norm0squared = norm0*norm0;
        norm1squared = norm1*norm1;
        norm2squared = norm2*norm2;
        C_constant0 = (B00*v02 + B01*v02 + B10*v12 + B11*v12);
        C_constant1 = (B00*v01 + B10*v11);
        C_constant2 = (B01*v00 + B11*v10);

        /* fill Th. Note that this is the transpose of the usual definition. */
        Th[ORDER(0, 0, 3, 3)] = B00*B00;
        Th[ORDER(1, 0, 3, 3)] = 2*B00*B10;
        Th[ORDER(2, 0, 3, 3)] = B10*B10;
        Th[ORDER(0, 1, 3, 3)] = B00*B01;
        Th[ORDER(1, 1, 3, 3)] = B00*B11 + B01*B10;
        Th[ORDER(2, 1, 3, 3)] = B10*B11;
        Th[ORDER(0, 2, 3, 3)] = B01*B01;
        Th[ORDER(1, 2, 3, 3)] = 2*B01*B11;
        Th[ORDER(2, 2, 3, 3)] = B11*B11;

        /* fill C. Note that this is the transpose of the usual definition. */
        C[ORDER(0, 0, 21, 21)]   = 1;
        C[ORDER(1, 1, 21, 21)]   = 1;
        C[ORDER(2, 2, 21, 21)]   = 1;
        C[ORDER(3, 3, 21, 21)]   = B00;
        C[ORDER(4, 3, 21, 21)]   = B10;
        C[ORDER(3, 4, 21, 21)]   = B01;
        C[ORDER(4, 4, 21, 21)]   = B11;
        C[ORDER(5, 5, 21, 21)]   = B00;
        C[ORDER(6, 5, 21, 21)]   = B10;
        C[ORDER(5, 6, 21, 21)]   = B01;
        C[ORDER(6, 6, 21, 21)]   = B11;
        C[ORDER(7, 7, 21, 21)]   = B00;
        C[ORDER(8, 7, 21, 21)]   = B10;
        C[ORDER(7, 8, 21, 21)]   = B01;
        C[ORDER(8, 8, 21, 21)]   = B11;
        C[ORDER(9, 9, 21, 21)]   = B00*B00;
        C[ORDER(10, 9, 21, 21)]  = 2*B00*B10;
        C[ORDER(11, 9, 21, 21)]  = B10*B10;
        C[ORDER(9, 10, 21, 21)]  = B00*B01;
        C[ORDER(10, 10, 21, 21)] = B00*B11 + B01*B10;
        C[ORDER(11, 10, 21, 21)] = B10*B11;
        C[ORDER(9, 11, 21, 21)]  = B01*B01;
        C[ORDER(10, 11, 21, 21)] = 2*B01*B11;
        C[ORDER(11, 11, 21, 21)] = B11*B11;
        C[ORDER(12, 12, 21, 21)] = B00*B00;
        C[ORDER(13, 12, 21, 21)] = 2*B00*B10;
        C[ORDER(14, 12, 21, 21)] = B10*B10;
        C[ORDER(12, 13, 21, 21)] = B00*B01;
        C[ORDER(13, 13, 21, 21)] = B00*B11 + B01*B10;
        C[ORDER(14, 13, 21, 21)] = B10*B11;
        C[ORDER(12, 14, 21, 21)] = B01*B01;
        C[ORDER(13, 14, 21, 21)] = 2*B01*B11;
        C[ORDER(14, 14, 21, 21)] = B11*B11;
        C[ORDER(15, 15, 21, 21)] = B00*B00;
        C[ORDER(16, 15, 21, 21)] = 2*B00*B10;
        C[ORDER(17, 15, 21, 21)] = B10*B10;
        C[ORDER(15, 16, 21, 21)] = B00*B01;
        C[ORDER(16, 16, 21, 21)] = B00*B11 + B01*B10;
        C[ORDER(17, 16, 21, 21)] = B10*B11;
        C[ORDER(15, 17, 21, 21)] = B01*B01;
        C[ORDER(16, 17, 21, 21)] = 2*B01*B11;
        C[ORDER(17, 17, 21, 21)] = B11*B11;

        /* row 18 */
        C[ORDER(0, 18, 21, 21)]  = -15.0/8.0*C_constant2/norm0squared;
        C[ORDER(1, 18, 21, 21)]  = 15.0/8.0*C_constant2/norm0squared;
        C[ORDER(3, 18, 21, 21)]  = -7.0/16.0*C_constant2*v00/norm0squared;
        C[ORDER(4, 18, 21, 21)]  = -7.0/16.0*C_constant2*v10/norm0squared;
        C[ORDER(5, 18, 21, 21)]  = -7.0/16.0*C_constant2*v00/norm0squared;
        C[ORDER(6, 18, 21, 21)]  = -7.0/16.0*C_constant2*v10/norm0squared;
        C[ORDER(9, 18, 21, 21)]  = -1.0/32.0*C_constant2*w00/norm0squared;
        C[ORDER(10, 18, 21, 21)] = -1.0/32.0*C_constant2*w01/norm0squared;
        C[ORDER(11, 18, 21, 21)] = -1.0/32.0*C_constant2*w02/norm0squared;
        C[ORDER(12, 18, 21, 21)] = 1.0/32.0*C_constant2*w00/norm0squared;
        C[ORDER(13, 18, 21, 21)] = 1.0/32.0*C_constant2*w01/norm0squared;
        C[ORDER(14, 18, 21, 21)] = 1.0/32.0*C_constant2*w02/norm0squared;
        C[ORDER(18, 18, 21, 21)] = -(B01*v10 - B11*v00)/norm0;

        /* row 19 */
        C[ORDER(0, 19, 21, 21)]  = 15.0/8.0*C_constant1/norm1squared;
        C[ORDER(2, 19, 21, 21)]  = -15.0/8.0*C_constant1/norm1squared;
        C[ORDER(3, 19, 21, 21)]  = 7.0/16.0*C_constant1*v01/norm1squared;
        C[ORDER(4, 19, 21, 21)]  = 7.0/16.0*C_constant1*v11/norm1squared;
        C[ORDER(7, 19, 21, 21)]  = 7.0/16.0*C_constant1*v01/norm1squared;
        C[ORDER(8, 19, 21, 21)]  = 7.0/16.0*C_constant1*v11/norm1squared;
        C[ORDER(9, 19, 21, 21)]  = 1.0/32.0*C_constant1*w10/norm1squared;
        C[ORDER(10, 19, 21, 21)] = 1.0/32.0*C_constant1*w11/norm1squared;
        C[ORDER(11, 19, 21, 21)] = 1.0/32.0*C_constant1*w12/norm1squared;
        C[ORDER(15, 19, 21, 21)] = -1.0/32.0*C_constant1*w10/norm1squared;
        C[ORDER(16, 19, 21, 21)] = -1.0/32.0*C_constant1*w11/norm1squared;
        C[ORDER(17, 19, 21, 21)] = -1.0/32.0*C_constant1*w12/norm1squared;
        C[ORDER(19, 19, 21, 21)] = (B00*v11 - B10*v01)/norm1;

        /* row 20 */
        C[ORDER(1, 20, 21, 21)]  = 15.0/16.0*C_constant0*SQRT2/norm2squared;
        C[ORDER(2, 20, 21, 21)]  = -15.0/16.0*C_constant0*SQRT2/norm2squared;
        C[ORDER(5, 20, 21, 21)]  = 7.0/32.0*C_constant0*SQRT2*v02/norm2squared;
        C[ORDER(6, 20, 21, 21)]  = 7.0/32.0*C_constant0*SQRT2*v12/norm2squared;
        C[ORDER(7, 20, 21, 21)]  = 7.0/32.0*C_constant0*SQRT2*v02/norm2squared;
        C[ORDER(8, 20, 21, 21)]  = 7.0/32.0*C_constant0*SQRT2*v12/norm2squared;
        C[ORDER(12, 20, 21, 21)] = 1.0/64.0*C_constant0*SQRT2*w20/norm2squared;
        C[ORDER(13, 20, 21, 21)] = 1.0/64.0*C_constant0*SQRT2*w21/norm2squared;
        C[ORDER(14, 20, 21, 21)] = 1.0/64.0*C_constant0*SQRT2*w22/norm2squared;
        C[ORDER(15, 20, 21, 21)] = -1.0/64.0*C_constant0*SQRT2*w20/norm2squared;
        C[ORDER(16, 20, 21, 21)] = -1.0/64.0*C_constant0*SQRT2*w21/norm2squared;
        C[ORDER(17, 20, 21, 21)] = -1.0/64.0*C_constant0*SQRT2*w22/norm2squared;
        C[ORDER(20, 20, 21, 21)] = 0.5*(B00*v12 + B01*v12 - B10*v02 - B11*v02)
                *SQRT2/norm2;
}
