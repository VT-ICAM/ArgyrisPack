void ap_affine_transformation(double* restrict B, double* restrict b,
                              double* restrict ref_x, double* restrict ref_y,
                              int length, double* restrict physical_x,
                              double* restrict physical_y)
{
/*
 * Use the affine map described by B and b to transform reference coordinates
 * into physical coordinates.
 */
        int i;
        for (i = 0; i < length; i++) {
                physical_x[i] = B[ORDER(0, 0, 2, 2)]*ref_x[i]
                              + B[ORDER(0, 1, 2, 2)]*ref_y[i] + b[0];
                physical_y[i] = B[ORDER(1, 0, 2, 2)]*ref_x[i]
                              + B[ORDER(1, 1, 2, 2)]*ref_y[i] + b[1];
        }
}

void ap_inverse_affine_transformation(double* restrict B, double* restrict b,
                                      double* restrict physical_x,
                                      double* restrict physical_y, int length,
                                      double* restrict ref_x,
                                      double* restrict ref_y)
{
/*
 * Use the inverse of the affine map described by B and b to transform physical
 * coordinates into reference coordinates.
 */
        int i;
        double determinant = B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                             B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)];
        for (i = 0; i < length; i++) {
                ref_x[i] = (B[ORDER(1, 1, 2, 2)]*(physical_x[i] - b[0])
                     - B[ORDER(0, 1, 2, 2)]*(physical_y[i] - b[1]))/determinant;
                ref_y[i] = (-1.0*B[ORDER(1, 0, 2, 2)]*(physical_x[i] - b[0])
                     + B[ORDER(0, 0, 2, 2)]*(physical_y[i] - b[1]))/determinant;
        }
}
