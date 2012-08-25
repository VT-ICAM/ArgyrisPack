/* a test file for ArgyrisPack */
#include <stdio.h>
#include "argyris_pack.h"
#include "order_logic.h"

void print_matrix(double* matrix, int rows, int cols)
{
        int i, j;
        for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                        printf("%.16f ", matrix[ORDER(i,j,rows,cols)]);
                }
                printf("\n");
        }
}

void print_matrix_flat(double* matrix, int rows, int cols)
{
        int i;
        for (i = 0; i < rows*cols; i++) {
                printf("%.16f ", matrix[i]);
        }
        printf("\n");
}

int main(int argc, char *argv[])
{
        double x[3] = {0.1, 0.2, 0.3};
        double y[3] = {0.2, 0.2, 0.3};

        double ref_dx[21*3];
        double ref_dy[21*3];

        ap_local_gradients(x, y, 3, ref_dx, ref_dy);

        printf("ref_dx: ");
        print_matrix(ref_dx, 21, 3);
        printf("ref_dy: ");
        print_matrix(ref_dy, 21, 3);


        return 0;
}
