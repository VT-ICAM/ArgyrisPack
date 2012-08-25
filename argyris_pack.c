#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 * dgemm_ appears to be more portable than dgemm (the version of BLAS with OSX
 * and netlib blas both support this naming convention).
 */
#define dgemm dgemm_
void dgemm(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*,
           double*, double*, int*);

/*
 * Most of these functions are short. Include them in this order to avoid a lot
 * of hard work with header files and dependencies.
 */

#include "order_logic.h"
#include "multiply_by_diagonal.c"

#include "local_functions.c"
#include "local_gradients.c"
#include "local_hessians.c"

#include "global_maps.c"
#include "global_functions.c"
#include "global_gradients.c"
#include "global_hessians.c"

#include "matrix_mass.c"
#include "matrix_stiffness.c"
#include "matrix_biharmonic.c"
