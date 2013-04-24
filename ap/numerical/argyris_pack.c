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

#include "ref_values.c"
#include "ref_gradients.c"
#include "ref_hessians.c"

#include "physical_maps.c"
#include "physical_values.c"
#include "physical_gradients.c"
#include "physical_hessians.c"

#include "matrix_mass.c"
#include "matrix_stiffness.c"
#include "matrix_biharmonic.c"
