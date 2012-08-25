#define __ORDER_LOGIC
#define SQRT2 1.4142135623730950

/*
 * Some users need column order (e.g. MATLAB, Fortran) while others (e.g. Python)
 * need row order. Therefore handle all array access and multiplication via macros
 * to 'do the right thing'. The code handling dgemm is ugly because all arguements
 * must be passed by reference (and therefore explicitly declared at the calling
 * scope).
 */

/* constants for DGEMM_WRAPPER */
char __c_N = 'N';
char __c_T = 'T';
double __d_zero = 0;
double __d_one = 1;

/* The DGEMM_WRAPPERs calls vary based on whether or not there is a transpose
 * going on; the dimensions supplied in the first three arguments are the
 * dimensions of op(A) and op(B), just like the original dgemm() call.
 */
#ifdef USE_ROW_MAJOR
        #define ORDER(row,col,nrows,ncols) (row)*(ncols) + (col)
        /* wrap C.T := B.T * A.T */
        #define DGEMM_WRAPPER(m, n, k, D, E, F) \
        dgemm(&(__c_N), &(__c_N), &(n), &(m), &(k), &(__d_one), E, &(n), D, \
             &(k), &(__d_zero), F, &(n))
        /* wrap C.T := B * A.T (or C = A * B.T).
         *
         * m : number of rows of A; n : number of columns of B.T; k : number of
         * columns of A / rows of B.T. For sanity, define D, E, and F to be the
         * transposes of A, B, and C.
         *
         * Therefore E.T is n x k, D is k x m, and F is n x m
         */
        #define DGEMM_WRAPPER_NT(m, n, k, D, E, F) \
        dgemm(&(__c_T), &(__c_N), &(n), &(m), &(k), &(__d_one), E, &(k), D, \
              &(k), &(__d_zero), F, &(n))
        #define DGEMM_WRAPPER_NT_ADD_C(m, n, k, D, E, F) \
        dgemm(&(__c_T), &(__c_N), &(n), &(m), &(k), &(__d_one), E, &(k), D, \
              &(k), &(__d_one), F, &(n))
#endif

#ifdef USE_COL_MAJOR
        #define ORDER(row,col,nrows,ncols) (col)*(nrows) + (row)
        #define DGEMM_WRAPPER(m, n, k, A, B, C) \
        dgemm(&(__c_N), &(__c_N), &(m), &(n), &(k), &(__d_one), A, &(m), B, \
              &(k), &(__d_zero), C, &(m))
        /* wrap C := A*B.T */
        #define DGEMM_WRAPPER_NT(m, n, k, A, B, C) \
        dgemm(&(__c_N), &(__c_T), &(m), &(n), &(k), &(__d_one), A, &(m), B, \
              &(n), &(__d_zero), C, &(m))
        /* wrap C.T := A*B.T + C */
        #define DGEMM_WRAPPER_NT_ADD_C(m, n, k, A, B, C) \
        dgemm(&(__c_N), &(__c_T), &(m), &(n), &(k), &(__d_one), A, &(m), B, \
              &(n), &(__d_one), C, &(m))
#endif

#ifdef USE_ROW_MAJOR
        #ifdef USE_COL_MAJOR
                #error "Must define either USE_COL_MAJOR or USE_ROW_MAJOR for "
                        "storage order, not both"
        #endif
#endif

#ifndef USE_ROW_MAJOR
        #ifndef USE_COL_MAJOR
                #error "Must define either USE_COL_MAJOR or USE_ROW_MAJOR for "
                       "storage order"
        #endif
#endif
