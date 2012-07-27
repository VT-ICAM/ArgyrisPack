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
double __d_zero = 0;
double __d_one = 1;

#ifdef USE_ROW_MAJOR
        #define ORDER(row,col,nrows,ncols) (row)*(ncols) + (col)
        #define DGEMM_WRAPPER(m, n, k, A, B, C) \
        dgemm(&__c_N, &__c_N, &n, &m, &k, &__d_one, B, &n, A, &k, &__d_zero, C, &n)
#endif

#ifdef USE_COL_MAJOR
        #define ORDER(row,col,nrows,ncols) (col)*(nrows) + (row)
        #define DGEMM_WRAPPER(m, n, k, A, B, C) \
        dgemm(&__c_N, &__c_N, &m, &n, &k, &__d_one, A, &m, B, &k, &__d_zero, C, &m)
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
