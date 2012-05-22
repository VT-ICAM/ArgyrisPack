mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='$CFLAGS -Wall -O3 -std=c99' ArgyrisMapsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='$CFLAGS -Wall -O3 -std=c99' ArgyrisFunctionsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='$CFLAGS -Wall -O3 -std=c99' ArgyrisGradientsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='$CFLAGS -Wall -O3 -std=c99' ArgyrisHessiansMex.c -lmwblas