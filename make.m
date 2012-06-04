mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' ArgyrisMapsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' ArgyrisFunctionsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' ArgyrisGradientsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64'ArgyrisHessiansMex.c -lmwblas
