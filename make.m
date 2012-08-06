mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apGlobalFunctionsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apGlobalGradientsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apGlobalHessiansMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apGlobalMapsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apMatrixBiharmonicMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apMatrixMassMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-O3 -mtune=native -DUSE_COL_MAJOR -ansi -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64' apMatrixStiffnessMex.c -lmwblas
