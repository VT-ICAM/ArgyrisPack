% Build file for MEX files.
cd ./ap/mex/
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apLocalFunctionsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apLocalGradientsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apLocalHessiansMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apGlobalFunctionsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apGlobalGradientsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apGlobalHessiansMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apGlobalMapsMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apMatrixMassMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apMatrixStiffnessMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apMatrixBetaplaneMex.c -lmwblas
mex -largeArrayDims -DLAPACKINDEX=mwSignedIndex CFLAGS='-fPIC -O3 -DUSE_COL_MAJOR -std=c99 -D_GNU_SOURCE -pthread -m64 -fexceptions -D_FILE_OFFSET_BITS=64 -I../numerical/' apMatrixBiharmonicMex.c -lmwblas
cd ../../
