mex -largeArrayDims CFLAGS='$CFLAGS -Wall -O3' ArgyrisMaps.cpp -lmwblas
mex -largeArrayDims CFLAGS='$CFLAGS -Wall -O3' ArgyrisGradients.cpp -lmwblas
mex -largeArrayDims CFLAGS='$CFLAGS -Wall -O3' ArgyrisHessians.cpp -lmwblas