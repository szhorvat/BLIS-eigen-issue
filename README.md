# BLIS-eigen-issue

A small C example that computes the eigenvalues and eigenvectors of the
specified 4-by-4 symmetric matrix with LAPACK's `dsyevr` routine.

## Build and run

```sh
cmake -S . -B build
cmake --build build
./build/eigen_example
```

CMake locates and links the LAPACK implementation (and its BLAS dependency).
To choose an implementation supported by CMake's `FindLAPACK`, pass its vendor
name when configuring, for example:

```sh
cmake -S . -B build -DBLA_VENDOR=OpenBLAS
```

The output lists the eigenvalues and eigenvectors, with each eigenvector in a
column.