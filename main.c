#include <stdio.h>
#include <stdlib.h>

#define MATRIX_SIZE 4
/* Set to 0 to disable dividing the matrix by four. */
#ifndef DIVIDE_MATRIX_BY_FOUR
#define DIVIDE_MATRIX_BY_FOUR 1
#endif

extern void dsyevr_(char *jobz, char *range, char *uplo, int *n,
                    double *a, int *lda, double *vl, double *vu, int *il,
                    int *iu, double *abstol, int *m, double *w, double *z,
                    int *ldz, int *isuppz, double *work, int *lwork,
                    int *iwork, int *liwork, int *info);

enum compute_status {
    COMPUTE_SUCCESS,
    COMPUTE_WORKSPACE_QUERY_ERROR,
    COMPUTE_WORKSPACE_ALLOCATION_ERROR,
    COMPUTE_EIGENSOLVE_ERROR,
    COMPUTE_EIGENVALUE_COUNT_ERROR
};

static enum compute_status compute_eigensystem(
    const double input[MATRIX_SIZE * MATRIX_SIZE],
    double eigenvalues[MATRIX_SIZE],
    double eigenvectors[MATRIX_SIZE * MATRIX_SIZE],
    int *info,
    int *found)
{
    int n = MATRIX_SIZE;
    int lda = MATRIX_SIZE;
    int ldz = MATRIX_SIZE;
    /* LAPACK stores matrices column-major. */
    double matrix[MATRIX_SIZE * MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        matrix[i] = input[i];
    }

    int isuppz[2 * MATRIX_SIZE];
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    double vl = 0.0;
    double vu = 0.0;
    int il = 0;
    int iu = 0;
    double abstol = 0.0;
    double work_query = 0.0;
    int iwork_query = 0;
    int lwork = -1;
    int liwork = -1;

    /* Query workspace sizes first; -1 asks for sizes without computing. */
    dsyevr_(&jobz, &range, &uplo, &n, matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, found, eigenvalues, eigenvectors, &ldz, isuppz,
            &work_query, &lwork, &iwork_query, &liwork, info);
    if (*info != 0) {
        return COMPUTE_WORKSPACE_QUERY_ERROR;
    }

    lwork = (int)work_query;
    liwork = iwork_query;
    double *work = malloc((size_t)lwork * sizeof(*work));
    int *iwork = malloc((size_t)liwork * sizeof(*iwork));
    if (work == NULL || iwork == NULL) {
        free(work);
        free(iwork);
        return COMPUTE_WORKSPACE_ALLOCATION_ERROR;
    }

    /* Call again with the allocated workspaces to compute the eigensystem. */
    dsyevr_(&jobz, &range, &uplo, &n, matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, found, eigenvalues, eigenvectors, &ldz, isuppz,
            work, &lwork, iwork, &liwork, info);
    free(work);
    free(iwork);
    if (*info != 0) {
        return COMPUTE_EIGENSOLVE_ERROR;
    }
    if (*found != n) {
        return COMPUTE_EIGENVALUE_COUNT_ERROR;
    }

    return COMPUTE_SUCCESS;
}

static void report_dsyevr_error(const char *stage, int info)
{
    if (info < 0) {
        fprintf(stderr, "dsyevr %s: argument %d has an illegal value\n",
                stage, -info);
    } else if (info > 0) {
        fprintf(stderr,
                "dsyevr %s: internal error in the DSTEMR eigensolver "
                "(info=%d)\n",
                stage, info);
    }
}

static void report_compute_error(enum compute_status status, int info,
                                 int found)
{
    switch (status) {
    case COMPUTE_WORKSPACE_QUERY_ERROR:
        report_dsyevr_error("workspace query", info);
        break;
    case COMPUTE_WORKSPACE_ALLOCATION_ERROR:
        fprintf(stderr, "could not allocate LAPACK workspace\n");
        break;
    case COMPUTE_EIGENSOLVE_ERROR:
        report_dsyevr_error("eigensolve", info);
        break;
    case COMPUTE_EIGENVALUE_COUNT_ERROR:
        fprintf(stderr, "dsyevr returned %d eigenvalues; expected %d\n",
                found, MATRIX_SIZE);
        break;
    case COMPUTE_SUCCESS:
        break;
    }
}

static int next_permutation(int permutation[MATRIX_SIZE])
{
    int pivot = MATRIX_SIZE - 2;
    while (pivot >= 0 && permutation[pivot] >= permutation[pivot + 1]) {
        --pivot;
    }
    if (pivot < 0) {
        return 0;
    }

    int successor = MATRIX_SIZE - 1;
    while (permutation[successor] <= permutation[pivot]) {
        --successor;
    }
    int temporary = permutation[pivot];
    permutation[pivot] = permutation[successor];
    permutation[successor] = temporary;

    for (int left = pivot + 1, right = MATRIX_SIZE - 1;
         left < right; ++left, --right) {
        temporary = permutation[left];
        permutation[left] = permutation[right];
        permutation[right] = temporary;
    }
    return 1;
}

int main(void)
{
    double matrix[MATRIX_SIZE * MATRIX_SIZE] = {
         9.0,  3.0, -3.0, -9.0,
         3.0,  1.0, -1.0, -3.0,
        -3.0, -1.0,  1.0,  3.0,
        -9.0, -3.0,  3.0,  9.0
    };

#if DIVIDE_MATRIX_BY_FOUR
    for (int i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        matrix[i] /= 4.0;
    }
#endif

    double eigenvalues[MATRIX_SIZE];
    double eigenvectors[MATRIX_SIZE * MATRIX_SIZE];
    int info = 0;
    int found = 0;
    enum compute_status status = compute_eigensystem(
        matrix, eigenvalues, eigenvectors, &info, &found);
    if (status != COMPUTE_SUCCESS) {
        report_compute_error(status, info, found);
        return EXIT_FAILURE;
    }

    printf("Eigenvalues:\n");
    for (int i = 0; i < found; ++i) {
        printf("  %.12g\n", eigenvalues[i]);
    }

    printf("Eigenvectors (columns):\n");
    for (int row = 0; row < MATRIX_SIZE; ++row) {
        for (int column = 0; column < found; ++column) {
            printf("  % .12g", eigenvectors[row + column * MATRIX_SIZE]);
        }
        putchar('\n');
    }

    int permutation[MATRIX_SIZE] = {0, 1, 2, 3};
    int permutation_count = 0;
    do {
        double permuted_matrix[MATRIX_SIZE * MATRIX_SIZE];
        /* Apply each permutation to both rows and columns to keep symmetry. */
        for (int row = 0; row < MATRIX_SIZE; ++row) {
            for (int column = 0; column < MATRIX_SIZE; ++column) {
                permuted_matrix[row + column * MATRIX_SIZE] =
                    matrix[permutation[row] + permutation[column] * MATRIX_SIZE];
            }
        }

        status = compute_eigensystem(
            permuted_matrix, eigenvalues, eigenvectors, &info, &found);
        if (status != COMPUTE_SUCCESS) {
            report_compute_error(status, info, found);
            return EXIT_FAILURE;
        }
        ++permutation_count;
    } while (next_permutation(permutation));

    if (permutation_count != 24) {
        fprintf(stderr, "generated %d matrix permutations; expected 24\n",
                permutation_count);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
