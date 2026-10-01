#include <stdio.h>
#include <stdlib.h>

extern void dsyevr_(char *jobz, char *range, char *uplo, int *n,
                    double *a, int *lda, double *vl, double *vu, int *il,
                    int *iu, double *abstol, int *m, double *w, double *z,
                    int *ldz, int *isuppz, double *work, int *lwork,
                    int *iwork, int *liwork, int *info);

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

int main(void)
{
    int n = 4;
    int lda = 4;
    int ldz = 4;
    /* LAPACK stores matrices column-major. */
    double matrix[] = {
         9.0,  3.0, -3.0, -9.0,
         3.0,  1.0, -1.0, -3.0,
        -3.0, -1.0,  1.0,  3.0,
        -9.0, -3.0,  3.0,  9.0
    };
    double eigenvalues[4];
    double eigenvectors[16];
    int isuppz[8];
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    double vl = 0.0;
    double vu = 0.0;
    int il = 0;
    int iu = 0;
    double abstol = 0.0;
    int found = 0;
    int info = 0;
    double work_query = 0.0;
    int iwork_query = 0;
    int lwork = -1;
    int liwork = -1;

    /* Query workspace sizes first; -1 asks for sizes without computing. */
    dsyevr_(&jobz, &range, &uplo, &n, matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, &found, eigenvalues, eigenvectors, &ldz, isuppz,
            &work_query, &lwork, &iwork_query, &liwork, &info);
    if (info != 0) {
        report_dsyevr_error("workspace query", info);
        return EXIT_FAILURE;
    }

    lwork = (int)work_query;
    liwork = iwork_query;
    double *work = malloc((size_t)lwork * sizeof(*work));
    int *iwork = malloc((size_t)liwork * sizeof(*iwork));
    if (work == NULL || iwork == NULL) {
        fprintf(stderr, "could not allocate LAPACK workspace\n");
        free(work);
        free(iwork);
        return EXIT_FAILURE;
    }

    /* Call again with the allocated workspaces to compute the eigensystem. */
    dsyevr_(&jobz, &range, &uplo, &n, matrix, &lda, &vl, &vu, &il, &iu,
            &abstol, &found, eigenvalues, eigenvectors, &ldz, isuppz,
            work, &lwork, iwork, &liwork, &info);
    free(work);
    free(iwork);
    if (info != 0) {
        report_dsyevr_error("eigensolve", info);
        return EXIT_FAILURE;
    }
    if (found != n) {
        fprintf(stderr, "dsyevr returned %d eigenvalues; expected %d\n",
                found, n);
        return EXIT_FAILURE;
    }

    printf("Eigenvalues:\n");
    for (int i = 0; i < found; ++i) {
        printf("  %.12g\n", eigenvalues[i]);
    }

    printf("Eigenvectors (columns):\n");
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < found; ++column) {
            printf("  % .12g", eigenvectors[row + column * ldz]);
        }
        putchar('\n');
    }

    return EXIT_SUCCESS;
}
