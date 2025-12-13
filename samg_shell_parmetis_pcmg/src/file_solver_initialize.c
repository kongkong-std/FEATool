#include "../include/main.h"

int SolverPetscResidualCheck(MySolver *mysolver)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "            || b ||_2 = %021.16le\n", b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "|| r ||_2 / || b ||_2 = %021.16le\n", r_norm_2 / b_norm_2));

    return 0;
}

int FileProcessCSRMatrix(const char *path_mat /*path to matrix file*/,
                         CSRMatrix *data_mat)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    char path[PETSC_MAX_PATH_LEN];
    sprintf(path, "%s/A-%d.txt", path_mat, my_rank);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        if (index == my_rank)
        {
            printf(">>>> in rank %d:\n", my_rank);
            printf("path to matrix file: %s\n", path);
        }
    }
#endif // test path

    FILE *fp = fopen(path, "rb");
    assert(fp);

    int nrows = 0, ncols = 0, nnz = 0;
    fscanf(fp, " %d %d %d ", &nrows, &ncols, &nnz);
    data_mat->nrows = nrows;
    data_mat->ncols = ncols;
    data_mat->nnz = nnz;

    data_mat->row_idx = (int *)malloc(nrows * sizeof(int));
    data_mat->row_ptr = (int *)malloc((nrows + 1) * sizeof(int));
    data_mat->col_idx = (int *)malloc(nnz * sizeof(int));
    data_mat->val = (double *)malloc(nnz * sizeof(double));
    assert(data_mat->row_idx &&
           data_mat->row_ptr &&
           data_mat->col_idx &&
           data_mat->val);

    for (int index = 0; index < nrows; ++index)
    {
        fscanf(fp, " %d ", data_mat->row_idx + index);
    }
    for (int index = 0; index < nrows + 1; ++index)
    {
        fscanf(fp, " %d ", data_mat->row_ptr + index);
    }
    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, " %d ", data_mat->col_idx + index);
    }
    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, " %lf ", data_mat->val + index);
    }
#if 0
    for (int rank = 0; rank < nprocs; ++rank)
    {
        if (rank == my_rank)
        {
            printf(">>>> CSR matrix data in rank %d:\n", rank);
            for (int index = 0; index < nrows; ++index)
            {
                printf("data_mat->row_idx[%d] = %d\n", index, data_mat->row_idx[index]);
            }
            for (int index = 0; index < nrows + 1; ++index)
            {
                printf("data_mat->row_ptr[%d] = %d\n", index, data_mat->row_ptr[index]);
            }
            for (int index = 0; index < nnz; ++index)
            {
                printf("data_mat->col_idx[%d] = %d\n", index, data_mat->col_idx[index]);
            }
            for (int index = 0; index < nnz; ++index)
            {
                printf("data_mat->val[%d] = %021.16le\n", index, data_mat->val[index]);
            }
        }
    }
#endif // debug csr data

    fclose(fp);

    return 0;
}

int FileProcessCSRVector(const char *path_rhs /*path to rhs file*/,
                         CSRVector *data_rhs)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    char path[PETSC_MAX_PATH_LEN];
    sprintf(path, "%s/b-%d.txt", path_rhs, my_rank);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        if (index == my_rank)
        {
            printf(">>>> in rank %d:\n", my_rank);
            printf("path to rhs file: %s\n", path);
        }
    }
#endif // test path

    FILE *fp = fopen(path, "rb");
    assert(fp);

    int nrows = 0;
    fscanf(fp, " %d ", &nrows);
    data_rhs->nrows = nrows;

    data_rhs->row_idx = (int *)malloc(nrows * sizeof(int));
    data_rhs->val = (double *)malloc(nrows * sizeof(double));
    assert(data_rhs->row_idx && data_rhs->val);

    for (int index = 0; index < nrows; ++index)
    {
        fscanf(fp, " %d %lf ", data_rhs->row_idx + index,
               data_rhs->val + index);
    }

    fclose(fp);

    return 0;
}

int SolverInitializeWithFile(const char *path_mat /*path to matrix file*/,
                             const char *path_rhs /*path to rhs file*/,
                             MySolver *mysolver /*solver data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    CSRMatrix data_mat;
    CSRVector data_rhs;

    FileProcessCSRMatrix(path_mat, &data_mat);
    FileProcessCSRVector(path_rhs, &data_rhs);

#if 1
    // petsc matrix
    PetscCall(MatCreate(comm, &mysolver->solver_a));
    PetscCall(MatSetSizes(mysolver->solver_a,
                          data_mat.nrows,
                          data_mat.ncols,
                          PETSC_DETERMINE,
                          PETSC_DETERMINE));
    PetscCall(MatSetType(mysolver->solver_a, MATAIJ));
    PetscCall(MatSetUp(mysolver->solver_a));

    for (int index = 0; index < data_mat.nrows; ++index)
    {
        int index_start = data_mat.row_ptr[index];
        int index_end = data_mat.row_ptr[index + 1];

        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscCall(MatSetValue(mysolver->solver_a,
                                  data_mat.row_idx[index],
                                  data_mat.col_idx[index_j],
                                  data_mat.val[index_j],
                                  INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(mysolver->solver_a, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(mysolver->solver_a, MAT_FINAL_ASSEMBLY));

    // size of matrix
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==== basic information of linear system ====\n"));
    PetscInt m_mat = 0, n_mat = 0; // nnz_mat = 0;
    PetscCall(MatGetSize(mysolver->solver_a, &m_mat, &n_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix Row = %" PetscInt_FMT ", matrix Column = %" PetscInt_FMT "\n", m_mat, n_mat));
    MatInfo info_mat;
    PetscCall(MatGetInfo(mysolver->solver_a, MAT_GLOBAL_SUM, &info_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix nz_allocated = %ld, matrix nz_used = %ld, matrix nz_unneeded = %ld\n",
                          (long)(info_mat.nz_allocated), (long)(info_mat.nz_used), (long)(info_mat.nz_unneeded)));
#endif // petsc matrix assembling

#if 1
    // petsc vector
    PetscCall(VecCreate(comm, &mysolver->solver_b));
    PetscCall(VecSetSizes(mysolver->solver_b, data_rhs.nrows, PETSC_DETERMINE));
    PetscCall(VecSetType(mysolver->solver_b, VECMPI));

    for (int index = 0; index < data_rhs.nrows; ++index)
    {
        PetscScalar val_tmp = data_rhs.val[index];

        PetscCall(VecSetValues(mysolver->solver_b,
                               1,
                               data_rhs.row_idx + index,
                               &val_tmp,
                               INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(mysolver->solver_b));
    PetscCall(VecAssemblyEnd(mysolver->solver_b));

    // sol vector and residual vector
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_x)));
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_r)));
#endif // petsc vector assembling

    PetscInt n_vec = 0;
    PetscCall(VecGetSize(mysolver->solver_b, &n_vec));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "vector Row = %" PetscInt_FMT "\n", n_vec));

    PetscReal b_norm_2 = 0.;
    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "initial     || b ||_2 = %021.16le\n", b_norm_2));

#if 1
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mysolver->ksp)));
#endif

    // free memory
    free(data_mat.row_idx);
    free(data_mat.row_ptr);
    free(data_mat.col_idx);
    free(data_mat.val);
    free(data_rhs.row_idx);
    free(data_rhs.val);

    return 0;
}
