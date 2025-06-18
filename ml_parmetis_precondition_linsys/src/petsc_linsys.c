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

static int FileProcessCSRMatrix(const char *path, CSRMatrix *mat_data)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    fscanf(fp, " %d %d ", &mat_data->nrows, &mat_data->ncols);
    fscanf(fp, "%d", &mat_data->nnz);

    mat_data->row_ptr = (int *)malloc((mat_data->nrows + 1) * sizeof(int));
    mat_data->col_idx = (int *)malloc(mat_data->nnz * sizeof(int));
    mat_data->val = (double *)malloc(mat_data->nnz * sizeof(double));
    assert(mat_data->row_ptr && mat_data->col_idx && mat_data->val);

    for (int index = 0; index < mat_data->nrows + 1; ++index)
    {
        fscanf(fp, "%d", mat_data->row_ptr + index);
    }

    for (int index = 0; index < mat_data->nnz; ++index)
    {
        fscanf(fp, "%d", mat_data->col_idx + index);
    }

    for (int index = 0; index < mat_data->nnz; ++index)
    {
        fscanf(fp, "%lf", mat_data->val + index);
    }

    fclose(fp);

    return 0;
}

static int FileProcessCSRVector(const char *path, CSRVector *vec_data)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    fscanf(fp, "%d", &vec_data->nrows);

    vec_data->val = (double *)malloc(vec_data->nrow * sizeof(double));
    assert(vec_data->val);

    for (int index = 0; index < vec_data->nrows; ++index)
    {
        fscanf(fp, "%lf", vec_data->val + index);
    }

    fclose(fp);

    return 0;
}

int SolverPetscInitialize(const char *path_mat, const char *path_rhs, const int *node_vtxdist,
                          MySolver *mysolver)
{
    int my_rank, nprocs;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    int local_index_start = node_vtxdist[my_rank],
        local_index_end = node_vtxdist[my_rank + 1];

    // matrix file process
    CSRMatrix mat_data;
    mat_data.local_nrows = (local_index_end - local_index_start) * 6; // 6 dofs each node

    if (my_rank == 0)
    {
        FileProcessCSRMatrix(path_mat, &mat_data);
    }

    // vector file process
    CSRVector rhs_data;
    mat_data.local_nrows = (local_index_end - local_index_start) * 6; // 6 dofs each node

    if (my_rank == 0)
    {
        FileProcessCSRVector(path_rhs, &rhs_data);
    }

#if 0
    PetscViewer fd;

    // matrix file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_mat, FILE_MODE_READ, &fd));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &(mysolver->solver_a)));
    PetscCall(MatLoad(mysolver->solver_a, fd));
    PetscCall(PetscViewerDestroy(&fd));

    // rhs file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_rhs, FILE_MODE_READ, &fd));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(mysolver->solver_b)));
    PetscCall(VecLoad(mysolver->solver_b, fd));
    PetscCall(PetscViewerDestroy(&fd));
#endif // petsc binary matrix

    // sol vector and residual vector
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_x)));
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_r)));

    // size of matrix
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==== basic information of linear system ====\n"));
    PetscInt m_mat = 0, n_mat = 0; // nnz_mat = 0;
    PetscCall(MatGetSize(mysolver->solver_a, &m_mat, &n_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix Row = %" PetscInt_FMT ", matrix Column = %" PetscInt_FMT "\n", m_mat, n_mat));
    MatInfo info_mat;
    PetscCall(MatGetInfo(mysolver->solver_a, MAT_GLOBAL_SUM, &info_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix nz_allocated = %ld, matrix nz_used = %ld, matrix nz_unneeded = %ld\n",
                          (long)(info_mat.nz_allocated), (long)(info_mat.nz_used), (long)(info_mat.nz_unneeded)));
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
    if (my_rank == 0)
    {
        // csr matrix data
        free(mat_data.row_ptr);
        free(mat_data.col_idx);
        free(mat_data.val);

        // csr rhs data
        free(rhs_data.val);
    }

    return 0;
}

#if 0
void SolverPetscPreprocess(int argc, char **argv, MySolver *mysolver)
{
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));
    PetscCall(KSPSetFromOptions(mysolver->ksp));
}

void SolverPetscSolve(int argc, char **argv, MySolver *mysolver)
{
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
}
#endif // uncomment
