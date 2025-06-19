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

    vec_data->val = (double *)malloc(vec_data->nrows * sizeof(double));
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

    int local_node_start = node_vtxdist[my_rank],
        local_node_end = node_vtxdist[my_rank + 1];

    int local_index_start = local_node_start * 6,
        local_index_end = local_node_end * 6;

    // matrix file process
    CSRMatrix mat_data;
    mat_data.local_nrows = (local_node_end - local_node_start) * 6; // 6 dofs each node

    mat_data.local_row_ptr = (int *)malloc((mat_data.local_nrows + 1) * sizeof(int));
    assert(mat_data.local_row_ptr);

    if (my_rank == 0)
    {
        FileProcessCSRMatrix(path_mat, &mat_data);

        memcpy(mat_data.local_row_ptr,
               mat_data.row_ptr,
               (mat_data.local_nrows + 1) * sizeof(int));

        int tmp_csr_shift = mat_data.local_row_ptr[0];
        for (int index = 0; index < mat_data.local_nrows + 1; ++index)
        {
            mat_data.local_row_ptr[index] -= tmp_csr_shift; // local row_ptr start from 0
        }

        mat_data.local_nnz = mat_data.local_row_ptr[mat_data.local_nrows] - mat_data.local_row_ptr[0];

        mat_data.local_col_idx = (int *)malloc(mat_data.local_nnz * sizeof(int));
        mat_data.local_val = (double *)malloc(mat_data.local_nnz * sizeof(double));
        assert(mat_data.local_col_idx && mat_data.local_val);

        memcpy(mat_data.local_col_idx,
               mat_data.col_idx,
               mat_data.local_nnz * sizeof(int));

        memcpy(mat_data.local_val,
               mat_data.val,
               mat_data.local_nnz * sizeof(double));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            int tmp_local_node_start = node_vtxdist[index_p],
                tmp_local_node_end = node_vtxdist[index_p + 1];

            int tmp_local_nrows = (tmp_local_node_end - tmp_local_node_start) * 6;

            int tmp_local_index_start = tmp_local_node_start * 6;

            (void)MPI_Send(mat_data.row_ptr + tmp_local_index_start, tmp_local_nrows + 1, MPI_INT,
                           index_p, 0, comm);

            int tmp_local_nnz = mat_data.row_ptr[tmp_local_index_start + tmp_local_nrows] -
                                mat_data.row_ptr[tmp_local_index_start];

            (void)MPI_Send(mat_data.col_idx + mat_data.row_ptr[tmp_local_index_start],
                           tmp_local_nnz, MPI_INT, index_p, 1, comm); // local_col_idx

            (void)MPI_Send(mat_data.val + mat_data.row_ptr[tmp_local_index_start],
                           tmp_local_nnz, MPI_DOUBLE, index_p, 2, comm); // local_val
        }
    }
    else
    {
        (void)MPI_Recv(mat_data.local_row_ptr, mat_data.local_nrows + 1, MPI_INT,
                       0, 0, comm, MPI_STATUS_IGNORE);

        int tmp_csr_shift = mat_data.local_row_ptr[0];
        for (int index = 0; index < mat_data.local_nrows + 1; ++index)
        {
            mat_data.local_row_ptr[index] -= tmp_csr_shift;
        }

        mat_data.local_nnz = mat_data.local_row_ptr[mat_data.local_nrows];

        mat_data.local_col_idx = (int *)malloc(mat_data.local_nnz * sizeof(int));
        mat_data.local_val = (double *)malloc(mat_data.local_nnz * sizeof(double));
        assert(mat_data.local_col_idx && mat_data.local_val);

        (void)MPI_Recv(mat_data.local_col_idx, mat_data.local_nnz, MPI_INT,
                       0, 1, comm, MPI_STATUS_IGNORE);

        (void)MPI_Recv(mat_data.local_val, mat_data.local_nnz, MPI_DOUBLE,
                       0, 2, comm, MPI_STATUS_IGNORE);
    }

    (void)MPI_Bcast(&mat_data.nrows, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&mat_data.ncols, 1, MPI_INT, 0, comm);

    // vector file process
    CSRVector rhs_data;
    rhs_data.local_nrows = (local_node_end - local_node_start) * 6; // 6 dofs each node

    rhs_data.local_val = (double *)malloc(rhs_data.local_nrows * sizeof(double));
    assert(rhs_data.local_val);

    double *global_rhs_val = NULL;
    int *rhs_sendcounts = NULL, *rhs_displs = NULL;

    if (my_rank == 0)
    {
        FileProcessCSRVector(path_rhs, &rhs_data);

        global_rhs_val = rhs_data.val;

        rhs_sendcounts = (int *)malloc(nprocs * sizeof(int));
        rhs_displs = (int *)malloc(nprocs * sizeof(int));

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int tmp_local_node_start = node_vtxdist[index_p],
                tmp_local_node_end = node_vtxdist[index_p + 1];
            rhs_sendcounts[index_p] = (tmp_local_node_end - tmp_local_node_start) * 6;
            rhs_displs[index_p] = tmp_local_node_start * 6;
        }
#if 0
        int local_node_start = node_vtxdist[my_rank],
            local_node_end = node_vtxdist[my_rank + 1];

        int local_nrows = (local_node_end - local_node_start) * 6;

        int local_rhs_row_start = local_node_start * 6;

        memcpy(rhs_data.local_val,
               rhs_data.val + local_rhs_row_start,
               rhs_data.local_nrows * sizeof(double));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            local_node_start = node_vtxdist[index_p];
            local_node_end = node_vtxdist[index_p + 1];
            local_nrows = (local_node_end - local_node_start) * 6;
            local_rhs_row_start = local_node_start * 6;
            (void)MPI_Send(rhs_data.val + local_rhs_row_start, local_nrows, MPI_DOUBLE,
                           index_p, 10, comm);
        }
#endif // mpi_send()
    }
#if 0
    else
    {
        (void)MPI_Recv(rhs_data.local_val, rhs_data.local_nrows, MPI_DOUBLE,
                       0, 10, comm, MPI_STATUS_IGNORE);
    }
#endif // mpi_recv()

    (void)MPI_Bcast(&rhs_data.nrows, 1, MPI_INT, 0, comm);
    (void)MPI_Scatterv(global_rhs_val, rhs_sendcounts, rhs_displs, MPI_DOUBLE,
                       rhs_data.local_val, rhs_data.local_nrows, MPI_DOUBLE, 0, comm);

    // distributed matrix
    PetscCall(MatCreate(comm, &mysolver->solver_a));
    PetscCall(MatSetSizes(mysolver->solver_a,
                          mat_data.local_nrows, rhs_data.local_nrows,
                          mat_data.nrows, mat_data.ncols));
    PetscCall(MatSetType(mysolver->solver_a, MATAIJ));
    PetscCall(MatSetUp(mysolver->solver_a));

    for (int index = local_index_start; index < local_index_end; ++index)
    {
        int index_start = mat_data.local_row_ptr[index - local_index_start],
            index_end = mat_data.local_row_ptr[index - local_index_start + 1];

        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscCall(MatSetValue(mysolver->solver_a,
                                  index,
                                  mat_data.local_col_idx[index_j],
                                  mat_data.local_val[index_j],
                                  INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(mysolver->solver_a, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(mysolver->solver_a, MAT_FINAL_ASSEMBLY));

    // distributed vector
    PetscCall(VecCreate(comm, &mysolver->solver_b));
    PetscCall(VecSetSizes(mysolver->solver_b, rhs_data.local_nrows, rhs_data.nrows));
    PetscCall(VecSetType(mysolver->solver_b, VECMPI));

    for (int index = local_index_start; index < local_index_end; ++index)
    {
        PetscScalar val_tmp = rhs_data.local_val[index - local_index_start];

        PetscCall(VecSetValues(mysolver->solver_b,
                               1,
                               &index,
                               &val_tmp,
                               INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(mysolver->solver_b));
    PetscCall(VecAssemblyEnd(mysolver->solver_b));

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
    free(mat_data.local_val);
    free(mat_data.local_col_idx);
    free(mat_data.local_row_ptr);
    free(rhs_data.local_val);
    if (my_rank == 0)
    {
        // csr matrix data
        free(mat_data.row_ptr);
        free(mat_data.col_idx);
        free(mat_data.val);

        // csr rhs data
        free(rhs_data.val);
        free(rhs_sendcounts);
        free(rhs_displs);
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
