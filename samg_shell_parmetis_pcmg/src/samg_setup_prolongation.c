#include "../include/main.h"

#define MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))

int SAMGLevelQMatrix(SAMGCtx **samg_ctx /*samg context data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    SAMGCtx *data_samg_ctx = *samg_ctx;

    int cnt_level = 0;
    for (cnt_level = 0; cnt_level < data_samg_ctx->num_level; ++cnt_level)
    {
        PetscCall(SAMGLevelKQMatrix(&data_samg_ctx->levels[cnt_level].data_f_mesh,
                                    &data_samg_ctx->levels[cnt_level].data_nullspace_levelk,
                                    &data_samg_ctx->levels[cnt_level].data_agg,
                                    &data_samg_ctx->levels[cnt_level].data_q_levelk));
    }

    return 0;
}

int SAMGLevelKQMatrix(MeshData *data_mesh_f /*fine-level mesh data*/,
                      NearNullSpaceLevelK *data_nullspace_levelk /*fine-level near null space data*/,
                      AggData *data_agg /*aggregation data*/,
                      QLevelK *data_q_levelk /*level k Q matrix*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    int *vtxdist_f = data_mesh_f->vtxdist;
    int size_global = data_mesh_f->nv;
    int vtx_start_f = vtxdist_f[my_rank];
    int size_local = vtxdist_f[my_rank + 1] - vtxdist_f[my_rank];

    data_q_levelk->size_global = size_global;
    data_q_levelk->size_local = size_local;
    data_q_levelk->data_q = (QDataVertexLevelK *)malloc(size_local * sizeof(QDataVertexLevelK));
    assert(data_q_levelk->data_q);

    for (int index = 0; index < size_local; ++index)
    {
        data_q_levelk->data_q[index].idx = data_nullspace_levelk->data_nullspace[index].idx;
        data_q_levelk->data_q[index].nrow = data_nullspace_levelk->data_nullspace[index].nrow;
        data_q_levelk->data_q[index].ncol = data_nullspace_levelk->data_nullspace[index].ncol;
        for (int index_q = 0; index_q < 36; ++index_q)
        {
            data_q_levelk->data_q[index].val[index_q] = 0.;
        }
    }

    const int BLOCK_CNT = 36;
    int *send_counts = (int *)calloc(nprocs, sizeof(int));
    int *recv_counts = (int *)calloc(nprocs, sizeof(int));
    int *sdispls = (int *)calloc(nprocs, sizeof(int));
    int *rdispls = (int *)calloc(nprocs, sizeof(int));
    assert(send_counts && recv_counts && sdispls && rdispls);

    // sending counts
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        if (data_agg->owner[index_p] == my_rank)
        {
            for (int index_r = 0; index_r < nprocs; ++index_r)
            {
                send_counts[index_r] += data_agg->nlocal_all[index_r * data_agg->np + index_p] * BLOCK_CNT;
            }
        }
    }

    // receving counts
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        int owner_rank = data_agg->owner[index_p];
        recv_counts[owner_rank] += data_agg->nlocal_all[my_rank * data_agg->np + index_p] * BLOCK_CNT;
    }

    // displacements
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        sdispls[index_r] = (index_r == 0) ? 0 : sdispls[index_r - 1] + send_counts[index_r - 1];
        rdispls[index_r] = (index_r == 0) ? 0 : rdispls[index_r - 1] + recv_counts[index_r - 1];
    }

    // buffer
    int total_send = sdispls[nprocs - 1] + send_counts[nprocs - 1];
    int total_recv = rdispls[nprocs - 1] + recv_counts[nprocs - 1];

    double *sendbuf = (double *)calloc((total_send > 0 ? total_send : 1), sizeof(double));
    double *recvbuf = (double *)calloc((total_recv > 0 ? total_recv : 1), sizeof(double));
    assert(sendbuf && recvbuf);

    // data packaging
    int *s_offset = (int *)malloc(nprocs * sizeof(int));
    assert(s_offset);
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        s_offset[index_r] = sdispls[index_r];
    }

    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        if (data_agg->owner[index_p] == my_rank)
        {
            GhostAggData *ghost = &data_agg->data_ghost_agg[index_p];
            NearNullSpaceDataVertexLevelK *ns_list = data_agg->fine_global_nullspace[index_p];

            int current_row_idx = 0;
            for (int index_v = 0; index_v < ghost->n_total; ++index_v)
            {
                int gid = ghost->flid2fgid[index_v];
                int target_rank = -1;

                for (int index_r = 0; index_r < nprocs; ++index_r)
                {
                    if (gid >= vtxdist_f[index_r] && gid < vtxdist_f[index_r + 1])
                    {
                        target_rank = index_r;
                        break;
                    }
                }

                int q_nrow = ns_list[index_v].nrow;
                int q_ncol = ns_list[index_v].ncol;
                int M = ghost->nrow; // nrow of mat_t or mat_q

                double *dst = &sendbuf[s_offset[target_rank]];
                int val_idx = 0;
                for (int index_q_c = 0; index_q_c < q_ncol; ++index_q_c)
                {
                    for (int index_q_r = 0; index_q_r < q_nrow; ++index_q_r)
                    {
                        int src_idx = (current_row_idx + index_q_r) + (index_q_c * M);
                        dst[val_idx] = ghost->mat_q[src_idx];
                        ++val_idx;
                    }
                }

                s_offset[target_rank] += BLOCK_CNT;
                current_row_idx += q_nrow;
            }
        }
    }

    // mpi-alltoallv
    (void)MPI_Alltoallv(sendbuf, send_counts, sdispls, MPI_DOUBLE,
                        recvbuf, recv_counts, rdispls, MPI_DOUBLE, comm);

    // data unpackaging
    int *r_offset = (int *)malloc(nprocs * sizeof(int));
    assert(r_offset);
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        r_offset[index_r] = rdispls[index_r];
    }

    for (int index = 0; index < size_local; ++index)
    {
        int gid = index + vtx_start_f;
        int part_id = data_agg->fgid2part[gid];
        int owner_rank = data_agg->owner[part_id];

        double *src = &recvbuf[r_offset[owner_rank]];
        memcpy(data_q_levelk->data_q[index].val, src, BLOCK_CNT * sizeof(double));

        r_offset[owner_rank] += BLOCK_CNT;
    }

#if 0
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, level k data_q, global size: %d, local size: %d\n", index_r, data_q_levelk->size_global, data_q_levelk->size_local);
            printf("local_vtx_id\t global_vtx_id\t nrow_q\t ncol_q\t data_q\n");
            for (int index_v = 0; index_v < size_local; ++index_v)
            {
                int nrow_q = data_q_levelk->data_q[index_v].nrow;
                int ncol_q = data_q_levelk->data_q[index_v].ncol;
                printf("%d\t %d\t %d\t %d\n", index_v, index_v + vtx_start_f, nrow_q, ncol_q);
                for (int index_q_r = 0; index_q_r < nrow_q; ++index_q_r)
                {
                    printf(" \t \t \t \t ");
                    for (int index_q_c = 0; index_q_c < ncol_q; ++index_q_c)
                    {
                        printf("%021.16le\t ", data_q_levelk->data_q[index_v].val[MAT_COL_MAJOR(index_q_r, index_q_c, nrow_q)]);
                    }
                    printf("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check level k data_q

    // free memory
    free(r_offset);
    free(s_offset);
    free(sendbuf);
    free(recvbuf);
    free(send_counts);
    free(recv_counts);
    free(sdispls);
    free(rdispls);

    return 0;
}

void QRFactorizationOpenBLAS(int m /*nrow*/, int n /*ncol*/,
                             double *Q /*m x n*/, double *R /*n x n*/)
{
    double *tau = (double *)malloc(n * sizeof(double));
    assert(tau);

    int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, n, Q, m, tau);
    assert(info == 0);

    // store R, upper triangular part of Q
    for (int index_j = 0; index_j < n; ++index_j)
    {
        for (int index_i = 0; index_i <= index_j; ++index_i)
        {
            R[MAT_COL_MAJOR(index_i, index_j, n)] = Q[MAT_COL_MAJOR(index_i, index_j, m)];
        }
    }

    // store Q
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, n, n, Q, m, tau);
    assert(info == 0);

    // free memory
    free(tau);
}

int SAMGTentativeProlongationOperator(SAMGCtx **samg_ctx /*samg context data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    SAMGCtx *data_samg_ctx = *samg_ctx;

    int cnt_level = 0;
    int num_level = data_samg_ctx->num_level;

    for (cnt_level = 0; cnt_level < num_level; ++cnt_level)
    {
        PetscCall(SAMGLevelKTentativeProlongationOperator(cnt_level, &data_samg_ctx->levels[cnt_level]));

        // calling PtAP
        PetscCall(MatPtAP(data_samg_ctx->levels[cnt_level].op_f,
                          data_samg_ctx->levels[cnt_level].op_ua_p,
                          MAT_INITIAL_MATRIX,
                          PETSC_DETERMINE,
                          &(data_samg_ctx->levels[cnt_level].op_c)));

        // next level
        PetscCall(MatDuplicate(data_samg_ctx->levels[cnt_level].op_c,
                               MAT_COPY_VALUES,
                               &(data_samg_ctx->levels[cnt_level + 1].op_f)));
    }

    return 0;
}

int SAMGLevelKTentativeProlongationOperator(int level /*current level*/,
                                            MGLevel *data_level /*level data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    // size of op_f
    int nrow_op_f = 0, ncol_op_f = 0;
    int local_row_op_f = 0, local_col_op_f = 0;
    PetscCall(MatGetSize(data_level->op_f, &nrow_op_f, &ncol_op_f));
    PetscCall(MatGetLocalSize(data_level->op_f, &local_row_op_f, &local_col_op_f));
#if 0
    PetscCall(PetscPrintf(comm, "in level %d, size of op_f: (%d, %d)\n", level, nrow_op_f, ncol_op_f));
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, local size of op_f: (%d, %d)\n", index_r, local_row_op_f, local_col_op_f);
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check fine-level operator size

    int *vtxdist_c = data_level->data_c_mesh.vtxdist;
    int local_nv_c = vtxdist_c[my_rank + 1] - vtxdist_c[my_rank];
    int global_nv_c = data_level->data_c_mesh.nv;
    int nrow_op_ua_p = nrow_op_f, ncol_op_ua_p = global_nv_c * 6;
    int local_row_op_ua_p = local_row_op_f, local_col_op_ua_p = local_nv_c * 6;
    PetscCall(MatCreate(comm, &data_level->op_ua_p));
    PetscCall(MatSetSizes(data_level->op_ua_p,
                          local_row_op_ua_p,
                          local_col_op_ua_p,
                          nrow_op_ua_p,
                          ncol_op_ua_p));
    PetscCall(MatSetType(data_level->op_ua_p, MATAIJ));
    PetscCall(MatSetUp(data_level->op_ua_p));

    int *vtxdist_f = data_level->data_f_mesh.vtxdist;
    int vtx_f_start = vtxdist_f[my_rank];
    int vtx_f_end = vtxdist_f[my_rank + 1];
    int local_nv_f = vtx_f_end - vtx_f_start;

    // local vertex id mapping to block start row and col of P
    int *flid2prow = (int *)calloc(local_nv_f, sizeof(int));
    int *flid2pcol = (int *)calloc(local_nv_f, sizeof(int));
    assert(flid2prow && flid2pcol);

    for (int index = 0; index < local_nv_f; ++index)
    {
        int gid = data_level->data_q_levelk.data_q[index].idx;
        int part_id = data_level->data_agg.fgid2part[gid];

        flid2pcol[index] = 6 * part_id;
    }

    for (int index = 0; index < local_nv_f - 1; ++index)
    {
        int q_nrow = data_level->data_q_levelk.data_q[index].nrow;
        flid2prow[index + 1] = flid2prow[index] + q_nrow;
    }

    int rstart_global;
    PetscCall(MatGetOwnershipRange(data_level->op_ua_p, &rstart_global, NULL));
    for (int index = 0; index < local_nv_f; ++index)
    {

        int index_row_start = flid2prow[index];
        int index_col_start = flid2pcol[index];
        int q_nrow = data_level->data_q_levelk.data_q[index].nrow;
        int q_ncol = data_level->data_q_levelk.data_q[index].ncol;
        for (int index_q_r = 0; index_q_r < q_nrow; ++index_q_r)
        {
            for (int index_q_c = 0; index_q_c < q_ncol; ++index_q_c)
            {
                PetscCall(MatSetValue(data_level->op_ua_p,
                                      rstart_global + index_row_start + index_q_r,
                                      index_col_start + index_q_c,
                                      data_level->data_q_levelk.data_q[index].val[MAT_COL_MAJOR(index_q_r, index_q_c, q_nrow)],
                                      INSERT_VALUES));
            }
        }
    }
    PetscCall(MatAssemblyBegin(data_level->op_ua_p, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(data_level->op_ua_p, MAT_FINAL_ASSEMBLY));
#if 0
    int tmp_nrow_op_ua_p = 0, tmp_ncol_op_ua_p = 0;
    int tmp_local_row_op_ua_p = 0, tmp_local_col_op_ua_p = 0;
    PetscCall(MatGetSize(data_level->op_ua_p, &tmp_nrow_op_ua_p, &tmp_ncol_op_ua_p));
    PetscCall(MatGetLocalSize(data_level->op_ua_p, &tmp_local_row_op_ua_p, &tmp_local_col_op_ua_p));

    PetscCall(PetscPrintf(comm, "in level %d, size of op_ua_p: (%d, %d)\n", level, tmp_nrow_op_ua_p, tmp_ncol_op_ua_p));
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, local size of op_ua_p: (%d, %d)\n", index_r, tmp_local_row_op_ua_p, tmp_local_col_op_ua_p);
            puts("\n");
        }
        fflush(stdout);
    }
    PetscCall(MatView(data_level->op_ua_p, PETSC_VIEWER_STDOUT_WORLD));
#endif // check tentative prolongation operator size

    // free memory
    free(flid2prow);
    free(flid2pcol);

    return 0;
}

int SAMGApplyProlongationSmoother(int n /*column size of prolongation operator*/,
                                  double omega /*scaling weight parameter*/,
                                  PSmoother *p_s /*prolongation operator smoother*/,
                                  Mat *p_sa /*smoothed prolongation operator*/,
                                  Mat *p_ua /*unsmoothed prolongation operator*/,
                                  Mat *A /*fine-level operator*/)
{
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;

    Vec c_vec_p_ua, c_vec_p_sa; // column vector of prolongation operator
    Vec r_tmp;                  // r_tmp = A * c_vec_p_ua
    PetscScalar *c_vec_array;
    PetscInt *row_indices;
    PetscInt m, M, N;

    // Create vectors compatible with matrix A (older PETSc method)
    PetscCall(MatGetSize(*A, &M, &N));
    PetscCall(MatGetLocalSize(*A, &m, NULL));

    // Create first vector
    PetscCall(VecCreate(comm, &c_vec_p_ua));
    PetscCall(VecSetSizes(c_vec_p_ua, m, M));
    PetscCall(VecSetFromOptions(c_vec_p_ua));

    // Duplicate for other vectors
    PetscCall(VecDuplicate(c_vec_p_ua, &r_tmp));
    PetscCall(VecDuplicate(r_tmp, &c_vec_p_sa));

    // Get matrix dimensions and create row indices array
    PetscCall(PetscMalloc1(m, &row_indices));
    PetscCall(MatGetOwnershipRange(*A, &row_indices[0], &row_indices[m - 1]));
    for (int i = 0; i < m; ++i)
    {
        row_indices[i] = row_indices[0] + i;
    }

    for (int index = 0; index < n; ++index)
    {
        PetscCall(MatGetColumnVector(*p_ua, c_vec_p_ua, index));
        PetscCall(MatMult(*A, c_vec_p_ua, r_tmp));
        PetscCall(PCApply(p_s->pc, r_tmp, c_vec_p_sa));
        PetscCall(VecAXPY(c_vec_p_sa, -omega, c_vec_p_ua));

        // Insert c_vec_p_sa as column 'index' into p_sa
        PetscCall(VecGetArrayRead(c_vec_p_sa, &c_vec_array));
        PetscCall(MatSetValues(*p_sa, m, row_indices, 1, &index, c_vec_array, INSERT_VALUES));
        PetscCall(VecRestoreArrayRead(c_vec_p_sa, &c_vec_array));
    }

    // Assemble the matrix after all insertions
    PetscCall(MatAssemblyBegin(*p_sa, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*p_sa, MAT_FINAL_ASSEMBLY));

    PetscCall(PetscFree(row_indices));
    return 0;
}

int SAMGSmoothedProlongation(MGLevel *level /*level hierarchy data*/)
{
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;

    int ps_num_steps = level->op_s.num_steps;
    int smoother_type = level->op_s.smoother_type;
    double omega = level->op_s.smoother_scale;
    int m_p_ua, n_p_ua;
    int local_m_p_ua, local_n_p_ua;

    KSP ksp;
    PetscCall(KSPCreate(comm, &ksp));
    PetscCall(KSPSetOperators(ksp, level->op_f, level->op_f));
    PetscCall(KSPGetPC(ksp, &level->op_s.pc));

    switch (smoother_type)
    {
    case 0:
        // jacobi
        PetscCall(PCSetType(level->op_s.pc, PCJACOBI));
        break;
    case 1:
        // sor
        PetscCall(PCSetType(level->op_s.pc, PCSOR));
        break;
    default:
        break;
    }

    // unsmoothed prolongation size
    PetscCall(MatGetSize(level->op_ua_p, &m_p_ua, &n_p_ua));
    PetscCall(MatGetLocalSize(level->op_ua_p, &local_m_p_ua, &local_n_p_ua));

    // create smoothed prolongation size
    PetscCall(MatCreate(comm, &level->op_sa_p));
    PetscCall(MatSetSizes(level->op_sa_p,
                          local_m_p_ua,
                          local_n_p_ua,
                          m_p_ua,
                          n_p_ua));
    PetscCall(MatSetType(level->op_sa_p, MATAIJ));
    PetscCall(MatSetUp(level->op_sa_p));

    Mat tmp;
    PetscCall(MatDuplicate(level->op_ua_p, MAT_COPY_VALUES, &tmp));

    for (int index = 0; index < ps_num_steps; ++index)
    {
        PetscCall(SAMGApplyProlongationSmoother(n_p_ua,
                                                omega,
                                                &level->op_s,
                                                &level->op_sa_p,
                                                &tmp,
                                                &level->op_f));

        PetscCall(MatCopy(level->op_sa_p, tmp, DIFFERENT_NONZERO_PATTERN));
    }

    // free
    PetscCall(MatDestroy(&tmp));
    PetscCall(KSPDestroy(&ksp));

    return 0;
}
