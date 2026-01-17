#include "../include/main.h"

#define MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))

int SAMGLevelNearNullSpace(SAMGCtx **samg_ctx /*samg context data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    SAMGCtx *data_samg_ctx = *samg_ctx;

    int cnt_level = 0;
    PetscCall(SAMGLevel0NearNullSpace(&data_samg_ctx->data_nullspace_level0,
                                      &data_samg_ctx->levels[cnt_level].data_nullspace_levelk));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, level 0 fine near null space data, global nv = %d, local nv = %d\n", index,
                   data_samg_ctx->levels[cnt_level].data_nullspace_levelk.size_global,
                   data_samg_ctx->levels[cnt_level].data_nullspace_levelk.size_local);
            printf("local_vtx_id\t global_vtx_id\t near_null_space\n");
            for (int index_v = 0; index_v < data_samg_ctx->levels[cnt_level].data_nullspace_levelk.size_local; ++index_v)
            {
                printf("%d\t %d\t ", index_v,
                       data_samg_ctx->levels[cnt_level].data_nullspace_levelk.data_nullspace[index_v].idx);
                for (int index_ns_r = 0; index_ns_r < data_samg_ctx->levels[cnt_level].data_nullspace_levelk.data_nullspace[index_v].nrow; ++index_ns_r)
                {
                    printf(" \t \t ");
                    for (int index_ns_c = 0; index_ns_c < data_samg_ctx->levels[cnt_level].data_nullspace_levelk.data_nullspace[index_v].ncol; ++index_ns_c)
                    {
                        printf("%021.16le\t ", data_samg_ctx->levels[cnt_level].data_nullspace_levelk.data_nullspace[index_v].val[MAT_COL_MAJOR(index_ns_r, index_ns_c, data_samg_ctx->levels[cnt_level].data_nullspace_levelk.data_nullspace[index_v].nrow)]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check level 0 fine near null space data

    for (cnt_level = 0; cnt_level < data_samg_ctx->num_level; ++cnt_level)
    {
        PetscCall(SAMGLevelKNearNullSpace(&data_samg_ctx->levels[cnt_level].data_f_mesh,
                                          &data_samg_ctx->levels[cnt_level].data_c_mesh,
                                          &data_samg_ctx->levels[cnt_level].data_nullspace_levelk,
                                          &data_samg_ctx->levels[cnt_level].data_agg,
                                          &data_samg_ctx->levels[cnt_level + 1].data_nullspace_levelk));
    }

    return 0;
}

int SAMGLevelKNearNullSpace(MeshData *data_mesh_f /*fine-level mesh data*/,
                            MeshData *data_mesh_c /*coarse-level mesh data*/,
                            NearNullSpaceLevelK *data_nullspace_f /*fine-level near null space*/,
                            AggData *data_agg /*aggregation data*/,
                            NearNullSpaceLevelK *data_nullspace_c /*coarse-level near null space*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    int *vtxdist_f = data_mesh_f->vtxdist;
    int *vtxdist_c = data_mesh_c->vtxdist;

    data_nullspace_c->size_global = data_agg->np;
    data_nullspace_c->size_local = vtxdist_c[my_rank + 1] - vtxdist_c[my_rank];
    int size_local_c = data_nullspace_c->size_local;

    PetscCall(SAMGLevelKGhostDataNearNullSpace(vtxdist_f, data_nullspace_f, data_agg));
#if 0
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, near null space of ghost data\n", index_r);
            printf("local_vtx_id\t global_vtx_id\t near_null_space\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == my_rank)
                {
                    for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
                    {
                        printf("%d\t %d\t ", index_v, data_agg->fine_global_nullspace[index_p][index_v].idx);
                        int nrow = data_agg->fine_global_nullspace[index_p][index_v].nrow;
                        int ncol = data_agg->fine_global_nullspace[index_p][index_v].ncol;
                        for (int index_ns_r = 0; index_ns_r < nrow; ++index_ns_r)
                        {
                            printf(" \t \t ");
                            for (int index_ns_c = 0; index_ns_c < ncol; ++index_ns_c)
                            {
                                printf("%021.16le\t ", data_agg->fine_global_nullspace[index_p][index_v].val[MAT_COL_MAJOR(index_ns_r, index_ns_c, nrow)]);
                            }
                            printf("\n");
                        }
                    }
                }
            }

            puts("\n");
        }
        fflush(stdout);
    }
#endif // check ghost data near null space

    // QR factorization for each partition
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        if (data_agg->owner[index_p] == my_rank)
        {
            NearNullSpaceDataVertexLevelK *fine_global_nullspace = data_agg->fine_global_nullspace[index_p];
            int tmp_ncol = fine_global_nullspace[0].ncol;
            int tmp_nrow = 0;
            for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
            {
                tmp_nrow += fine_global_nullspace[index_v].nrow;
            }
            data_agg->data_ghost_agg[index_p].nrow = tmp_nrow;
            data_agg->data_ghost_agg[index_p].ncol = tmp_ncol;

            data_agg->data_ghost_agg[index_p].mat_t = (double *)calloc(tmp_nrow * tmp_ncol, sizeof(double));
            assert(data_agg->data_ghost_agg[index_p].mat_t);
            double *mat_t = data_agg->data_ghost_agg[index_p].mat_t;
            int cnt_row = 0;
            for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
            {
                int ns_nrow = fine_global_nullspace[index_v].nrow;
                int ns_ncol = fine_global_nullspace[index_v].ncol;
                for (int index_ns_r = 0; index_ns_r < ns_nrow; ++index_ns_r)
                {
                    for (int index_ns_c = 0; index_ns_c < ns_ncol; ++index_ns_c)
                    {
                        mat_t[MAT_COL_MAJOR(cnt_row, index_ns_c, tmp_nrow)] = fine_global_nullspace[index_v].val[MAT_COL_MAJOR(index_ns_r, index_ns_c, ns_nrow)];
                    }
                    ++cnt_row;
                }
            }
            assert(cnt_row == tmp_nrow);

            data_agg->data_ghost_agg[index_p].mat_q = (double *)malloc(tmp_nrow * tmp_ncol * sizeof(double));
            data_agg->data_ghost_agg[index_p].mat_r = (double *)calloc(tmp_ncol * tmp_ncol, sizeof(double));
            assert(data_agg->data_ghost_agg[index_p].mat_q && data_agg->data_ghost_agg[index_p].mat_r);
            double *mat_q = data_agg->data_ghost_agg[index_p].mat_q;
            double *mat_r = data_agg->data_ghost_agg[index_p].mat_r;

            memcpy(mat_q, mat_t, tmp_nrow * tmp_ncol * sizeof(double));
            QRFactorizationOpenBLAS(tmp_nrow, tmp_ncol, mat_q, mat_r);
        }
    }
#if 0
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, block near null space matrix\n", index_r);
            printf("local_vtx_id\t global_vtx_id\t block_mat\n");
            int coarse_vtx_start = vtxdist_c[my_rank];
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == my_rank)
                {
                    printf("%d\t %d\t ", index_p - coarse_vtx_start, index_p);
                    int block_mat_nrow = data_agg->data_ghost_agg[index_p].nrow;
                    int block_mat_ncol = data_agg->data_ghost_agg[index_p].ncol;
                    // printf("%d\t %d\n", block_mat_nrow, block_mat_ncol);
                    for (int index_block_mat_r = 0; index_block_mat_r < block_mat_nrow; ++index_block_mat_r)
                    {
                        printf(" \t \t ");
                        for (int index_block_mat_c = 0; index_block_mat_c < block_mat_ncol; ++index_block_mat_c)
                        {
                            printf("%021.16le\t ", data_agg->data_ghost_agg[index_p].mat_t[MAT_COL_MAJOR(index_block_mat_r, index_block_mat_c, block_mat_nrow)]);
                        }
                        printf("\n");
                    }
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check partition near null space block matrix

    data_nullspace_c->data_nullspace = (NearNullSpaceDataVertexLevelK *)malloc(size_local_c * sizeof(NearNullSpaceDataVertexLevelK));
    assert(data_nullspace_c->data_nullspace);

    for (int index = 0; index < size_local_c; ++index)
    {
        data_nullspace_c->data_nullspace[index].idx = data_mesh_c->data_vtx.idx[index];
        data_nullspace_c->data_nullspace[index].nrow = 6;
        data_nullspace_c->data_nullspace[index].ncol = 6;

        // copying R from QR factorization of mat_t
        int coarse_vtx_start = vtxdist_c[my_rank];
        memcpy(data_nullspace_c->data_nullspace[index].val,
               data_agg->data_ghost_agg[index + coarse_vtx_start].mat_r,
               36 * sizeof(double));
    }
#if 0
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        (void)MPI_Barrier(comm);
        if (index_r == my_rank)
        {
            printf(">>>> in rank %d, coarse level near null space\n", index_r);
            printf("local_vtx_id\t global_vtx_id\t near_null_space\n");
            for (int index_v = 0; index_v < size_local_c; ++index_v)
            {
                int coarse_vtx_start = vtxdist_c[my_rank];
                printf("%d\t %d\t ", index_v, index_v + coarse_vtx_start);
                int ns_nrow = data_nullspace_c->data_nullspace[index_v].nrow;
                int ns_ncol = data_nullspace_c->data_nullspace[index_v].ncol;
                for (int index_ns_r = 0; index_ns_r < ns_nrow; ++index_ns_r)
                {
                    printf(" \t \t ");
                    for (int index_ns_c = 0; index_ns_c < ns_ncol; ++index_ns_c)
                    {
                        printf("%021.16le\t ", data_nullspace_c->data_nullspace[index_v].val[MAT_COL_MAJOR(index_ns_r, index_ns_c, ns_nrow)]);
                    }
                    printf("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // coarse level near null space data

    return 0;
}

int SAMGLevelKGhostDataNearNullSpace(const int *vtxdist_f /*fine-lelve mesh vtxdist array*/,
                                     NearNullSpaceLevelK *data_nullspace_f /*fine-level near null space*/,
                                     AggData *data_agg /*aggregation data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    int vtx_start_f = vtxdist_f[my_rank];
    size_t unit_size = sizeof(NearNullSpaceDataVertexLevelK);

    data_agg->fine_global_nullspace = (NearNullSpaceDataVertexLevelK **)malloc(data_agg->np * sizeof(NearNullSpaceDataVertexLevelK *));
    assert(data_agg->fine_global_nullspace);

    int *part_fill_offset = (int *)calloc(data_agg->np, sizeof(int));

    // local near null space data
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            data_agg->fine_global_nullspace[index] = (NearNullSpaceDataVertexLevelK *)malloc(data_agg->n_fine[index] * sizeof(NearNullSpaceDataVertexLevelK));
            assert(data_agg->fine_global_nullspace[index]);

            for (int index_v = 0; index_v < data_agg->nlocal[index]; ++index_v)
            {
                int gid = data_agg->local_gids[index][index_v];
                assert(gid == data_agg->data_ghost_agg[index].flid2fgid[index_v]);
                NearNullSpaceDataVertexLevelK *src = &(data_nullspace_f->data_nullspace[gid - vtx_start_f]);
                NearNullSpaceDataVertexLevelK *dest = &(data_agg->fine_global_nullspace[index][index_v]);

                memcpy(dest, src, unit_size);
            }

            part_fill_offset[index] = data_agg->nlocal[index];
        }
        else
        {
            data_agg->fine_global_nullspace[index] = NULL;
            part_fill_offset[index] = 0;
        }
    }

    // // redo nlocal_all in fine_gids constructor
    // int *nlocal_all = (int *)malloc(nprocs * data_agg->np * sizeof(int));
    // assert(nlocal_all);
    // (void)MPI_Allgather(data_agg->nlocal, data_agg->np, MPI_INT,
    //                     nlocal_all, data_agg->np, MPI_INT, comm);
    int *nlocal_all = data_agg->nlocal_all;

    int *send_counts = (int *)calloc(nprocs, sizeof(int));
    int *recv_counts = (int *)calloc(nprocs, sizeof(int));
    int *sdispls = (int *)malloc(nprocs * sizeof(int));
    int *rdispls = (int *)malloc(nprocs * sizeof(int));
    assert(send_counts && recv_counts && sdispls && rdispls);

    // sending counts
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        int target = data_agg->owner[index_p];
        if (target != my_rank && data_agg->nlocal[index_p] > 0)
        {
            send_counts[target] += data_agg->nlocal[index_p] * unit_size;
        }
    }

    // receving counts
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        if (index_r == my_rank)
        {
            continue;
        }
        for (int index_p = 0; index_p < data_agg->np; ++index_p)
        {
            if (data_agg->owner[index_p] == my_rank)
            {
                recv_counts[index_r] += nlocal_all[index_r * data_agg->np + index_p] * unit_size;
            }
        }
    }

    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int index = 1; index < nprocs; ++index)
    {
        sdispls[index] = sdispls[index - 1] + send_counts[index - 1];
        rdispls[index] = rdispls[index - 1] + recv_counts[index - 1];
    }

    // data packaging
    long long total_send_bytes = sdispls[nprocs - 1] + send_counts[nprocs - 1];
    char *send_buf = (char *)malloc((total_send_bytes > 0 ? total_send_bytes : 1));
    int *tmp_offset = (int *)calloc(nprocs, sizeof(int));
    assert(send_buf && tmp_offset);

    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        int target = data_agg->owner[index_p];
        if (target != my_rank && data_agg->nlocal[index_p] > 0)
        {
            for (int index_v = 0; index_v < data_agg->nlocal[index_p]; ++index_v)
            {
                int gid = data_agg->local_gids[index_p][index_v];
                NearNullSpaceDataVertexLevelK *src = &(data_nullspace_f->data_nullspace[gid - vtx_start_f]);

                memcpy(send_buf + sdispls[target] + tmp_offset[target],
                       src,
                       unit_size);

                tmp_offset[target] += unit_size;
            }
        }
    }

    // mpi-all2allv
    long long total_recv_bytes = rdispls[nprocs - 1] + recv_counts[nprocs - 1];
    char *recv_buf = (char *)malloc((total_recv_bytes > 0 ? total_recv_bytes : 1));
    assert(recv_buf);

    (void)MPI_Alltoallv(send_buf, send_counts, sdispls, MPI_BYTE,
                        recv_buf, recv_counts, rdispls, MPI_BYTE, comm);

    // data unpackaging
    memset(tmp_offset, 0, nprocs * sizeof(int));
    for (int index_r = 0; index_r < nprocs; ++index_r)
    {
        if (index_r == my_rank)
        {
            continue;
        }

        for (int index_p = 0; index_p < data_agg->np; ++index_p)
        {
            int cnt = nlocal_all[index_r * data_agg->np + index_p];
            if (cnt > 0 && data_agg->owner[index_p] == my_rank)
            {
                if (part_fill_offset[index_p] + cnt > data_agg->n_fine[index_p])
                {
                    printf("Rank %d CRITICAL: Partition %d overflow! Max: %d, Current: %d, Adding: %d\n",
                           my_rank, index_p, data_agg->n_fine[index_p], part_fill_offset[index_p], cnt);
                    MPI_Abort(comm, 1);
                }

                char *src_ptr = recv_buf + rdispls[index_r] + tmp_offset[index_r];

                NearNullSpaceDataVertexLevelK *dest_ptr = &(data_agg->fine_global_nullspace[index_p][part_fill_offset[index_p]]);

                memcpy(dest_ptr, src_ptr, cnt * unit_size);

                tmp_offset[index_r] += cnt * unit_size;
                part_fill_offset[index_p] += cnt;
            }
        }
    }

    // free memory
    free(recv_buf);
    free(send_buf);
    free(tmp_offset);
    free(send_counts);
    free(recv_counts);
    free(sdispls);
    free(rdispls);
    // free(nlocal_all);
    free(part_fill_offset);

    return 0;
}

int SAMGLevel0NearNullSpace(NearNullSpaceLevel0 *data_nullspace_level0 /*initial near null space*/,
                            NearNullSpaceLevelK *data_nullspace_levelk /*level 0 near null space*/)
{
    data_nullspace_levelk->size_global = data_nullspace_level0->size_global;
    data_nullspace_levelk->size_local = data_nullspace_level0->size_local;

    int size_local = data_nullspace_levelk->size_local;
    data_nullspace_levelk->data_nullspace = (NearNullSpaceDataVertexLevelK *)malloc(size_local * sizeof(NearNullSpaceDataVertexLevelK));
    assert(data_nullspace_levelk);

    for (int index = 0; index < size_local; ++index)
    {
        data_nullspace_levelk->data_nullspace[index].idx = data_nullspace_level0->data_nullspace[index].idx;
        data_nullspace_levelk->data_nullspace[index].nrow = data_nullspace_level0->data_nullspace[index].nrow;
        data_nullspace_levelk->data_nullspace[index].ncol = data_nullspace_level0->data_nullspace[index].ncol;

        memcpy(data_nullspace_levelk->data_nullspace[index].val,
               data_nullspace_level0->data_nullspace[index].val,
               6 * 6 * sizeof(double));
    }

    return 0;
}

int SAMGInitialNearNullSpace(MeshData *data_f_mesh /*fine level 0 mesh data*/,
                             NearNullSpaceLevel0 *data_nullspace_level0 /*level 0 near null space*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    data_nullspace_level0->size_global = data_f_mesh->nv;
    data_nullspace_level0->size_local = data_f_mesh->vtxdist[my_rank + 1] - data_f_mesh->vtxdist[my_rank];

    int size_local = data_nullspace_level0->size_local;
    assert(size_local == data_f_mesh->data_vtx.local_nv);

    data_nullspace_level0->data_nullspace = (NearNullSpaceDataVertexLevel0 *)malloc(size_local * sizeof(NearNullSpaceDataVertexLevel0));
    assert(data_nullspace_level0);

    for (int index = 0; index < size_local; ++index)
    {
        data_nullspace_level0->data_nullspace[index].idx = data_f_mesh->data_vtx.idx[index];
        data_nullspace_level0->data_nullspace[index].type = data_f_mesh->data_vtx.type[index];
        int type = data_nullspace_level0->data_nullspace[index].type;
        CoorData tmp_coor_data = data_f_mesh->data_vtx.data_coor[index];
        if (type == 0)
        {
            // shell type
            data_nullspace_level0->data_nullspace[index].nrow = 6;
            data_nullspace_level0->data_nullspace[index].ncol = 6;

            for (int index_tmp = 0; index_tmp < 36; ++index_tmp)
            {
                data_nullspace_level0->data_nullspace[index].val[index_tmp] = 0.;
            }
            double *val = data_nullspace_level0->data_nullspace[index].val;

            // 6c + r, c是列号, r是行号
            // MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))
            /*    0  1  2  3     4     5
             * 0 [1  0  0  0     z     -y  ]
             * 1 [0  1  0  -z    0     x   ]
             * 2 [0  0  1  y     -x    0   ]
             * 3 [0  0  0  0     n_z   -n_y]
             * 4 [0  0  0  -n_z  0     n_x ]
             * 5 [0  0  0  n_y   -n_x  0   ]
             */
            val[MAT_COL_MAJOR(0, 0, 6)] = 1.;
            val[MAT_COL_MAJOR(0, 4, 6)] = tmp_coor_data.z;
            val[MAT_COL_MAJOR(0, 5, 6)] = -tmp_coor_data.y;

            val[MAT_COL_MAJOR(1, 1, 6)] = 1.;
            val[MAT_COL_MAJOR(1, 3, 6)] = -tmp_coor_data.z;
            val[MAT_COL_MAJOR(1, 5, 6)] = tmp_coor_data.x;

            val[MAT_COL_MAJOR(2, 2, 6)] = 1.;
            val[MAT_COL_MAJOR(2, 3, 6)] = tmp_coor_data.y;
            val[MAT_COL_MAJOR(2, 4, 6)] = -tmp_coor_data.x;

            val[MAT_COL_MAJOR(3, 4, 6)] = tmp_coor_data.nz;
            val[MAT_COL_MAJOR(3, 5, 6)] = -tmp_coor_data.ny;

            val[MAT_COL_MAJOR(4, 3, 6)] = -tmp_coor_data.nz;
            val[MAT_COL_MAJOR(4, 5, 6)] = tmp_coor_data.nx;

            val[MAT_COL_MAJOR(5, 3, 6)] = tmp_coor_data.ny;
            val[MAT_COL_MAJOR(5, 4, 6)] = -tmp_coor_data.nx;
        }
        else if (type == 1)
        {
            // solid type
            data_nullspace_level0->data_nullspace[index].nrow = 3;
            data_nullspace_level0->data_nullspace[index].ncol = 6;

            for (int index_tmp = 0; index_tmp < 36; ++index_tmp)
            {
                data_nullspace_level0->data_nullspace[index].val[index_tmp] = 0.;
            }
            double *val = data_nullspace_level0->data_nullspace[index].val;

            // 6c + r, c是列号, r是行号
            // MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))
            /*    0  1  2  3     4     5
             * 0 [1  0  0  0     z     -y  ]
             * 1 [0  1  0  -z    0     x   ]
             * 2 [0  0  1  y     -x    0   ]
             */
            val[MAT_COL_MAJOR(0, 0, 3)] = 1.;
            val[MAT_COL_MAJOR(0, 4, 3)] = tmp_coor_data.z;
            val[MAT_COL_MAJOR(0, 5, 3)] = -tmp_coor_data.y;

            val[MAT_COL_MAJOR(1, 1, 3)] = 1.;
            val[MAT_COL_MAJOR(1, 3, 3)] = -tmp_coor_data.z;
            val[MAT_COL_MAJOR(1, 5, 3)] = tmp_coor_data.x;

            val[MAT_COL_MAJOR(2, 2, 3)] = 1.;
            val[MAT_COL_MAJOR(2, 3, 3)] = tmp_coor_data.y;
            val[MAT_COL_MAJOR(2, 4, 3)] = -tmp_coor_data.x;
        }
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, level 0 near null space data, global size = %d, local size = %d\n", index,
                   data_nullspace_level0->size_global,
                   data_nullspace_level0->size_local);
            printf("local_vtx_id\t global_vtx_id\t type\t near_null_space\n");
            for (int index_v = 0; index_v < size_local; ++index_v)
            {
                printf("%d\t %d\t %d\t ", index_v,
                       data_nullspace_level0->data_nullspace[index_v].idx,
                       data_nullspace_level0->data_nullspace[index_v].type);
                for (int index_ns_r = 0; index_ns_r < data_nullspace_level0->data_nullspace[index_v].nrow; ++index_ns_r)
                {
                    printf(" \t \t \t ");
                    for (int index_ns_c = 0; index_ns_c < data_nullspace_level0->data_nullspace[index_v].ncol; ++index_ns_c)
                    {
                        printf("%021.16le\t ", data_nullspace_level0->data_nullspace[index_v].val[MAT_COL_MAJOR(index_ns_r, index_ns_c, data_nullspace_level0->data_nullspace[index_v].nrow)]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            puts("\n");
        }
        fflush(stdout);
    }

#endif // check level 0 near null space data

    return 0;
}
