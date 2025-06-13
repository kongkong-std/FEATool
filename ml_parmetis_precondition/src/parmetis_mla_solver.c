#include "../include/main.h"

int CoarseLevelGenerator(const AdjDataMesh *fine_graph_data /*fine level graph data*/,
                         AdjDataMesh *coarse_graph_data /*coarse level graph data*/)
{
    coarse_graph_data->nn = fine_graph_data->nparts;
    coarse_graph_data->dim = fine_graph_data->dim;
    coarse_graph_data->nparts = coarse_graph_data->nn / 4;

    coarse_graph_data->coordinates = (double *)calloc(coarse_graph_data->dim * coarse_graph_data->nn, sizeof(double));
    assert(coarse_graph_data->coordinates);

    int *cnt_node_partition = (int *)calloc(coarse_graph_data->nn, sizeof(int));
    assert(cnt_node_partition);

    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t id_part = fine_graph_data->part[index];
        for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
        {
            coarse_graph_data->coordinates[coarse_graph_data->dim * id_part + index_i] += fine_graph_data->coordinates[fine_graph_data->dim * index + index_i];
        }
        ++(cnt_node_partition[id_part]);
    }

    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
        {
            coarse_graph_data->coordinates[coarse_graph_data->dim * index + index_i] /= cnt_node_partition[index];
        }
    }

    // adjacency list generator
    coarse_graph_data->xadj = (idx_t *)calloc(coarse_graph_data->nn + 1, sizeof(idx_t));

    bool **mat_adj = NULL; // adjacency matrix, coarse.nn x coarse.nn
    mat_adj = (bool **)malloc(coarse_graph_data->nn * sizeof(bool *));
    assert(mat_adj);

    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_graph_data->nn, sizeof(bool));
        assert(mat_adj[index]);
    }

    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t index_start = fine_graph_data->xadj[index];
        idx_t index_end = fine_graph_data->xadj[index + 1];
        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            idx_t fine_node_i = fine_graph_data->adjncy[index_i];
            idx_t coarse_node_i = fine_graph_data->part[fine_node_i];
            for (idx_t index_j = index_i + 1; index_j < index_end; ++index_j)
            {
                idx_t fine_node_j = fine_graph_data->adjncy[index_j];
                idx_t coarse_node_j = fine_graph_data->part[fine_node_j];

                mat_adj[coarse_node_i][coarse_node_j] = true;
                mat_adj[coarse_node_j][coarse_node_i] = true;
            }
        }
    }

    // set mat_adj diagonal to 0
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index][index] = false;
    }

    // assign coarse xadj
    for (int index_i = 0; index_i < coarse_graph_data->nn; ++index_i)
    {
        int cnt_tmp = 0;
        for (int index_j = 0; index_j < coarse_graph_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                ++cnt_tmp;
            }
        }
        coarse_graph_data->xadj[index_i + 1] = coarse_graph_data->xadj[index_i] + cnt_tmp;
    }

    // assign coarse adjncy
    coarse_graph_data->adjncy = (idx_t *)malloc(coarse_graph_data->xadj[coarse_graph_data->nn] * sizeof(idx_t));
    assert(coarse_graph_data->adjncy);
    int pos_tmp = 0;
    for (int index_i = 0; index_i < coarse_graph_data->nn; ++index_i)
    {
        for (int index_j = 0; index_j < coarse_graph_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                coarse_graph_data->adjncy[pos_tmp] = index_j;
                ++pos_tmp;
            }
        }
    }

    coarse_graph_data->part = (idx_t *)malloc(coarse_graph_data->nn * sizeof(idx_t));
    assert(coarse_graph_data->part);

    // free memeory
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
    free(cnt_node_partition);

    return 0;
}

static int DeepCopyDataMesh2Level0Fine(const AdjDataMesh *data_mesh /*mesh data*/,
                                       AdjDataMesh *data_fine_level /*level 0 fine*/)
{
    data_fine_level->nn = data_mesh->nn;
    data_fine_level->dim = data_mesh->dim;
    data_fine_level->nparts = data_mesh->nparts;

    data_fine_level->coordinates = (double *)malloc(data_fine_level->dim *
                                                    data_fine_level->nn *
                                                    sizeof(double));
    assert(data_fine_level);
    memcpy(data_fine_level->coordinates,
           data_mesh->coordinates,
           data_mesh->dim * data_mesh->nn * sizeof(double));

    data_fine_level->xadj = (idx_t *)malloc((data_fine_level->nn + 1) *
                                            sizeof(idx_t));
    assert(data_fine_level->xadj);
    memcpy(data_fine_level->xadj,
           data_mesh->xadj,
           (data_mesh->nn + 1) * sizeof(idx_t));

    data_fine_level->adjncy = (idx_t *)malloc(data_fine_level->xadj[data_fine_level->nn] *
                                              sizeof(idx_t));
    assert(data_fine_level->adjncy);
    memcpy(data_fine_level->adjncy,
           data_mesh->adjncy,
           data_mesh->xadj[data_mesh->nn] * sizeof(idx_t));

    return 0;
}

int ParMetisMLASolverSetupPhase(MLAContext *mla_ctx /*mla context data*/)
{
    if (mla_ctx->setup == 1)
    {
        return 0;
    }

    mla_ctx->setup = 1;

    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    idx_t *local_xadj = NULL, *local_adjncy = NULL;
    idx_t wgtflag = 0;
    idx_t numflag = 0, ncon = 1;

    real_t *tpwgts = NULL, ubvec = 1.05;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    idx_t local_index_start, local_index_end;

    int config_num_level = mla_ctx->config.mla_config.mla_level;
    int cnt_num_level = 0;

    mla_ctx->metis_mla[cnt_num_level].fine = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
    mla_ctx->metis_mla[cnt_num_level].coarse = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
    assert(mla_ctx->metis_mla[cnt_num_level].fine &&
           mla_ctx->metis_mla[cnt_num_level].coarse);

    mla_ctx->metis_mla[cnt_num_level].fine->vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(mla_ctx->metis_mla[cnt_num_level].fine->vtxdist);

    if (my_rank == 0)
    {
        // copy mla_ctx->data_mesh to mla_ctx->metis_mla[0].fine
        DeepCopyDataMesh2Level0Fine(mla_ctx->data_mesh,
                                    mla_ctx->metis_mla[cnt_num_level].fine);
        mla_ctx->metis_mla[cnt_num_level].fine->part = (idx_t *)malloc(mla_ctx->metis_mla[cnt_num_level].fine->nn * sizeof(idx_t));
        assert(mla_ctx->metis_mla[cnt_num_level].fine->part);
#if 0
        puts("\n>>>> information of graph >>>>\n");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->nn; ++index)
        {
            printf("node %d:\t", index);
            for (int index_i = 0; index_i < mla_ctx->metis_mla[cnt_num_level].fine->dim; ++index_i)
            {
                printf("%021.16le\t",
                       mla_ctx->metis_mla[cnt_num_level].fine->coordinates[mla_ctx->metis_mla[cnt_num_level].fine->dim * index + index_i]);
            }
            putchar('\n');
        }

        puts("\n\ngraph_data xadj value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->nn + 1; ++index)
        {
            printf("%" PRIDX "\t", mla_ctx->metis_mla[cnt_num_level].fine->xadj[index]);
        }
        putchar('\n');

        puts("\n\ngraph_data adjncy value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->xadj[mla_ctx->metis_mla[cnt_num_level].fine->nn]; ++index)
        {
            printf("%" PRIDX "\t", mla_ctx->metis_mla[cnt_num_level].fine->adjncy[index]);
        }
        putchar('\n');

        puts("\n>>>>>>>>\n\n");
#endif // graph csr data

        int base_num = mla_ctx->metis_mla[cnt_num_level].fine->nn / nprocs;
        int remainder_num = mla_ctx->metis_mla[cnt_num_level].fine->nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p + 1] = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p] + count;
        }

        local_index_start = mla_ctx->metis_mla[cnt_num_level].fine->xadj[mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[my_rank]];
        local_index_end = mla_ctx->metis_mla[cnt_num_level].fine->xadj[mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[my_rank + 1]];

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].fine->xadj + mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p], sizeof(idx_t), MPI_BYTE, index_p, 0, comm);     // index_start
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].fine->xadj + mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p + 1], sizeof(idx_t), MPI_BYTE, index_p, 1, comm); // index_end
        }
    }
    else
    {
        (void)MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE, 0, 0, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
    }

    (void)MPI_Bcast(&mla_ctx->metis_mla[cnt_num_level].fine->nn, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(mla_ctx->metis_mla[cnt_num_level].fine->vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    idx_t nparts = mla_ctx->metis_mla[cnt_num_level].fine->nn / 4;
    idx_t edgecut = 0;
    idx_t num_local_node = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[my_rank + 1] -
                           mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[my_rank];
    idx_t nnz_local_node = local_index_end - local_index_start;
    idx_t *local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(local_part && tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / ncon / nparts;
    }

    local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(local_xadj && local_adjncy);

    if (my_rank == 0)
    {
        memcpy(local_xadj, mla_ctx->metis_mla[cnt_num_level].fine->xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(local_adjncy, mla_ctx->metis_mla[cnt_num_level].fine->adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p + 1] -
                                       mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p];
            idx_t tmp_local_index_start = mla_ctx->metis_mla[cnt_num_level].fine->xadj[mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p]];
            idx_t tmp_local_index_end = mla_ctx->metis_mla[cnt_num_level].fine->xadj[mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].fine->xadj + mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p],
                           (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                           index_p, 1, comm); // xadj data
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].fine->adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                           index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        (void)MPI_Recv(local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    // local xadj start from 0
    idx_t tmp_shift = local_xadj[0];
    for (int i = 0; i <= num_local_node; ++i)
    {
        local_xadj[i] -= tmp_shift;
    }

    // partition
    int metis_status = ParMETIS_V3_PartKway(mla_ctx->metis_mla[cnt_num_level].fine->vtxdist,
                                            local_xadj, local_adjncy,
                                            NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, local_part, &comm);

    // gather to root rank
    int *recvcounts = NULL, *displs = NULL;
    if (my_rank == 0)
    {
        recvcounts = (int *)malloc(nprocs * sizeof(int));
        displs = (int *)malloc(nprocs * sizeof(int));
        assert(recvcounts && displs);

        for (int index = 0; index < nprocs; ++index)
        {
#if 1
            recvcounts[index] = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index + 1] -
                                mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index];
            displs[index] = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index];
#endif // int
        }
    }

#if 1
    (void)MPI_Gatherv(local_part, num_local_node, MPI_INT,
                      mla_ctx->metis_mla[cnt_num_level].fine->part,
                      recvcounts, displs, MPI_INT, 0, comm);
#endif // int

    // coarse level
    mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist);

    if (my_rank == 0)
    {
        CoarseLevelGenerator(mla_ctx->metis_mla[cnt_num_level].fine,
                             mla_ctx->metis_mla[cnt_num_level].coarse);

        int base_num = mla_ctx->metis_mla[cnt_num_level].coarse->nn / nprocs;
        int remainder_num = mla_ctx->metis_mla[cnt_num_level].coarse->nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p + 1] = mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p] + count;
        }

        local_index_start = mla_ctx->metis_mla[cnt_num_level].coarse->xadj[mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[my_rank]];
        local_index_end = mla_ctx->metis_mla[cnt_num_level].coarse->xadj[mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[my_rank + 1]];
        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].coarse->xadj + mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p],
                           sizeof(idx_t), MPI_BYTE, index_p, 10, comm); // coarse level local index_start
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].coarse->xadj + mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p + 1],
                           sizeof(idx_t), MPI_BYTE, index_p, 11, comm); // coarse level local index_end
        }
    }
    else
    {
        (void)MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE,
                       0, 10, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE,
                       0, 11, comm, MPI_STATUS_IGNORE);
    }

    (void)MPI_Bcast(&mla_ctx->metis_mla[cnt_num_level].coarse->nn, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    nparts = mla_ctx->metis_mla[cnt_num_level].coarse->nn / 4;
    edgecut = 0;
    num_local_node = mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[my_rank + 1] -
                     mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[my_rank];
    nnz_local_node = local_index_end - local_index_start;
    ncon = 1;
    numflag = 0;
    wgtflag = 0;
    idx_t *coarse_local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    real_t *coarse_tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    ubvec = 1.05;
    METIS_SetDefaultOptions(options);
    assert(coarse_local_part && coarse_tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        coarse_tpwgts[index] = 1. / ncon / nparts;
    }

    idx_t *coarse_local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    idx_t *coarse_local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(coarse_local_xadj && coarse_local_adjncy);

    if (my_rank == 0)
    {
        memcpy(coarse_local_xadj, mla_ctx->metis_mla[cnt_num_level].coarse->xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(coarse_local_adjncy, mla_ctx->metis_mla[cnt_num_level].coarse->adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p + 1] - mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p];
            idx_t tmp_local_index_start = mla_ctx->metis_mla[cnt_num_level].coarse->xadj[mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p]];
            idx_t tmp_local_index_end = mla_ctx->metis_mla[cnt_num_level].coarse->xadj[mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].coarse->xadj + mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index_p], (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                           index_p, 1, comm); // xadj data
            (void)MPI_Send(mla_ctx->metis_mla[cnt_num_level].coarse->adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                           index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        (void)MPI_Recv(coarse_local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(coarse_local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    tmp_shift = coarse_local_xadj[0];
    for (int index = 0; index <= num_local_node; ++index)
    {
        coarse_local_xadj[index] -= tmp_shift;
    }

    int coarse_metis_status = ParMETIS_V3_PartKway(mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist,
                                                   coarse_local_xadj, coarse_local_adjncy,
                                                   NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                                   coarse_tpwgts, &ubvec, options,
                                                   &edgecut, coarse_local_part, &comm);

    int *coarse_recvcounts = NULL, *coarse_displs = NULL;
    if (my_rank == 0)
    {
        coarse_recvcounts = (int *)malloc(nprocs * sizeof(int));
        coarse_displs = (int *)malloc(nprocs * sizeof(int));
        assert(coarse_displs && coarse_displs);

        for (int index = 0; index < nprocs; ++index)
        {
#if 1
            coarse_recvcounts[index] = mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index + 1] - mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index];
            coarse_displs[index] = mla_ctx->metis_mla[cnt_num_level].coarse->vtxdist[index];
#endif // int
#if 0
            recvcounts[index] = (fine_graph_data.vtxdist[index + 1] - fine_graph_data.vtxdist[index]) * sizeof(idx_t);
            displs[index] = (fine_graph_data.vtxdist[index]) * sizeof(idx_t);
#endif // byte
        }
    }

#if 1
    MPI_Gatherv(coarse_local_part, num_local_node, MPI_INT,
                mla_ctx->metis_mla[cnt_num_level].coarse->part, coarse_recvcounts, coarse_displs, MPI_INT, 0, comm);
#endif // int

#if 1
    if (my_rank == 0)
    {
        puts("\ncoarse global partition value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].coarse->part[index]);
        }

        puts("\n========\n");
    }
#endif // print global partition value, coarse level

#if 0
    if (my_rank == 0)
    {
        puts("\nglobal partition value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].fine->part[index]);
        }

        puts("\n========\n");
    }
#endif // print global partition value

#if 0
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            puts("\n==== mla_ctx information ====\n");
            printf("in rank %d/%d\n", my_rank, nprocs);
            printf("setup = %d\n", mla_ctx->setup);
            printf("path_config = %s\n", mla_ctx->path_config);
            printf("order_rbm = %d\n", mla_ctx->order_rbm);
            printf("angle_type = %d\n", mla_ctx->angle_type);

            puts("\n--------\n");

            printf("file_mat = %s\n", mla_ctx->config.file_config.file_mat);
            printf("file_rhs = %s\n", mla_ctx->config.file_config.file_rhs);
            printf("file_mesh = %s\n", mla_ctx->config.file_config.file_mesh);
            printf("pre_smooth_v = %d\n", mla_ctx->config.mla_config.pre_smooth_v);
            printf("post_smooth_v = %d\n", mla_ctx->config.mla_config.post_smooth_v);
            printf("mla_max_it = %d\n", mla_ctx->config.mla_config.mla_max_it);
            printf("mla_rtol = %021.1le\n", mla_ctx->config.mla_config.mla_rtol);
            printf("mla_level = %d\n", mla_ctx->config.mla_config.mla_level);
            printf("mla_phase = %d\n", mla_ctx->config.mla_config.mla_phase);
            printf("coarse_restart = %d\n", mla_ctx->config.mla_config.coarse_restart);
            printf("label_bound = %d\n", mla_ctx->config.mesh_label_config.label_bound);
            printf("label_omega = %d\n", mla_ctx->config.mesh_label_config.label_omega);

            puts("\n--------\n");
        }
    }
#endif // mla_ctx data information

    return 0;
}

int ParMetisMLASolver(MLAContext *mla_ctx /*mla context data*/,
                      int mla_phase /*multilevel phase*/)
{
    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        ParMetisMLASolverSetupPhase(mla_ctx);
    }

    // solve phase
    if (mla_phase == 1 || mla_phase == 2)
    {
    }

    return 0;
}
