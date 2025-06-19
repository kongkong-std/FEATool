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

int Level0FinePartition(const AdjDataMesh *data_mesh /*initial data mesh*/,
                        AdjDataMesh *fine_graph /*level 0 fine graph data*/)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    idx_t wgtflag = 0;
    idx_t numflag = 0, ncon = 1;

    real_t *tpwgts = NULL, ubvec = 1.05;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    idx_t local_index_start, local_index_end;

    fine_graph->vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(fine_graph->vtxdist);

    if (my_rank == 0)
    {
        DeepCopyDataMesh2Level0Fine(data_mesh, fine_graph);
        fine_graph->part = (idx_t *)malloc(fine_graph->nn * sizeof(idx_t));
        assert(fine_graph->part);
#if 0
        puts("\n>>>> information of graph >>>>\n");
        for (int index = 0; index < fine_graph->nn; ++index)
        {
            printf("node %d:\t", index);
            for (int index_i = 0; index_i < fine_graph->dim; ++index_i)
            {
                printf("%021.16le\t",
                       fine_graph->coordinates[fine_graph->dim * index + index_i]);
            }
            putchar('\n');
        }

        puts("\n\ngraph_data xadj value:");
        for (int index = 0; index < fine_graph->nn + 1; ++index)
        {
            printf("%" PRIDX "\t", fine_graph->xadj[index]);
        }
        putchar('\n');

        puts("\n\ngraph_data adjncy value:");
        for (int index = 0; index < fine_graph->xadj[fine_graph->nn]; ++index)
        {
            printf("%" PRIDX "\t", fine_graph->adjncy[index]);
        }
        putchar('\n');

        puts("\n>>>>>>>>\n\n");
#endif // graph csr data

        int base_num = fine_graph->nn / nprocs;
        int remainder_num = fine_graph->nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            fine_graph->vtxdist[index_p + 1] = fine_graph->vtxdist[index_p] + count;
        }

        local_index_start = fine_graph->xadj[fine_graph->vtxdist[my_rank]];
        local_index_end = fine_graph->xadj[fine_graph->vtxdist[my_rank + 1]];

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            (void)MPI_Send(fine_graph->xadj + fine_graph->vtxdist[index_p], sizeof(idx_t), MPI_BYTE, index_p, 0, comm);     // index_start
            (void)MPI_Send(fine_graph->xadj + fine_graph->vtxdist[index_p + 1], sizeof(idx_t), MPI_BYTE, index_p, 1, comm); // index_end
        }
    }
    else
    {
        (void)MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE, 0, 0, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
    }

    (void)MPI_Bcast(&fine_graph->nn, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&fine_graph->dim, 1, MPI_INT, 0, comm);

    if (my_rank != 0)
    {
        fine_graph->coordinates = (double *)malloc(fine_graph->dim * fine_graph->nn * sizeof(double));
        assert(fine_graph->coordinates);
    }

    (void)MPI_Bcast(fine_graph->coordinates, fine_graph->dim * fine_graph->nn, MPI_DOUBLE, 0, comm);
    (void)MPI_Bcast(fine_graph->vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    idx_t nparts = fine_graph->nn / 4;
    idx_t edgecut = 0;
    idx_t num_local_node = fine_graph->vtxdist[my_rank + 1] -
                           fine_graph->vtxdist[my_rank];
    idx_t nnz_local_node = local_index_end - local_index_start;
    fine_graph->local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(fine_graph->local_part && tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / ncon / nparts;
    }

    fine_graph->local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    fine_graph->local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(fine_graph->local_xadj && fine_graph->local_adjncy);

    if (my_rank == 0)
    {
        memcpy(fine_graph->local_xadj, fine_graph->xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(fine_graph->local_adjncy, fine_graph->adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = fine_graph->vtxdist[index_p + 1] -
                                       fine_graph->vtxdist[index_p];
            idx_t tmp_local_index_start = fine_graph->xadj[fine_graph->vtxdist[index_p]];
            idx_t tmp_local_index_end = fine_graph->xadj[fine_graph->vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            (void)MPI_Send(fine_graph->xadj + fine_graph->vtxdist[index_p],
                           (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                           index_p, 1, comm); // xadj data
            (void)MPI_Send(fine_graph->adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                           index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        (void)MPI_Recv(fine_graph->local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(fine_graph->local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    // local xadj start from 0
    idx_t tmp_shift = fine_graph->local_xadj[0];
    for (int i = 0; i <= num_local_node; ++i)
    {
        fine_graph->local_xadj[i] -= tmp_shift;
    }

    // partition
    int metis_status = ParMETIS_V3_PartKway(fine_graph->vtxdist,
                                            fine_graph->local_xadj, fine_graph->local_adjncy,
                                            NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, fine_graph->local_part, &comm);

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
            recvcounts[index] = fine_graph->vtxdist[index + 1] -
                                fine_graph->vtxdist[index];
            displs[index] = fine_graph->vtxdist[index];
#endif // int
        }
    }

#if 1
    (void)MPI_Gatherv(fine_graph->local_part, num_local_node, MPI_INT,
                      fine_graph->part,
                      recvcounts, displs, MPI_INT, 0, comm);
#endif // int

    return 0;
}

int LevelKCoarsePartition(const AdjDataMesh *fine_graph /*level k fine data mesh*/,
                          AdjDataMesh *coarse_graph /*level k coarse data mesh*/)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    idx_t local_index_start, local_index_end;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    coarse_graph->vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(coarse_graph->vtxdist);

    if (my_rank == 0)
    {
        CoarseLevelGenerator(fine_graph,
                             coarse_graph);

        int base_num = coarse_graph->nn / nprocs;
        int remainder_num = coarse_graph->nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            coarse_graph->vtxdist[index_p + 1] = coarse_graph->vtxdist[index_p] + count;
        }

        local_index_start = coarse_graph->xadj[coarse_graph->vtxdist[my_rank]];
        local_index_end = coarse_graph->xadj[coarse_graph->vtxdist[my_rank + 1]];
        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            (void)MPI_Send(coarse_graph->xadj + coarse_graph->vtxdist[index_p],
                           sizeof(idx_t), MPI_BYTE, index_p, 10, comm); // coarse level local index_start
            (void)MPI_Send(coarse_graph->xadj + coarse_graph->vtxdist[index_p + 1],
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

    (void)MPI_Bcast(&coarse_graph->nn, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&coarse_graph->dim, 1, MPI_INT, 0, comm);

    if (my_rank != 0)
    {
        coarse_graph->coordinates = (double *)malloc(coarse_graph->dim * coarse_graph->nn * sizeof(double));
        assert(coarse_graph->coordinates);
    }

    (void)MPI_Bcast(coarse_graph->coordinates, coarse_graph->dim * coarse_graph->nn, MPI_DOUBLE, 0, comm);
    (void)MPI_Bcast(coarse_graph->vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    idx_t nparts = coarse_graph->nn / 4;
    idx_t edgecut = 0;
    idx_t num_local_node = coarse_graph->vtxdist[my_rank + 1] -
                           coarse_graph->vtxdist[my_rank];
    idx_t nnz_local_node = local_index_end - local_index_start;
    idx_t ncon = 1;
    idx_t numflag = 0;
    idx_t wgtflag = 0;

    coarse_graph->local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    real_t *coarse_tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    real_t ubvec = 1.05;

    assert(coarse_graph->local_part && coarse_tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        coarse_tpwgts[index] = 1. / ncon / nparts;
    }

    coarse_graph->local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    coarse_graph->local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(coarse_graph->local_xadj && coarse_graph->local_adjncy);

    if (my_rank == 0)
    {
        memcpy(coarse_graph->local_xadj, coarse_graph->xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(coarse_graph->local_adjncy, coarse_graph->adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = coarse_graph->vtxdist[index_p + 1] - coarse_graph->vtxdist[index_p];
            idx_t tmp_local_index_start = coarse_graph->xadj[coarse_graph->vtxdist[index_p]];
            idx_t tmp_local_index_end = coarse_graph->xadj[coarse_graph->vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            (void)MPI_Send(coarse_graph->xadj + coarse_graph->vtxdist[index_p], (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                           index_p, 1, comm); // xadj data
            (void)MPI_Send(coarse_graph->adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                           index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        (void)MPI_Recv(coarse_graph->local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(coarse_graph->local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    idx_t tmp_shift = coarse_graph->local_xadj[0];
    for (int index = 0; index <= num_local_node; ++index)
    {
        coarse_graph->local_xadj[index] -= tmp_shift;
    }

    int coarse_metis_status = ParMETIS_V3_PartKway(coarse_graph->vtxdist,
                                                   coarse_graph->local_xadj, coarse_graph->local_adjncy,
                                                   NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                                   coarse_tpwgts, &ubvec, options,
                                                   &edgecut, coarse_graph->local_part, &comm);

    int *coarse_recvcounts = NULL, *coarse_displs = NULL;
    if (my_rank == 0)
    {
        coarse_recvcounts = (int *)malloc(nprocs * sizeof(int));
        coarse_displs = (int *)malloc(nprocs * sizeof(int));
        assert(coarse_displs && coarse_displs);

        for (int index = 0; index < nprocs; ++index)
        {
#if 1
            coarse_recvcounts[index] = coarse_graph->vtxdist[index + 1] - coarse_graph->vtxdist[index];
            coarse_displs[index] = coarse_graph->vtxdist[index];
#endif // int
        }
    }

#if 1
    (void)MPI_Gatherv(coarse_graph->local_part, num_local_node, MPI_INT,
                      coarse_graph->part, coarse_recvcounts, coarse_displs, MPI_INT, 0, comm);
#endif // int

    return 0;
}

int DeepCopyCoarse2NextLevelFine(const AdjDataMesh *coarse_graph /*coarse level graph*/,
                                 AdjDataMesh *fine_graph /*next level fine graph*/)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    idx_t local_index_start, local_index_end;

    fine_graph->vtxdist = (idx_t *)malloc((nprocs + 1) * sizeof(idx_t));
    assert(fine_graph->vtxdist);

    if (my_rank == 0)
    {
        fine_graph->nn = coarse_graph->nn;
        fine_graph->dim = coarse_graph->dim;
        fine_graph->nparts = coarse_graph->nparts;

        fine_graph->coordinates = (double *)malloc(fine_graph->dim * fine_graph->nn * sizeof(double));
        assert(fine_graph->coordinates);
        memcpy(fine_graph->coordinates,
               coarse_graph->coordinates,
               fine_graph->dim * fine_graph->nn * sizeof(double));

        fine_graph->xadj = (idx_t *)malloc((fine_graph->nn + 1) * sizeof(idx_t));
        assert(fine_graph->xadj);
        memcpy(fine_graph->xadj,
               coarse_graph->xadj,
               (coarse_graph->nn + 1) * sizeof(idx_t));

        fine_graph->adjncy = (idx_t *)malloc(fine_graph->xadj[fine_graph->nn] * sizeof(idx_t));
        assert(fine_graph->adjncy);
        memcpy(fine_graph->adjncy,
               coarse_graph->adjncy,
               coarse_graph->xadj[coarse_graph->nn] * sizeof(idx_t));

        memcpy(fine_graph->vtxdist,
               coarse_graph->vtxdist,
               (nprocs + 1) * sizeof(idx_t));

        local_index_start = fine_graph->xadj[fine_graph->vtxdist[my_rank]];
        local_index_end = fine_graph->xadj[fine_graph->vtxdist[my_rank + 1]];

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            (void)MPI_Send(fine_graph->xadj + fine_graph->vtxdist[index_p], sizeof(idx_t), MPI_BYTE, index_p, 0, comm);     // index_start
            (void)MPI_Send(fine_graph->xadj + fine_graph->vtxdist[index_p + 1], sizeof(idx_t), MPI_BYTE, index_p, 1, comm); // index_end
        }

        fine_graph->part = (idx_t *)malloc(fine_graph->nn * sizeof(idx_t));
        assert(fine_graph->part);
        memcpy(fine_graph->part,
               coarse_graph->part,
               coarse_graph->nn * sizeof(idx_t));
    }
    else
    {
        (void)MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE, 0, 0, comm, MPI_STATUS_IGNORE);
        (void)MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
    }

    (void)MPI_Bcast(&fine_graph->nn, 1, MPI_INT, 0, comm);
    (void)MPI_Bcast(&fine_graph->dim, 1, MPI_INT, 0, comm);

    if (my_rank != 0)
    {
        fine_graph->coordinates = (double *)malloc(fine_graph->dim * fine_graph->nn * sizeof(double));
        assert(fine_graph->coordinates);
    }

    (void)MPI_Bcast(fine_graph->coordinates, fine_graph->dim * fine_graph->nn, MPI_DOUBLE, 0, comm);
    (void)MPI_Bcast(fine_graph->vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    idx_t num_local_node = fine_graph->vtxdist[my_rank + 1] -
                           fine_graph->vtxdist[my_rank];
    idx_t nnz_local_node = local_index_end - local_index_start;
    fine_graph->local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    fine_graph->local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    fine_graph->local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(fine_graph->local_part &&
           fine_graph->local_xadj &&
           fine_graph->local_adjncy);

    memcpy(fine_graph->local_part,
           coarse_graph->local_part,
           num_local_node * sizeof(idx_t));
    memcpy(fine_graph->local_xadj,
           coarse_graph->local_xadj,
           (num_local_node + 1) * sizeof(idx_t));
    memcpy(fine_graph->local_adjncy,
           coarse_graph->local_adjncy,
           nnz_local_node * sizeof(idx_t));

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

    int config_num_level = mla_ctx->config.mla_config.mla_level;
    int cnt_num_level = 0;

    mla_ctx->metis_mla[cnt_num_level].fine = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
    mla_ctx->metis_mla[cnt_num_level].coarse = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
    assert(mla_ctx->metis_mla[cnt_num_level].fine &&
           mla_ctx->metis_mla[cnt_num_level].coarse);

    Level0FinePartition(mla_ctx->data_mesh,
                        mla_ctx->metis_mla[cnt_num_level].fine);

    LevelKCoarsePartition(mla_ctx->metis_mla[cnt_num_level].fine,
                          mla_ctx->metis_mla[cnt_num_level].coarse);

    // multilevel aggreagtion
    for (cnt_num_level = 1; cnt_num_level < config_num_level; ++cnt_num_level)
    {
        if (mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn < 100)
        {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n\n>>>>nodes in level %" PetscInt_FMT " is %" PetscInt_FMT ", less than 100, setup DONE!>>>>\n\n\n",
                                  cnt_num_level - 1,
                                  mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn));

            break;
        }

        mla_ctx->metis_mla[cnt_num_level].fine = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
        mla_ctx->metis_mla[cnt_num_level].coarse = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
        assert(mla_ctx->metis_mla[cnt_num_level].fine &&
               mla_ctx->metis_mla[cnt_num_level].coarse);

        DeepCopyCoarse2NextLevelFine(mla_ctx->metis_mla[cnt_num_level - 1].coarse,
                                     mla_ctx->metis_mla[cnt_num_level].fine);

        LevelKCoarsePartition(mla_ctx->metis_mla[cnt_num_level].fine,
                              mla_ctx->metis_mla[cnt_num_level].coarse);
    }

    mla_ctx->num_level = cnt_num_level;
    mla_ctx->true_num_level = cnt_num_level;

#if 0
    for (int index = 0; index < mla_ctx->num_level; ++index)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "fine level coordinates, nn = %d, dim = %d\n",
                              mla_ctx->metis_mla[index].fine->nn,
                              mla_ctx->metis_mla[index].fine->dim));

        for (int index_i = 0; index_i < mla_ctx->metis_mla[index].fine->nn; ++index_i)
        {
            printf("node %d:\t", index_i);
            for (int index_j = 0; index_j < mla_ctx->metis_mla[index].fine->dim; ++index_j)
            {
                printf("%021.16le\t", mla_ctx->metis_mla[index].fine->coordinates[mla_ctx->metis_mla[index].fine->dim *
                                                                                      index_i +
                                                                                  index_j]);
            }
            putchar('\n');
        }

        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "coarse level coordinates, nn = %d, dim = %d\n",
                              mla_ctx->metis_mla[index].coarse->nn,
                              mla_ctx->metis_mla[index].coarse->dim));
    }
#endif // coordinates information

    // operator assembling
    PetscInt mat_m = 0, mat_n = 0;
    PetscCall(MatGetSize(mla_ctx->mysolver.solver_a, &mat_m, &mat_n));
#if 0
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix size: m = %" PetscInt_FMT ", n = %" PetscInt_FMT "\n", mat_m, mat_n));
#endif // print matrix size information

    PetscCall(MatDuplicate(mla_ctx->mysolver.solver_a,
                           MAT_COPY_VALUES,
                           &mla_ctx->metis_mla[0].operator_fine));

    if (mla_ctx->angle_type == 0)
    {
        // theta, rotational angle with axis
        if (mla_ctx->order_rbm == 1)
        {
            // rbm type
            /*
             * C \in R^{6 x 6}
             * [1  0  0   0   z  -y]
             * [0  1  0  -z   0   x]
             * [0  0  1   y  -x   0]
             * [0  0  0   1   0   0]
             * [0  0  0   0   1   0]
             * [0  0  0   0   0   1]
             *     x = fine_node_x - coarse_node_x
             *     y = fine_node_y - coarse_node_y
             *     z = fine_node_z - coarse_node_z
             */
            PetscInt prolongation_m = 0;
            for (int index_cnt_level = 0; index_cnt_level < cnt_num_level; ++index_cnt_level)
            {
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &prolongation_m, NULL));
                PetscInt prolongation_n = mla_ctx->metis_mla[index_cnt_level].coarse->nn * 6;
                int local_index_start = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank],
                    local_index_end = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank + 1];
                int coarse_local_index_start = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank],
                    coarse_local_index_end = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank + 1];
                PetscInt local_prolongation_m = local_index_end - local_index_start,
                         local_prolongation_n = coarse_local_index_end - coarse_local_index_start;
#if 0
                printf("in rank %d/%d,"
                       " local_prolonation_n = %d,"
                       " global_prolongation_n = %d\n",
                       my_rank, nprocs,
                       local_prolongation_n, prolongation_n);
#endif // prolongation size
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      6 * local_prolongation_m, 6 * local_prolongation_n,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                {
                    double p_loc[6][6] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                    double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                           fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                           fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                    double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                           coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                           coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                    p_loc[0][4] = fine_node_z - coarse_node_z; // z
                    p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                    p_loc[1][3] = coarse_node_z - fine_node_z; // -z
                    p_loc[1][5] = fine_node_x - coarse_node_x; // x
                    p_loc[2][3] = fine_node_y - coarse_node_y; // y
                    p_loc[2][4] = coarse_node_x - fine_node_x; // -x

                    for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                    {
                        for (int index_tmp_j = 0; index_tmp_j < 6; ++index_tmp_j)
                        {
                            PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                  6 * index_fine_node + index_tmp_i,
                                                  6 * index_coarse_node + index_tmp_j,
                                                  p_loc[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));

#if 0
                puts("prolongation operator done!!\n\n");
                PetscLayout rmap_A, rmap_P;
                PetscBool is_same;
                PetscCall(MatGetLayouts(mla_ctx->metis_mla[index_cnt_level].operator_fine, &rmap_A, NULL));
                PetscCall(MatGetLayouts(mla_ctx->metis_mla[index_cnt_level].prolongation, &rmap_P, NULL));
                PetscCall(PetscLayoutCompare(rmap_A, rmap_P, &is_same));
                if (is_same)
                {
                    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "same layout, can call matptap"));
                }
                else
                {
                    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "different layout, cannot call matptap"));
                }
#endif // information of parallel layouts

#if 1
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif
            }
        }

        if (mla_ctx->order_rbm == 2)
        {
            // kirchhoff type
            /*
             * C \in R^{6 x 9}
             * [1  0  0   0   z  -y   zx      0       zy/2]
             * [0  1  0  -z   0   x   0       zy      zx/2]
             * [0  0  1   y  -x   0  -x^2/2  -y^2/2  -xy/2]
             * [0  0  0   1   0   0   0      -y       -x/2]
             * [0  0  0   0   1   0   x       0        y/2]
             * [0  0  0   0   0   1   0       0          0]
             *     x = fine_node_x - coarse_node_x
             *     y = fine_node_y - coarse_node_y
             *     z = fine_node_z - coarse_node_z
             */
            PetscInt prolongation_m = 0;
            for (int index_cnt_level = 0; index_cnt_level < cnt_num_level; ++index_cnt_level)
            {
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &prolongation_m, NULL));
                PetscInt prolongation_n = mla_ctx->metis_mla[index_cnt_level].coarse->nn * 9;
                int local_index_start = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank],
                    local_index_end = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank + 1];
                int coarse_local_index_start = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank],
                    coarse_local_index_end = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank + 1];
                PetscInt local_prolongation_m = local_index_end - local_index_start,
                         local_prolongation_n = coarse_local_index_end - coarse_local_index_start;

                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                if (index_cnt_level == 0)
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                          local_prolongation_m * 6, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                else
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                          local_prolongation_m * 9, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                if (index_cnt_level == 0)
                {
                    for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                    {
                        double p_loc[6][9] = {0};
                        for (int index = 0; index < 6; ++index)
                        {
                            p_loc[index][index] = 1.;
                        }

                        idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                        double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                               fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                               fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                        double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                               coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                               coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                        p_loc[0][4] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2
                        p_loc[1][3] = coarse_node_z - fine_node_z;        // -z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2
                        p_loc[2][3] = fine_node_y - coarse_node_y;        // y
                        p_loc[2][4] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
                        p_loc[3][7] = coarse_node_y - fine_node_y;        // -y
                        p_loc[3][8] = (coarse_node_x - fine_node_x) / 2.; // -x/2
                        p_loc[4][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[4][8] = (fine_node_y - coarse_node_y) / 2.; // y/2

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                      6 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }
                else
                {
                    for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                    {
                        double p_loc[9][9] = {0};
                        for (int index = 0; index < 9; ++index)
                        {
                            p_loc[index][index] = 1.;
                        }

                        idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                        double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                               fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                               fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                        double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                               coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                               coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                        p_loc[0][4] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2
                        p_loc[1][3] = coarse_node_z - fine_node_z;        // -z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2
                        p_loc[2][3] = fine_node_y - coarse_node_y;        // y
                        p_loc[2][4] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
                        p_loc[3][7] = coarse_node_y - fine_node_y;        // -y
                        p_loc[3][8] = (coarse_node_x - fine_node_x) / 2.; // -x/2
                        p_loc[4][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[4][8] = (fine_node_y - coarse_node_y) / 2.; // y/2

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                      9 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
            }
        }
    }

    if (mla_ctx->angle_type == 1)
    {
        // phi, normal rotational angle with axis
        if (mla_ctx->order_rbm == 1)
        {
            // rbm type
            /*
             * C \in R^{6 x 6}
             * [1  0  0   z   0  -y]
             * [0  1  0   0   z   x]
             * [0  0  1  -x  -y   0]
             * [0  0  0   1   0   0]
             * [0  0  0   0   1   0]
             * [0  0  0   0   0   1]
             *     x = fine_node_x - coarse_node_x
             *     y = fine_node_y - coarse_node_y
             *     z = fine_node_z - coarse_node_z
             */
            PetscInt prolongation_m = 0;
            for (int index_cnt_level = 0; index_cnt_level < cnt_num_level; ++index_cnt_level)
            {
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &prolongation_m, NULL));
                PetscInt prolongation_n = mla_ctx->metis_mla[index_cnt_level].coarse->nn * 6;
                int local_index_start = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank],
                    local_index_end = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank + 1];
                int coarse_local_index_start = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank],
                    coarse_local_index_end = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank + 1];
                PetscInt local_prolongation_m = local_index_end - local_index_start,
                         local_prolongation_n = coarse_local_index_end - coarse_local_index_start;
#if 0
                printf("in rank %d/%d,"
                       " local_prolonation_n = %d,"
                       " global_prolongation_n = %d\n",
                       my_rank, nprocs,
                       local_prolongation_n, prolongation_n);
#endif // prolongation size
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      6 * local_prolongation_m, 6 * local_prolongation_n,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                {
                    double p_loc[6][6] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                    double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                           fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                           fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                    double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                           coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                           coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                    p_loc[0][3] = fine_node_z - coarse_node_z; // z
                    p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                    p_loc[1][4] = fine_node_z - coarse_node_z; // z
                    p_loc[1][5] = fine_node_x - coarse_node_x; // x
                    p_loc[2][3] = coarse_node_x - fine_node_x; // -x
                    p_loc[2][4] = coarse_node_y - fine_node_y; // -y

                    for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                    {
                        for (int index_tmp_j = 0; index_tmp_j < 6; ++index_tmp_j)
                        {
                            PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                  6 * index_fine_node + index_tmp_i,
                                                  6 * index_coarse_node + index_tmp_j,
                                                  p_loc[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));

#if 0
                puts("prolongation operator done!!\n\n");
                PetscLayout rmap_A, rmap_P;
                PetscBool is_same;
                PetscCall(MatGetLayouts(mla_ctx->metis_mla[index_cnt_level].operator_fine, &rmap_A, NULL));
                PetscCall(MatGetLayouts(mla_ctx->metis_mla[index_cnt_level].prolongation, &rmap_P, NULL));
                PetscCall(PetscLayoutCompare(rmap_A, rmap_P, &is_same));
                if (is_same)
                {
                    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "same layout, can call matptap"));
                }
                else
                {
                    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "different layout, cannot call matptap"));
                }
#endif // information of parallel layouts

#if 1
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif
            }
        }

        if (mla_ctx->order_rbm == 2)
        {
            // kirchhoff type
            /*
             * C \in R^{6 x 9}
             * [1  0  0   z   0  -y   zx      0       zy/2]
             * [0  1  0   0   z   x   0       zy      zx/2]
             * [0  0  1  -x  -y   0  -x^2/2  -y^2/2  -xy/2]
             * [0  0  0   1   0   0   x       0        y/2]
             * [0  0  0   0   1   0   0       y        x/2]
             * [0  0  0   0   0   1   0       0          0]
             *     x = fine_node_x - coarse_node_x
             *     y = fine_node_y - coarse_node_y
             *     z = fine_node_z - coarse_node_z
             */
            PetscInt prolongation_m = 0;
            for (int index_cnt_level = 0; index_cnt_level < cnt_num_level; ++index_cnt_level)
            {
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &prolongation_m, NULL));
                PetscInt prolongation_n = mla_ctx->metis_mla[index_cnt_level].coarse->nn * 9;
                int local_index_start = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank],
                    local_index_end = mla_ctx->metis_mla[index_cnt_level].fine->vtxdist[my_rank + 1];
                int coarse_local_index_start = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank],
                    coarse_local_index_end = mla_ctx->metis_mla[index_cnt_level].coarse->vtxdist[my_rank + 1];
                PetscInt local_prolongation_m = local_index_end - local_index_start,
                         local_prolongation_n = coarse_local_index_end - coarse_local_index_start;

                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                if (index_cnt_level == 0)
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                          local_prolongation_m * 6, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                else
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                          local_prolongation_m * 9, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                if (index_cnt_level == 0)
                {
                    for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                    {
                        double p_loc[6][9] = {0};
                        for (int index = 0; index < 6; ++index)
                        {
                            p_loc[index][index] = 1.;
                        }

                        idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                        double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                               fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                               fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                        double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                               coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                               coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2
                        p_loc[1][4] = fine_node_z - coarse_node_z;        // z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2
                        p_loc[2][3] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y;        // -y
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                      6 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }
                else
                {
                    for (int index_fine_node = local_index_start; index_fine_node < local_index_end; ++index_fine_node)
                    {
                        double p_loc[9][9] = {0};
                        for (int index = 0; index < 9; ++index)
                        {
                            p_loc[index][index] = 1.;
                        }

                        idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->local_part[index_fine_node - local_index_start];
                        double fine_node_x = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node],
                               fine_node_y = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 1],
                               fine_node_z = mla_ctx->metis_mla[index_cnt_level].fine->coordinates[3 * index_fine_node + 2];

                        double coarse_node_x = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node],
                               coarse_node_y = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 1],
                               coarse_node_z = mla_ctx->metis_mla[index_cnt_level].coarse->coordinates[3 * index_coarse_node + 2];

                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2
                        p_loc[1][4] = fine_node_z - coarse_node_z;        // z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2
                        p_loc[2][3] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y;        // -y
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                                      9 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].prolongation, MAT_FINAL_ASSEMBLY));

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
            }
        }
    }

#if 0
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        (void)MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            puts("\n==== level k fine mesh partition value ====\n");
            printf("in rank %d/%d\n", my_rank, nprocs);
            idx_t index_start = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p];
            idx_t index_end = mla_ctx->metis_mla[cnt_num_level].fine->vtxdist[index_p + 1];
            for(idx_t index = index_start; index < index_end; ++index)
            {
                printf("global node %" PRIDX "\tlocal_part[%" PRIDX "] = %" PRIDX "\n",
                index,
                index - index_start,
                mla_ctx->metis_mla[cnt_num_level].fine->local_part[index - index_start]);
            }

            puts("\n--------\n\n");
        }
    }
#endif // local partition value

#if 0
    if (my_rank == 0)
    {
        puts("\nglobal partition value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].fine->part[index]);
        }

        puts("\n========\n");

        puts("\ncoarse global partition value:");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].coarse->part[index]);
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

int MLASolverRelativeResidual(MySolver *mysolver, double *value)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    *value = r_norm_2 / b_norm_2;

    return 0;
}

int ParMetisMLASolverSolvePhase(const ConfigJSON *config,
                                MLAContext *mla_ctx,
                                int order_rbm,
                                MySolver *mysolver)
{
    return 0;
}

int ParMetisMLASolver(MLAContext *mla_ctx /*mla context data*/,
                      int mla_phase /*multilevel phase*/)
{
    PetscLogDouble time1, time2;

    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        PetscCall(PetscTime(&time1));
        ParMetisMLASolverSetupPhase(mla_ctx);
        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> setup time: %g (s)\n", time2 - time1));
    }

    // solve phase
    if (mla_phase == 1 || mla_phase == 2)
    {
        PetscCall(PetscTime(&time1));

        if (mla_ctx->setup == 0)
        {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase has not been done! failed\n"));
            exit(EXIT_FAILURE);
        }

        int iter_cnt = 0;
        double rela_resid = 0.;
        MLASolverRelativeResidual(&(mla_ctx->mysolver), &rela_resid);
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));

        while (iter_cnt < mla_ctx->config.mla_config.mla_max_it &&
               rela_resid > mla_ctx->config.mla_config.mla_rtol)
        {
            ParMetisMLASolverSolvePhase(&(mla_ctx->config),
                                        mla_ctx,
                                        mla_ctx->order_rbm,
                                        &(mla_ctx->mysolver));
            MLASolverRelativeResidual(&(mla_ctx->mysolver), &rela_resid);
            ++iter_cnt;
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));
        }

        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> solve time: %g (s)\n", time2 - time1));
    }

    return 0;
}
