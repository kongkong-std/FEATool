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

    // Calculate average coordinates for each coarse node (super node)
    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t id_part = fine_graph_data->part[index];
        for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
        {
            coarse_graph_data->coordinates[coarse_graph_data->dim * id_part + index_i] +=
                fine_graph_data->coordinates[fine_graph_data->dim * index + index_i];
        }
        ++(cnt_node_partition[id_part]);
    }

    // Average the coordinates
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        if (cnt_node_partition[index] > 0) // Check to avoid division by zero
        {
            for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
            {
                coarse_graph_data->coordinates[coarse_graph_data->dim * index + index_i] /= cnt_node_partition[index];
            }
        }
    }

    // Initialize adjacency list generator
    coarse_graph_data->xadj = (idx_t *)calloc(coarse_graph_data->nn + 1, sizeof(idx_t));
    assert(coarse_graph_data->xadj);

    // Create adjacency matrix for coarse graph
    bool **mat_adj = (bool **)malloc(coarse_graph_data->nn * sizeof(bool *));
    assert(mat_adj);

    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_graph_data->nn, sizeof(bool));
        assert(mat_adj[index]);
    }

    // Build adjacency matrix based on fine graph connections
    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t current_coarse_node = fine_graph_data->part[index];
        idx_t index_start = fine_graph_data->xadj[index];
        idx_t index_end = fine_graph_data->xadj[index + 1];

        // Check all neighbors of current fine node
        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            idx_t neighbor_fine_node = fine_graph_data->adjncy[index_i];
            idx_t neighbor_coarse_node = fine_graph_data->part[neighbor_fine_node];

            // If fine nodes belong to different coarse nodes, create edge between coarse nodes
            if (current_coarse_node != neighbor_coarse_node)
            {
                mat_adj[current_coarse_node][neighbor_coarse_node] = true;
                mat_adj[neighbor_coarse_node][current_coarse_node] = true;
            }
        }
    }

    // Ensure diagonal is false (no self-loops)
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index][index] = false;
    }

    // Build xadj array (cumulative adjacency counts)
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

    // Build adjncy array (actual adjacency list)
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

    // Initialize partition array for coarse graph
    coarse_graph_data->part = (idx_t *)malloc(coarse_graph_data->nn * sizeof(idx_t));
    assert(coarse_graph_data->part);

    // Free memory
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
    free(cnt_node_partition);

    return 0;
}

int DeepCopyDataMesh2Level0Fine(const AdjDataMesh *data_mesh /*mesh data*/,
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
    options[METIS_OPTION_SEED] = 42; // for any integer

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
    options[METIS_OPTION_SEED] = 42; // for any integer

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

int USAProlongation2SAProlongation(Mat A, Mat P_unsmooth, PetscScalar omega,
                                   Mat *P_smooth)
{
    PetscInt prolongation_m = 0, m_local_A, n_local_A;
    PetscCall(MatGetSize(A, &prolongation_m, NULL));
    PetscCall(MatGetLocalSize(A, &m_local_A, &n_local_A));

    Vec D_vec;
    PetscCall(MatCreateVecs(A, NULL, &D_vec));
    PetscCall(MatGetDiagonal(A, D_vec));

    PetscInt n_local;
    PetscScalar *d_array;
    PetscCall(VecGetLocalSize(D_vec, &n_local));
    PetscCall(VecGetArray(D_vec, &d_array));

    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, nprocs;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);
    printf("in rank %d/%d, m_local_A = %d, n_local = %d\n", my_rank, nprocs, m_local_A, n_local);

    for (PetscInt i = 0; i < n_local; i++)
    {
        if (PetscAbsScalar(d_array[i]) < PETSC_MACHINE_EPSILON)
        {
            d_array[i] = 1.0; // Set zero diagonal entries to 1
        }
    }

    PetscCall(VecRestoreArray(D_vec, &d_array));
    // Now compute D^{-1} (reciprocal of diagonal elements)
    PetscCall(VecReciprocal(D_vec));

    Mat D_inv;
    PetscCall(MatCreate(PETSC_COMM_WORLD, &D_inv));
    PetscCall(MatSetSizes(D_inv, m_local_A, n_local_A, prolongation_m, prolongation_m));
    PetscCall(MatSetType(D_inv, MATAIJ));
    PetscCall(MatSetUp(D_inv));
    PetscCall(MatDiagonalSet(D_inv, D_vec, INSERT_VALUES));
    PetscCall(MatAssemblyBegin(D_inv, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(D_inv, MAT_FINAL_ASSEMBLY));

    Mat Dinv_A;
    PetscCall(MatMatMult(D_inv, A, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &Dinv_A));
    PetscCall(MatScale(Dinv_A, -omega));

    PetscCall(MatShift(Dinv_A, 1.));

    // Mat A_times_usa;
    // PetscCall(MatMatMult(A, P_unsmooth, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &A_times_usa));

    // Mat Dinv_A_usa;
    // PetscCall(MatMatMult(D_inv, A_times_usa, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &Dinv_A_usa));

    // PetscCall(MatScale(Dinv_A_usa, omega));

    // PetscCall(MatDuplicate(P_unsmooth,
    //                        MAT_COPY_VALUES,
    //                        &P_smooth));
    // PetscCall(MatAXPY(P_unsmooth, -1.0, Dinv_A_usa, UNKNOWN_NONZERO_PATTERN));

    // PetscCall(MatDuplicate(P_unsmooth,
    //                        MAT_COPY_VALUES,
    //                        &P_smooth));

    PetscCall(MatMatMult(Dinv_A, P_unsmooth, MAT_INITIAL_MATRIX, PETSC_DETERMINE, P_smooth));

#if 1
    int m_dinv_a, n_dinv_a,
        m_p_unsmooth, n_p_unsmooth,
        m_p_smooth, n_p_smooth;
    MatGetSize(Dinv_A, &m_dinv_a, &n_dinv_a);
    MatGetSize(P_unsmooth, &m_p_unsmooth, &n_p_unsmooth);
    MatGetSize(*P_smooth, &m_p_smooth, &n_p_smooth);
    printf(">>>> size of dinv_a: m = %d, n = %d\n", m_dinv_a, n_dinv_a);
    printf(">>>> size of p_unsmooth: m = %d, n = %d\n", m_p_unsmooth, n_p_unsmooth);
    printf(">>>> size of p_smooth: m = %d, n = %d\n", m_p_smooth, n_p_smooth);
#endif

    PetscCall(MatDestroy(&Dinv_A));
    PetscCall(VecDestroy(&D_vec));
    PetscCall(MatDestroy(&D_inv));
    // PetscCall(MatDestroy(&A_times_usa));
    // PetscCall(MatDestroy(&Dinv_A_usa));

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

#if 0
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      6 * local_prolongation_m, 6 * local_prolongation_n,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].prolongation, MATAIJ)); // prolongation set type- mataij
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));
#endif // original implementation, usa-prolongation operator

                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].usa_prolongation))); // usa-prolongation operator
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                      6 * local_prolongation_m, 6 * local_prolongation_n,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

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
                            PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                  6 * index_fine_node + index_tmp_i,
                                                  6 * index_coarse_node + index_tmp_j,
                                                  p_loc[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));

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

#if 0
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // unsmoothed aggregation-based mg

#if 1
                // smoothed prolongation operator
                // PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].sa_prolongation))); // usa-prolongation operator
                // PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                       6 * local_prolongation_m, 6 * local_prolongation_n,
                //                       prolongation_m, prolongation_n));
                // PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, MATAIJ)); // prolongation set type- mataij
                // PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

                PetscScalar omega = 0.67;

                PetscCall(USAProlongation2SAProlongation(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                                         mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                         omega,
                                                         &mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

#if 1
                PetscInt m_sa_prolongation, n_sa_prolongation,
                    m_operator_fine, n_operator_fine;
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, &m_sa_prolongation, &n_sa_prolongation));
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &m_operator_fine, &n_operator_fine));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of sa_prolongation: m = %d, n = %d\n", m_sa_prolongation, n_sa_prolongation));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of operator_fine: m = %d, n = %d\n", m_operator_fine, n_operator_fine));
#endif // size

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // smoothed aggregation-based mg
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

                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].usa_prolongation)));
                if (index_cnt_level == 0)
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                          local_prolongation_m * 6, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                else
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                          local_prolongation_m * 9, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

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
                        p_loc[1][3] = coarse_node_z - fine_node_z; // -z
                        p_loc[1][5] = fine_node_x - coarse_node_x; // x
                        p_loc[2][3] = fine_node_y - coarse_node_y; // y
                        p_loc[2][4] = coarse_node_x - fine_node_x; // -x

#if 1
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2

                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2

                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
#endif                                                                    // kirchhoff 2nd order
#if 1
                        p_loc[3][7] = coarse_node_y - fine_node_y;        // -y
                        p_loc[3][8] = (coarse_node_x - fine_node_x) / 2.; // -x/2
                        p_loc[4][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[4][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
#endif                                                                    // suppose 0 in right-bottom

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
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
                        p_loc[1][3] = coarse_node_z - fine_node_z; // -z
                        p_loc[1][5] = fine_node_x - coarse_node_x; // x
                        p_loc[2][3] = fine_node_y - coarse_node_y; // y
                        p_loc[2][4] = coarse_node_x - fine_node_x; // -x

#if 1
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2

                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2

                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
#endif                                                                    // kirchhoff 2nd-order
#if 1
                        p_loc[3][7] = coarse_node_y - fine_node_y;        // -y
                        p_loc[3][8] = (coarse_node_x - fine_node_x) / 2.; // -x/2
                        p_loc[4][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[4][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
#endif                                                                    // suppose 0 in right-bottom

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                      9 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));

#if 0
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // unsmoothed aggregation-based mg

#if 1
                // smoothed prolongation operator
                // PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].sa_prolongation))); // usa-prolongation operator
                // if (index_cnt_level == 0)
                // {
                //     PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                           local_prolongation_m * 6, local_prolongation_n * 9,
                //                           prolongation_m, prolongation_n));
                // }
                // else
                // {
                //     PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                           local_prolongation_m * 9, local_prolongation_n * 9,
                //                           prolongation_m, prolongation_n));
                // }
                // PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                // PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

                PetscScalar omega = 0.67;

                PetscCall(USAProlongation2SAProlongation(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                                         mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                         omega,
                                                         &mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

#if 1
                PetscInt m_sa_prolongation, n_sa_prolongation,
                    m_operator_fine, n_operator_fine;
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, &m_sa_prolongation, &n_sa_prolongation));
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &m_operator_fine, &n_operator_fine));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of sa_prolongation: m = %d, n = %d\n", m_sa_prolongation, n_sa_prolongation));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of operator_fine: m = %d, n = %d\n", m_operator_fine, n_operator_fine));
#endif // size

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // smoothed aggregation-based mg
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
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].usa_prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                      6 * local_prolongation_m, 6 * local_prolongation_n,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

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
                            PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                  6 * index_fine_node + index_tmp_i,
                                                  6 * index_coarse_node + index_tmp_j,
                                                  p_loc[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));

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

#if 0
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // unsmoothed aggregation-based mg

#if 1
                // smoothed prolongation operator
                // PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].sa_prolongation))); // usa-prolongation operator
                // PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                       6 * local_prolongation_m, 6 * local_prolongation_n,
                //                       prolongation_m, prolongation_n));
                // PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, MATAIJ)); // prolongation set type- mataij
                // PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

                PetscScalar omega = 0.67;

                PetscCall(USAProlongation2SAProlongation(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                                         mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                         omega,
                                                         &mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

#if 1
                PetscInt m_sa_prolongation, n_sa_prolongation,
                    m_operator_fine, n_operator_fine;
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, &m_sa_prolongation, &n_sa_prolongation));
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &m_operator_fine, &n_operator_fine));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of sa_prolongation: m = %d, n = %d\n", m_sa_prolongation, n_sa_prolongation));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of operator_fine: m = %d, n = %d\n", m_operator_fine, n_operator_fine));
#endif // size

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // smoothed aggregation-based mg
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

                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].usa_prolongation)));
                if (index_cnt_level == 0)
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                          local_prolongation_m * 6, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                else
                {
                    PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                          local_prolongation_m * 9, local_prolongation_n * 9,
                                          prolongation_m, prolongation_n));
                }
                PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

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
#if 1
                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[1][4] = fine_node_z - coarse_node_z; // z
                        p_loc[1][5] = fine_node_x - coarse_node_x; // x
                        p_loc[2][3] = coarse_node_x - fine_node_x; // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y; // -y
#if 1
                        p_loc[0][7] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[0][8] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[1][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[1][8] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // -xz
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
#endif                                                               // nullspace 2nd-order

#if 0
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2

                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2

                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
#endif                                                                    // kirchhoff 2nd-order
#if 0
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2
#endif                                                                    // suppose 0 in right-bottom
#endif                                                                    // kirchhoff-love

#if 0
                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][7] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[0][8] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[1][4] = fine_node_z - coarse_node_z;        // z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[1][8] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
                        p_loc[2][3] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y;        // -y
                        p_loc[2][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
#if 0
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2
#endif // suppose 0 in right-bottom
#endif // nullspace

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
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

#if 1
                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[1][4] = fine_node_z - coarse_node_z; // z
                        p_loc[1][5] = fine_node_x - coarse_node_x; // x
                        p_loc[2][3] = coarse_node_x - fine_node_x; // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y; // -y
#if 1
                        p_loc[0][7] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[0][8] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[1][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[1][8] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // -xz
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
#endif                                                               // nullspace 2nd-order

#if 0
                        p_loc[0][6] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x); // zx
                        p_loc[0][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y) / 2.; // zy/2

                        p_loc[1][7] = (fine_node_z - coarse_node_z) *
                                      (fine_node_y - coarse_node_y); // zy
                        p_loc[1][8] = (fine_node_z - coarse_node_z) *
                                      (fine_node_x - coarse_node_x) / 2.; // zx/2

                        p_loc[2][6] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_x - coarse_node_x) / 2.; // -x^2/2
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_y - coarse_node_y) / 2.; // -y^2/2
                        p_loc[2][8] = -(fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y) / 2.; // -xy/2
#endif // kirchhoff 2nd-order
#if 0
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2
#endif // suppose 0 in right-bottom
#endif // krichhoff

#if 0
                        p_loc[0][3] = fine_node_z - coarse_node_z; // z
                        p_loc[0][5] = coarse_node_y - fine_node_y; // -y
                        p_loc[0][7] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[0][8] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[1][4] = fine_node_z - coarse_node_z;        // z
                        p_loc[1][5] = fine_node_x - coarse_node_x;        // x
                        p_loc[1][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_y - coarse_node_y); // xy
                        p_loc[1][8] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
                        p_loc[2][3] = coarse_node_x - fine_node_x;        // -x
                        p_loc[2][4] = coarse_node_y - fine_node_y;        // -y
                        p_loc[2][6] = (fine_node_x - coarse_node_x) *
                                      (fine_node_z - coarse_node_z); // xz
                        p_loc[2][7] = -(fine_node_y - coarse_node_y) *
                                      (fine_node_z - coarse_node_z); // -yz
#if 0
                        p_loc[3][6] = fine_node_x - coarse_node_x;        // x
                        p_loc[3][8] = (fine_node_y - coarse_node_y) / 2.; // y/2
                        p_loc[4][7] = fine_node_y - coarse_node_y;        // y
                        p_loc[4][8] = (fine_node_x - coarse_node_x) / 2.; // x/2
#endif // suppose 0 in right-bottom
#endif // nullspace

                        for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                        {
                            for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                            {
                                PetscCall(MatSetValue(mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                      9 * index_fine_node + index_tmp_i,
                                                      9 * index_coarse_node + index_tmp_j,
                                                      p_loc[index_tmp_i][index_tmp_j],
                                                      INSERT_VALUES));
                            }
                        }
                    }
                }

                PetscCall(MatAssemblyBegin(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MAT_FINAL_ASSEMBLY));

#if 0
                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // unsmoothed aggregation-based mg

#if 1
                // smoothed prolongation operator
                // PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].sa_prolongation))); // usa-prolongation operator
                // if (index_cnt_level == 0)
                // {
                //     PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                           local_prolongation_m * 6, local_prolongation_n * 9,
                //                           prolongation_m, prolongation_n));
                // }
                // else
                // {
                //     PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                //                           local_prolongation_m * 9, local_prolongation_n * 9,
                //                           prolongation_m, prolongation_n));
                // }
                // PetscCall(MatSetType(mla_ctx->metis_mla[index_cnt_level].usa_prolongation, MATAIJ)); // prolongation set type- mataij
                // PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].usa_prolongation));

                PetscScalar omega = 0.67;

                PetscCall(USAProlongation2SAProlongation(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                                         mla_ctx->metis_mla[index_cnt_level].usa_prolongation,
                                                         omega,
                                                         &mla_ctx->metis_mla[index_cnt_level].sa_prolongation));

#if 1
                PetscInt m_sa_prolongation, n_sa_prolongation,
                    m_operator_fine, n_operator_fine;
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].sa_prolongation, &m_sa_prolongation, &n_sa_prolongation));
                PetscCall(MatGetSize(mla_ctx->metis_mla[index_cnt_level].operator_fine, &m_operator_fine, &n_operator_fine));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of sa_prolongation: m = %d, n = %d\n", m_sa_prolongation, n_sa_prolongation));
                PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> size of operator_fine: m = %d, n = %d\n", m_operator_fine, n_operator_fine));
#endif // size

                // coarse operator
                PetscCall(MatPtAP(mla_ctx->metis_mla[index_cnt_level].operator_fine,
                                  mla_ctx->metis_mla[index_cnt_level].sa_prolongation,
                                  MAT_INITIAL_MATRIX,
                                  PETSC_DETERMINE,
                                  &(mla_ctx->metis_mla[index_cnt_level].operator_coarse)));

                // next level fine operator
                PetscCall(MatDuplicate(mla_ctx->metis_mla[index_cnt_level].operator_coarse,
                                       MAT_COPY_VALUES,
                                       &(mla_ctx->metis_mla[index_cnt_level + 1].operator_fine)));
#endif // smoothed aggregation-based mg
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

    for (int level = 0; level < cnt_num_level; ++level)
    {
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_presmooth)));
#if 0
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_presmooth,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));
#endif

        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_postsmooth)));
#if 0
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_postsmooth,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));
#endif

        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_coarse)));
#if 0
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_coarse,
                                  mla_ctx->metis_mla[level].operator_coarse,
                                  mla_ctx->metis_mla[level].operator_coarse));
#endif
    }

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

int ParMetisMLANestedProcedurePreSmooth(KSP ksp, PC pc,
                                        int level,
                                        MLAContext *mla_ctx,
                                        Vec *mg_recur_x,
                                        Vec *mg_recur_b,
                                        int v_pre_smooth)
{
#if 0
    // Duplicate vector if not at finest level
    if (level != 0)
    {
        PetscCall(VecDuplicate(mg_recur_b[level], mg_recur_x + level));
    }

    // Reset KSP to clear previous configurations
    // PetscCall(KSPReset(ksp));

    if (level == 0)
    {
        // Finest level - use ParaSails smoother
        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));

        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCHYPRE));
        PetscCall(PCHYPRESetType(pc, "parasails"));

        PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v_pre_smooth));
        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));

        PetscCall(KSPSetUp(ksp));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }
    else
    {
        // Coarser levels - use ParaSails with matrix shift
        double shift = 1e-12; // regularization for stability
        PetscCall(MatShift(mla_ctx->metis_mla[level].operator_fine, shift));

        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));

        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCHYPRE));
        PetscCall(PCHYPRESetType(pc, "parasails"));

        PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v_pre_smooth));
        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));

        PetscCall(KSPSetUp(ksp));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }

    return 0;
#endif // smoother-hypre

#if 0
    // PetscCall(VecDuplicate(mg_recur_b[level], mg_recur_x + level));
    if (level != 0)
    {
        PetscCall(VecDuplicate(mg_recur_b[level], mg_recur_x + level));
    }

    if (level == 0)
    {
#if 0
    // KSP ksp_loc;
    // PC pc_loc;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
#endif
        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));
        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCSOR));
        PetscCall(PCSORSetOmega(pc, 1.)); // gauss-seidel
        PetscCall(PCSORSetIterations(pc, v_pre_smooth, v_pre_smooth));
        PetscCall(PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
        // PetscCall(KSPSetFromOptions(ksp));

        PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }
    else
    {
#if 1
        double shift = 1e-12; // test metis aggregation implementation
        // double shift = 0.; // test metis aggregation implementation
        PetscCall(MatShift(mla_ctx->metis_mla[level].operator_fine, shift));
#endif
        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine)); // add shift to diagonal
        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCSOR));
        PetscCall(PCSORSetOmega(pc, 1.)); // gauss-seidel
        PetscCall(PCSORSetIterations(pc, v_pre_smooth, v_pre_smooth));
        PetscCall(PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
        // PetscCall(KSPSetFromOptions(ksp));

        PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }

#if 0
    PetscCall(KSPDestroy(&ksp_loc));
    PetscCall(PCDestroy(&pc_loc));
#endif

    return 0;
#endif // smoother-sor

#if 1
    // PetscCall(VecDuplicate(mg_recur_b[level], mg_recur_x + level));
    if (level != 0)
    {
        PetscCall(VecDuplicate(mg_recur_b[level], mg_recur_x + level));
    }

    if (level == 0)
    {
        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));
        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));

        // Configure ASM
        PetscCall(PCSetType(pc, PCASM));
        PetscCall(PCASMSetOverlap(pc, 2));

        // Set up the preconditioner first
        PetscCall(PCSetUp(pc));

        // Get subdomain KSP contexts
        KSP *subksp;
        PetscInt nlocal, first;
        PetscCall(PCASMGetSubKSP(pc, &nlocal, &first, &subksp));

        // Configure SOR on each subdomain
        for (PetscInt i = 0; i < nlocal; i++)
        {
            PC subpc;
            PetscCall(KSPGetPC(subksp[i], &subpc));
            PetscCall(PCSetType(subpc, PCSOR));
            PetscCall(PCSORSetOmega(subpc, 1.0)); // Gauss-Seidel
            PetscCall(PCSORSetIterations(subpc, v_pre_smooth, v_pre_smooth));
            PetscCall(PCSORSetSymmetric(subpc, SOR_SYMMETRIC_SWEEP));

            // Set subdomain KSP parameters
            PetscCall(KSPSetType(subksp[i], KSPRICHARDSON));
            PetscCall(KSPSetTolerances(subksp[i], 1e-10, 1e-10, PETSC_DEFAULT, 1));
            PetscCall(KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE));
            PetscCall(KSPSetNormType(subksp[i], KSP_NORM_UNPRECONDITIONED));
        }

        // Set main KSP parameters
        PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }
    else
    {
        double shift = 1e-12; // test metis aggregation implementation
        PetscCall(MatShift(mla_ctx->metis_mla[level].operator_fine, shift));

        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));
        PetscCall(KSPSetType(ksp, KSPRICHARDSON));
        PetscCall(KSPGetPC(ksp, &pc));

        // Configure ASM
        PetscCall(PCSetType(pc, PCASM));
        PetscCall(PCASMSetOverlap(pc, 2));

        // Set up the preconditioner first
        PetscCall(PCSetUp(pc));

        // Get subdomain KSP contexts
        KSP *subksp;
        PetscInt nlocal, first;
        PetscCall(PCASMGetSubKSP(pc, &nlocal, &first, &subksp));

        // Configure SOR on each subdomain
        for (PetscInt i = 0; i < nlocal; i++)
        {
            PC subpc;
            PetscCall(KSPGetPC(subksp[i], &subpc));
            PetscCall(PCSetType(subpc, PCSOR));
            PetscCall(PCSORSetOmega(subpc, 1.0)); // Gauss-Seidel
            PetscCall(PCSORSetIterations(subpc, v_pre_smooth, v_pre_smooth));
            PetscCall(PCSORSetSymmetric(subpc, SOR_SYMMETRIC_SWEEP));

            // Set subdomain KSP parameters
            PetscCall(KSPSetType(subksp[i], KSPRICHARDSON));
            PetscCall(KSPSetTolerances(subksp[i], 1e-10, 1e-10, PETSC_DEFAULT, 1));
            PetscCall(KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE));
            PetscCall(KSPSetNormType(subksp[i], KSP_NORM_UNPRECONDITIONED));
        }

        // Set main KSP parameters
        PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
        PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));
    }

    return 0;
#endif // smoother-asm
}

int ParMetisMLASolverCoarsetCorrectionPhase(int order_rbm, KSP ksp, PC pc,
                                            int level,
                                            MLAContext *mla_ctx,
                                            Vec *mg_recur_x,
                                            Vec *mg_recur_b)
{
    PetscCall(VecDuplicate(mg_recur_b[level + 1], mg_recur_x + level + 1));
    int gcr_restart = mla_ctx->config.mla_config.coarse_restart;

#if 0
    PetscCall(VecView(mg_recur_b[level + 1], PETSC_VIEWER_STDOUT_WORLD));
#endif // residual information

#if 0
    printf(">>>> in coarset level, level = %d\n", level);
    int vec_size = 0, m_a_H = 0, n_a_H = 0;
    PetscCall(MatGetSize((mla_ctx->mla + level)->operator_coarse, &m_a_H, &n_a_H));
    PetscCall(VecGetSize(mg_recur_b[level + 1], &vec_size));
    printf("coarsest matrix size: m = %d, n = %d\n", m_a_H, n_a_H);
    printf("coarsest vector size: m = %d\n", vec_size);
#endif // size information

#if 0
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
#endif

#if 1
    if (order_rbm == 1)
    {
        if (level == 0)
        {
            // rbm order 1
            PetscCall(KSPSetOperators(ksp,
                                      mla_ctx->metis_mla[level].operator_coarse,
                                      mla_ctx->metis_mla[level].operator_coarse));
            PetscCall(KSPSetType(ksp, KSPGCR));
            // PetscCall(KSPGetPC(ksp_H, &pc));
            // PetscCall(PCSetType(pc, PCLU));
            PetscCall(KSPGCRSetRestart(ksp, gcr_restart));
            PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
            // PetscCall(KSPSetFromOptions(ksp));
            PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1000));
            // PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 10));

            PetscCall(KSPSolve(ksp, mg_recur_b[level + 1], mg_recur_x[level + 1]));
        }
        else
        {
#if 1
            // double shift = 1e-12; // test metis aggregation implementation
            double shift = 0.; // test metis aggregation implementation
            PetscCall(MatShift(mla_ctx->metis_mla[level].operator_coarse, shift));
#endif
            PetscCall(KSPSetOperators(ksp,
                                      mla_ctx->metis_mla[level].operator_coarse,
                                      mla_ctx->metis_mla[level].operator_coarse)); // add shift to diagonal
            PetscCall(KSPSetType(ksp, KSPGCR));
#if 0
        PetscCall(KSPGetPC(ksp, pc));
        PetscCall(PCSetType(pc, PCSVD));
#endif
            PetscCall(KSPGCRSetRestart(ksp, gcr_restart));
            PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
            // PetscCall(KSPSetFromOptions(ksp));
            PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1000));
            // PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 10));

            PetscCall(KSPSolve(ksp, mg_recur_b[level + 1], mg_recur_x[level + 1]));
        }
    }
    else if (order_rbm == 2)
    {
#if 1
        // double shift = 1e-12; // test metis aggregation implementation
        double shift = 0.; // test metis aggregation implementation
        PetscCall(MatShift(mla_ctx->metis_mla[level].operator_coarse, shift));
#endif
        PetscCall(KSPSetOperators(ksp,
                                  mla_ctx->metis_mla[level].operator_coarse,
                                  mla_ctx->metis_mla[level].operator_coarse)); // add shift to diagonal
        PetscCall(KSPSetType(ksp, KSPGCR));
#if 0
        PetscCall(KSPGetPC(ksp, pc));
        PetscCall(PCSetType(pc, PCSVD));
#endif
        PetscCall(KSPGCRSetRestart(ksp, gcr_restart));
        PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
        // PetscCall(KSPSetFromOptions(ksp));
        PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1000));
        // PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 10));

        PetscCall(KSPSolve(ksp, mg_recur_b[level + 1], mg_recur_x[level + 1]));
    }
#endif

#if 0
    double norm_r_H = 0.;
    PetscCall(VecNorm(mg_recur_b[level + 1], NORM_2, &norm_r_H));
    Vec tmp_r_H;
    PetscCall(VecDuplicate(mg_recur_b[level + 1], &tmp_r_H));
    PetscCall(MatMult((mla_ctx->mla + level)->operator_coarse, mg_recur_x[level + 1], tmp_r_H));
    PetscCall(VecAXPY(tmp_r_H, -1., mg_recur_b[level + 1]));
    double norm_tmp_r_H = 0.;
    PetscCall(VecNorm(tmp_r_H, NORM_2, &norm_tmp_r_H));

    double norm_e_H = 0.;
    PetscCall(VecNorm(mg_recur_x[level + 1], NORM_2, &norm_e_H));

    printf(">>>>>>>> coarse correction: norm_r_H = %021.16le\n", norm_r_H);
    printf(">>>>>>>> coarse correction: norm_e_H = %021.16le\n", norm_e_H);
    printf(">>>>>>>> coarse correction: relative = %021.16le\n\n", norm_tmp_r_H / norm_r_H);
#endif // print coarse level correction information

    return 0;
}

int ParMetisMLANestedProcedurePostSmooth(KSP ksp, PC pc,
                                         int level,
                                         MLAContext *mla_ctx,
                                         Vec *mg_recur_x,
                                         Vec *mg_recur_b,
                                         int v_post_smooth)
{
#if 0
    // Alternative 1: Use HYPRE Jacobi smoother (simpler)
    // PetscCall(KSPReset(ksp));
    PetscCall(KSPSetOperators(ksp,
                              mla_ctx->metis_mla[level].operator_fine,
                              mla_ctx->metis_mla[level].operator_fine));

    PetscCall(KSPSetType(ksp, KSPRICHARDSON));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCHYPRE));
    PetscCall(PCHYPRESetType(pc, "pilut")); // or "pilut"

    PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, v_post_smooth));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
    PetscCall(KSPSetUp(ksp));
    PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));

    return 0;
#endif // smoother-hypre

#if 0
#if 0
    // KSP ksp_loc;
    // PC pc_loc;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
#endif
    PetscCall(KSPSetOperators(ksp,
                              mla_ctx->metis_mla[level].operator_fine,
                              mla_ctx->metis_mla[level].operator_fine));
    PetscCall(KSPSetType(ksp, KSPRICHARDSON));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCSOR));
    PetscCall(PCSORSetOmega(pc, 1.)); // gauss-seidel
    PetscCall(PCSORSetIterations(pc, v_post_smooth, v_post_smooth));
    PetscCall(PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP));
    PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
    // PetscCall(KSPSetFromOptions(ksp));

    PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
    PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));

    return 0;
#endif // smoother-sor

#if 1
    PetscCall(KSPSetOperators(ksp,
                              mla_ctx->metis_mla[level].operator_fine,
                              mla_ctx->metis_mla[level].operator_fine));
    PetscCall(KSPSetType(ksp, KSPRICHARDSON));
    PetscCall(KSPGetPC(ksp, &pc));

    // Configure ASM
    PetscCall(PCSetType(pc, PCASM));
    PetscCall(PCASMSetOverlap(pc, 2));

    // Set up the preconditioner first
    PetscCall(PCSetUp(pc));

    // Get subdomain KSP contexts
    KSP *subksp;
    PetscInt nlocal, first;
    PetscCall(PCASMGetSubKSP(pc, &nlocal, &first, &subksp));

    // Configure SOR on each subdomain
    for (PetscInt i = 0; i < nlocal; i++)
    {
        PC subpc;
        PetscCall(KSPGetPC(subksp[i], &subpc));
        PetscCall(PCSetType(subpc, PCSOR));
        PetscCall(PCSORSetOmega(subpc, 1.0)); // Gauss-Seidel
        PetscCall(PCSORSetIterations(subpc, v_post_smooth, v_post_smooth));
        PetscCall(PCSORSetSymmetric(subpc, SOR_SYMMETRIC_SWEEP));

        // Set subdomain KSP parameters
        PetscCall(KSPSetType(subksp[i], KSPRICHARDSON));
        PetscCall(KSPSetTolerances(subksp[i], 1e-10, 1e-10, PETSC_DEFAULT, 1));
        PetscCall(KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE));
        PetscCall(KSPSetNormType(subksp[i], KSP_NORM_UNPRECONDITIONED));
    }

    // Set main KSP parameters
    PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));
    PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
    PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));

    return 0;
#endif // smoother-asm
}

int ParMetisMLANestedProcedure(int level /*level*/,
                               int num_level /*number of levels*/,
                               MySolver *mysolver /*solver data*/,
                               MLAContext *mla_ctx /*mla context*/,
                               Vec *mg_recur_x /*x*/,
                               Vec *mg_recur_b /*b*/,
                               int v_pre_smooth /*pre-smooth times*/,
                               int v_post_smooth /*post-smooth times*/,
                               int order_rbm /*rbm order*/)
{
    Vec *r_h = NULL, *tmp_e_h = NULL;

    r_h = (Vec *)malloc(num_level * sizeof(Vec));
    tmp_e_h = (Vec *)malloc(num_level * sizeof(Vec));
    assert(r_h && tmp_e_h);

    // loop implementation
    /*
     * v-cycle downward direction, from fine mesh to coarse mesh
     */
    for (level = 0; level < num_level; ++level)
    {
        PetscCall(VecDuplicate(mg_recur_b[level], r_h + level));

        // pre-smooth procedure
        ParMetisMLANestedProcedurePreSmooth(mla_ctx->metis_mla[level].ksp_presmooth,
                                            mla_ctx->metis_mla[level].pc_presmooth,
                                            level,
                                            mla_ctx,
                                            mg_recur_x,
                                            mg_recur_b,
                                            v_pre_smooth);
        PetscCall(MatMult(mla_ctx->metis_mla[level].operator_fine,
                          mg_recur_x[level],
                          r_h[level]));
        PetscCall(VecAYPX(r_h[level], -1., mg_recur_b[level]));

        // restriction
        PetscInt m_prolongation = 0, n_prolongation = 0; // size of prolongation operator
        PetscCall(MatGetSize(mla_ctx->metis_mla[level].prolongation,
                             &m_prolongation,
                             &n_prolongation));
        PetscInt local_n_prolongation = 0;
        PetscCall(MatGetLocalSize(mla_ctx->metis_mla[level].prolongation, NULL, &local_n_prolongation));
        PetscCall(VecCreate(PETSC_COMM_WORLD, mg_recur_b + level + 1));
        PetscCall(VecSetSizes(mg_recur_b[level + 1], local_n_prolongation, n_prolongation));
        PetscCall(VecSetFromOptions(mg_recur_b[level + 1]));
        PetscCall(MatMultTranspose(mla_ctx->metis_mla[level].prolongation, r_h[level], mg_recur_b[level + 1]));
    }

    // coarsest level
    ParMetisMLASolverCoarsetCorrectionPhase(order_rbm,
                                            mla_ctx->metis_mla[level - 1].ksp_coarse,
                                            mla_ctx->metis_mla[level - 1].pc_coarse,
                                            level - 1,
                                            mla_ctx,
                                            mg_recur_x,
                                            mg_recur_b);

    /*
     * v-cycle upward direction, from coarse mesh to fine mesh
     */
    for (level = num_level - 1; level >= 0; --level)
    {
        PetscCall(VecDuplicate(mg_recur_x[level], tmp_e_h + level));

        PetscCall(MatMult(mla_ctx->metis_mla[level].prolongation, mg_recur_x[level + 1], tmp_e_h[level]));
        PetscCall(VecAXPY(mg_recur_x[level], 1., tmp_e_h[level]));

        // post-smooth procedure
        ParMetisMLANestedProcedurePostSmooth(mla_ctx->metis_mla[level].ksp_postsmooth,
                                             mla_ctx->metis_mla[level].pc_postsmooth,
                                             level,
                                             mla_ctx,
                                             mg_recur_x,
                                             mg_recur_b,
                                             v_post_smooth);
    }

    // free memory
    for (int index = 0; index < num_level; ++index)
    {
        PetscCall(VecDestroy(r_h + index));
        PetscCall(VecDestroy(tmp_e_h + index));
    }
    free(r_h);
    free(tmp_e_h);

    return 0;
}

int ParMetisMLASolverSolvePhase(const ConfigJSON *config,
                                MLAContext *mla_ctx,
                                int order_rbm,
                                MySolver *mysolver)
{
    // mg recursive implementation
    int v_pre_smooth = config->mla_config.pre_smooth_v;
    int v_post_smooth = config->mla_config.post_smooth_v;
    int num_level = mla_ctx->num_level;
    Vec *mg_recur_x, *mg_recur_b;

#if 1
    mg_recur_x = (Vec *)malloc((num_level + 1) * sizeof(Vec));
    mg_recur_b = (Vec *)malloc((num_level + 1) * sizeof(Vec));
    assert(mg_recur_x && mg_recur_b);
#endif

    PetscCall(VecDuplicate(mysolver->solver_x, &mg_recur_x[0]));
    PetscCall(VecDuplicate(mysolver->solver_b, &mg_recur_b[0]));
    PetscCall(VecCopy(mysolver->solver_b, mg_recur_b[0]));
    PetscCall(VecCopy(mysolver->solver_x, mg_recur_x[0]));

    ParMetisMLANestedProcedure(0, num_level,
                               mysolver,
                               mla_ctx,
                               mg_recur_x, mg_recur_b,
                               v_pre_smooth, v_post_smooth,
                               order_rbm);

    // updating solution after nested mg procedure
    PetscCall(VecCopy(mg_recur_x[0], mysolver->solver_x));

// free memeory
#if 1
    for (int index = 0; index < num_level + 1; ++index)
    // for (int index = 0; index < num_level; ++index)
    {
        PetscCall(VecDestroy(mg_recur_x + index));
        PetscCall(VecDestroy(mg_recur_b + index));
    }
#endif
    free(mg_recur_x);
    free(mg_recur_b);

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
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> mla setup time: %g (s)\n", time2 - time1));
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
        // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));

        while (iter_cnt < mla_ctx->config.mla_config.mla_max_it &&
               rela_resid > mla_ctx->config.mla_config.mla_rtol)
        {
            ParMetisMLASolverSolvePhase(&(mla_ctx->config),
                                        mla_ctx,
                                        mla_ctx->order_rbm,
                                        &(mla_ctx->mysolver));
            MLASolverRelativeResidual(&(mla_ctx->mysolver), &rela_resid);
            ++iter_cnt;
            // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));
        }

        PetscCall(PetscTime(&time2));
        // PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> mla solve time: %g (s)\n", time2 - time1));
    }

    return 0;
}
