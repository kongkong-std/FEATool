#include "../include/main.h"

#define MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))

int SAMGFineVertex2PartitionMap(AggData **agg /*aggregation data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    AggData *data_agg = *agg;

    int *fgid2part = (int *)malloc(data_agg->nv * sizeof(int));
    data_agg->fgid2part = (int *)malloc(data_agg->nv * sizeof(int));
    assert(fgid2part && data_agg->fgid2part);

    memset(fgid2part, -1, data_agg->nv * sizeof(int));

    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        if (data_agg->owner[index_p] == my_rank)
        {
            for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
            {
                int fine_vertex = data_agg->fine_gids[index_p][index_v];
                fgid2part[fine_vertex] = index_p;
            }
        }
    }

    (void)MPI_Allreduce(fgid2part,
                        data_agg->fgid2part,
                        data_agg->nv,
                        MPI_INT,
                        MPI_MAX,
                        comm);

    // free memory
    free(fgid2part);

    return 0;
}

int SAMGCoarseAdjacentListConstructor(AggData **agg,
                                      MeshData **mesh_f,
                                      MeshData **mesh_c)
{
    int my_rank, nprocs;
    MPI_Comm comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    MeshData *data_mesh_c = *mesh_c;
    AggData *data_agg = *agg;

    // Initialize coarse mesh adjacency structure
    data_mesh_c->data_adj.nv = data_mesh_c->nv;
    data_mesh_c->data_adj.local_nv = data_mesh_c->vtxdist[my_rank + 1] - data_mesh_c->vtxdist[my_rank];
    data_mesh_c->data_adj.idx = (int *)malloc(data_mesh_c->data_adj.local_nv * sizeof(int));
    data_mesh_c->data_adj.xadj = (int *)calloc(data_mesh_c->data_adj.local_nv + 1, sizeof(int));
    assert(data_mesh_c->data_adj.idx && data_mesh_c->data_adj.xadj);

    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        data_mesh_c->data_adj.idx[index] = index + data_mesh_c->vtxdist[my_rank];
    }

    int coarse_index_start = data_mesh_c->vtxdist[my_rank];
    int coarse_index_end = data_mesh_c->vtxdist[my_rank + 1];

    // Allocate marking arrays - one per coarse vertex to track unique adjacencies
    int **mark = (int **)malloc(data_mesh_c->data_adj.local_nv * sizeof(int *));
    assert(mark);
    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        mark[index] = (int *)calloc(data_agg->np, sizeof(int));
        assert(mark[index]);
    }

    // Build adjacency by iterating over fine vertices
    for (int index_v = 0; index_v < data_mesh_f->data_adj.local_nv; ++index_v)
    {
        int fine_v = data_mesh_f->data_adj.idx[index_v];
        int partition_v = data_agg->fgid2part[fine_v];

        // Skip if this fine vertex doesn't belong to a local coarse vertex
        if (partition_v < coarse_index_start || partition_v >= coarse_index_end)
        {
            continue;
        }

        int partition_local = partition_v - coarse_index_start;

        // Iterate over neighbors of fine vertex
        int index_v_start = data_mesh_f->data_adj.xadj[index_v];
        int index_v_end = data_mesh_f->data_adj.xadj[index_v + 1];

        for (int index_k = index_v_start; index_k < index_v_end; ++index_k)
        {
            int fine_u = data_mesh_f->data_adj.adjncy[index_k];
            int partition_u = data_agg->fgid2part[fine_u];

            // Mark connection to different partition (coarse edge)
            if (partition_u != partition_v)
            {
                mark[partition_local][partition_u] = 1;
            }
        }
    }

    // Count unique adjacencies and build adjacency lists
    int *cnt_nv_coarse = (int *)calloc(data_mesh_c->data_adj.local_nv, sizeof(int));
    int **adj_coarse = (int **)malloc(data_mesh_c->data_adj.local_nv * sizeof(int *));
    assert(cnt_nv_coarse && adj_coarse);

    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        // Count marked adjacencies
        int tmp_cnt = 0;
        for (int index_i = 0; index_i < data_agg->np; ++index_i)
        {
            if (mark[index][index_i] == 1)
            {
                ++tmp_cnt;
            }
        }
        cnt_nv_coarse[index] = tmp_cnt;

        // Allocate and fill adjacency list
        if (tmp_cnt > 0)
        {
            adj_coarse[index] = (int *)malloc(tmp_cnt * sizeof(int));
            assert(adj_coarse[index]);

            int pos = 0;
            for (int index_i = 0; index_i < data_agg->np; ++index_i)
            {
                if (mark[index][index_i] == 1)
                {
                    adj_coarse[index][pos] = index_i;
                    ++pos;
                }
            }
        }
        else
        {
            adj_coarse[index] = NULL;
        }
    }

    // Build CSR structure
    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        data_mesh_c->data_adj.xadj[index + 1] = data_mesh_c->data_adj.xadj[index] + cnt_nv_coarse[index];
    }

    int adj_nnv = data_mesh_c->data_adj.xadj[data_mesh_c->data_adj.local_nv];
    data_mesh_c->data_adj.adjncy = (int *)malloc(adj_nnv * sizeof(int));
    assert(data_mesh_c->data_adj.adjncy);

    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        if (cnt_nv_coarse[index] > 0)
        {
            int index_start = data_mesh_c->data_adj.xadj[index];
            memcpy(data_mesh_c->data_adj.adjncy + index_start,
                   adj_coarse[index],
                   cnt_nv_coarse[index] * sizeof(int));
        }
    }

    // Free memory
    for (int index = 0; index < data_mesh_c->data_adj.local_nv; ++index)
    {
        free(mark[index]);
        if (adj_coarse[index])
        {
            free(adj_coarse[index]);
        }
    }
    free(mark);
    free(cnt_nv_coarse);
    free(adj_coarse);

    return 0;
}

int SAMGCoarseVertexCoordinate(AggData **agg /*aggregation data*/,
                               MeshData **mesh_f /*fine-level mesh data*/,
                               MeshData **mesh_c /*coarse-level mesh data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    MeshData *data_mesh_c = *mesh_c;
    AggData *data_agg = *agg;

    data_mesh_c->data_vtx.nv = data_mesh_c->nv;
    data_mesh_c->data_vtx.local_nv = data_mesh_c->vtxdist[my_rank + 1] - data_mesh_c->vtxdist[my_rank];
    data_mesh_c->data_vtx.idx = (int *)malloc(data_mesh_c->data_vtx.local_nv * sizeof(int));
    data_mesh_c->data_vtx.type = (int *)malloc(data_mesh_c->data_vtx.local_nv * sizeof(int));
    assert(data_mesh_c->data_vtx.idx && data_mesh_c->data_vtx.type);

    for (int index = 0; index < data_mesh_c->data_vtx.local_nv; ++index)
    {
        data_mesh_c->data_vtx.idx[index] = index + data_mesh_c->vtxdist[my_rank];

        int fine_vtx_start = data_mesh_f->vtxdist[my_rank];
        int fine_vtx = data_agg->local_gids[index + data_mesh_c->vtxdist[my_rank]][0];
        int local_fine_vtx = fine_vtx - fine_vtx_start;
        data_mesh_c->data_vtx.type[index] = data_mesh_f->data_vtx.type[local_fine_vtx];
    }

    return 0;
}

#if 0
static int CheckValueInArray(int *a, int size, int val)
{
    for (int index = 0; index < size; ++index)
    {
        if (a[index] == val)
        {
            return 1;
        }
    }
    return 0;
}
#endif

int SAMGPartitionGhostDataMapping(AggData **agg /*aggregation data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    AggData *data_agg = *agg;

    data_agg->data_ghost_agg = (GhostAggData *)malloc(data_agg->np * sizeof(GhostAggData));
    assert(data_agg->data_ghost_agg);

    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        if (data_agg->owner[index_p] != my_rank)
        {
            continue;
        }

        GhostAggData *g = &data_agg->data_ghost_agg[index_p];
        int nlocal = data_agg->nlocal[index_p];
        int *local_gids = data_agg->local_gids[index_p];
        int n_fine = data_agg->n_fine[index_p];
        int *fine_gids = data_agg->fine_gids[index_p];

        g->n_owned = nlocal;
        g->n_total = n_fine;
        g->nghost = g->n_total - g->n_owned;

        g->flid2fgid = (int *)malloc(g->n_total * sizeof(int));
        // g->fgid2flid = (int *)malloc(g->n_total * sizeof(int));
        assert(g->flid2fgid);
        // assert(g->flid2fgid && g->fgid2flid);

        int lid_local = 0, lid_ghost = g->n_owned;
        for (int index_i = 0; index_i < n_fine; ++index_i)
        {
            int gid = fine_gids[index_i];
            g->flid2fgid[index_i] = gid;
#if 0
            if (index_i < nlocal)
            {
                g->fgid2flid[lid_local] = lid_local;
                ++lid_local;
            }
            else
            {
                g->fgid2flid[lid_ghost] = lid_ghost;
                ++lid_ghost;
            }
#endif
        }
    }

    return 0;
}

int SAMGPartitionRenumberID(AggData **agg /*aggregation data*/,
                            MeshData **mesh_f /*fine-level mesh data*/,
                            MeshData **mesh_c /*coarse-level mesh data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    MeshData *data_mesh_c = *mesh_c;
    AggData *data_agg = *agg;

    int *old2new = (int *)malloc(data_agg->np * sizeof(int)); // mapping old partition_id to new partition_id
    assert(old2new);
    memset(old2new, -1, data_agg->np * sizeof(int));

    for (int index = 0; index < nprocs; ++index)
    {
        int start = data_mesh_c->vtxdist[index];
        int tmp_cnt = 0;

        for (int index_p = 0; index_p < data_agg->np; ++index_p)
        {
            if (data_agg->owner[index_p] == index)
            {
                old2new[index_p] = start + tmp_cnt;
                ++tmp_cnt;
            }
        }
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, mapping from old to new\n", index);
            printf("old_partition_id\t ->\t new_partition_id\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                printf("%d\t ->\t %d\n", index_p, old2new[index_p]);
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check mapping from old to new

#if 0
    // updating partition_id of fine-level mesh
    for (int index = 0; index < data_mesh_f->data_vtx.local_nv; ++index)
    {
        int old_part = data_mesh_f->parts[index];
        data_mesh_f->parts[index] = old2new[old_part];
    }
#endif // partition_id updating

    // updating owner rank of data_agg
    int *new_owner = (int *)malloc(data_agg->np * sizeof(int));
    assert(new_owner);
    for (int index = 0; index < data_agg->np; ++index)
    {
        new_owner[old2new[index]] = data_agg->owner[index];
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, new/old owner rank\n", index);
            printf("new_partition_id\t new_owner_rank\t old_partition_id\t old_owner_rank\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == my_rank)
                {
                    printf("%d\t %d\t %d\t %d\n", old2new[index_p], new_owner[old2new[index_p]],
                           index_p, data_agg->owner[index_p]);
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check new_owner rank

    // updating nlocal and local_gids
    int *new_nlocal = (int *)malloc(data_agg->np * sizeof(int));
    assert(new_nlocal);
    for (int index = 0; index < data_agg->np; ++index)
    {
        new_nlocal[old2new[index]] = data_agg->nlocal[index];
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf("in rank %d, new/old nlocal\n", index);
            printf("new_partition_id\t new_nlocal\t old_partition_id\t old_nlocal\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == index)
                {
                    printf("%d\t %d\t %d\t %d\n", old2new[index_p], new_nlocal[old2new[index_p]],
                           index_p, data_agg->nlocal[index_p]);
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check new_nlocal

    int **new_local_gids = (int **)malloc(data_agg->np * sizeof(int *));
    assert(new_local_gids);
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (new_nlocal[index] > 0)
        {
            new_local_gids[index] = (int *)malloc(new_nlocal[index] * sizeof(int));
            assert(new_local_gids[index]);
        }
    }
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->nlocal[index] > 0)
        {
            memcpy(new_local_gids[old2new[index]],
                   data_agg->local_gids[index],
                   data_agg->nlocal[index] * sizeof(int));
        }
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, new/old local_gids\n", index);
            printf("new_partition_id\t vertex_id\n");
            printf("old_partition_id\t vertex_id\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == index)
                {
                    printf("%d\t ", old2new[index_p]);
                    for (int index_v = 0; index_v < new_nlocal[old2new[index_p]]; ++index_v)
                    {
                        printf("%d\t ", new_local_gids[old2new[index_p]][index_v]);
                    }
                    printf("\n");
                    printf("%d\t ", index_p);
                    for (int index_v = 0; index_v < data_agg->nlocal[index_p]; ++index_v)
                    {
                        printf("%d\t ", data_agg->local_gids[index_p][index_v]);
                    }
                    puts("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check new_local_gids

    // updating n_fine and fine_gids
    int *new_n_fine = (int *)malloc(data_agg->np * sizeof(int));
    assert(new_n_fine);
    for (int index = 0; index < data_agg->np; ++index)
    {
        new_n_fine[old2new[index]] = data_agg->n_fine[index];
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, new/old n_fine\n", index);
            printf("new_partition_id\t new_n_fine\t old_partition_id\t old_n_fine\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == index)
                {
                    printf("%d\t %d\t %d\t %d\n", old2new[index_p], new_n_fine[old2new[index_p]],
                           index_p, data_agg->n_fine[index_p]);
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check new_n_fine

    int **new_fine_gids = (int **)malloc(data_agg->np * sizeof(int *));
    assert(new_fine_gids);
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (new_owner[index] == my_rank)
        {
            new_fine_gids[index] = (int *)malloc(new_n_fine[index] * sizeof(int));
            assert(new_fine_gids);
        }
    }
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            memcpy(new_fine_gids[old2new[index]],
                   data_agg->fine_gids[index],
                   data_agg->n_fine[index] * sizeof(int));
        }
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, new/old fine_gids\n", index);
            printf("new_partition_id\t vertex_id\n");
            printf("old_partition_id\t vertex_id\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == my_rank)
                {
                    printf("%d\t ", old2new[index_p]);
                    for (int index_v = 0; index_v < new_n_fine[old2new[index_p]]; ++index_v)
                    {
                        printf("%d\t ", new_fine_gids[old2new[index_p]][index_v]);
                    }
                    printf("\n");
                    printf("%d\t ", index_p);
                    for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
                    {
                        printf("%d\t ", data_agg->fine_gids[index_p][index_v]);
                    }
                    puts("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check new_fine_gids

    // free memory
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            free(data_agg->fine_gids[index]);
        }
    }
    free(data_agg->fine_gids);
    free(data_agg->n_fine);
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->nlocal[index] > 0)
        {
            free(data_agg->local_gids[index]);
        }
    }
    free(data_agg->local_gids);
    free(data_agg->nlocal);
    free(data_agg->owner);
    free(old2new);

    // updating data agg value
    data_agg->owner = new_owner;
    data_agg->nlocal = new_nlocal;
    data_agg->local_gids = new_local_gids;
    data_agg->n_fine = new_n_fine;
    data_agg->fine_gids = new_fine_gids;

    return 0;
}

int SAMGPartitionOwnerVertexData(AggData **agg /*aggregation data*/,
                                 MeshData **mesh_f /*fine-level mesh data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    AggData *data_agg = *agg;

    data_agg->n_fine = (int *)malloc(data_agg->np * sizeof(int));
    data_agg->fine_gids = (int **)malloc(data_agg->np * sizeof(int *));
    assert(data_agg->n_fine && data_agg->fine_gids);

    (void)MPI_Allreduce(data_agg->nlocal,
                        data_agg->n_fine,
                        data_agg->np,
                        MPI_INT,
                        MPI_SUM,
                        comm);

    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            data_agg->fine_gids[index] = (int *)malloc(data_agg->n_fine[index] * sizeof(int));
            assert(data_agg->fine_gids[index]);

            memcpy(data_agg->fine_gids[index],
                   data_agg->local_gids[index],
                   data_agg->nlocal[index] * sizeof(int));
        }
        else
        {
            data_agg->fine_gids[index] = NULL;
        }
    }

    int *nlocal_all = (int *)malloc(nprocs * data_agg->np * sizeof(int));
    assert(nlocal_all);
    (void)MPI_Allgather(data_agg->nlocal,
                        data_agg->np,
                        MPI_INT,
                        nlocal_all,
                        data_agg->np,
                        MPI_INT,
                        comm);

    // all2all
    int *send_counts = (int *)calloc(nprocs, sizeof(int));
    int *recv_counts = (int *)calloc(nprocs, sizeof(int));
    int *sdispls = (int *)malloc(nprocs * sizeof(int));
    int *rdispls = (int *)malloc(nprocs * sizeof(int));
    assert(send_counts && recv_counts && sdispls && rdispls);

    // sending counts
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        int target = data_agg->owner[index_p];
        if (target != my_rank &&
            data_agg->nlocal[index_p] > 0)
        {
            send_counts[target] += data_agg->nlocal[index_p];
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
                recv_counts[index_r] += nlocal_all[index_r * data_agg->np + index_p];
            }
        }
    }

    // sending displacements
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int index_r = 1; index_r < nprocs; ++index_r)
    {
        sdispls[index_r] = sdispls[index_r - 1] + send_counts[index_r - 1];
        rdispls[index_r] = rdispls[index_r - 1] + recv_counts[index_r - 1];
    }

    // data packaging
    int total_send = sdispls[nprocs - 1] + send_counts[nprocs - 1];
    int *send_buf = (int *)malloc((total_send > 0 ? total_send : 1) * sizeof(int));
    int *tmp_offset = (int *)calloc(nprocs, sizeof(int));
    assert(send_buf && tmp_offset);

    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        int target = data_agg->owner[index_p];
        if (target != my_rank &&
            data_agg->nlocal[index_p] > 0)
        {
            memcpy(send_buf + sdispls[target] + tmp_offset[target],
                   data_agg->local_gids[index_p],
                   data_agg->nlocal[index_p] * sizeof(int));
            tmp_offset[target] += data_agg->nlocal[index_p];
        }
    }

    // mpi_all2allv
    int total_recv = rdispls[nprocs - 1] + recv_counts[nprocs - 1];
    int *recv_buf = (int *)malloc((total_recv > 0 ? total_recv : 1) * sizeof(int));
    assert(recv_buf);

    (void)MPI_Alltoallv(send_buf, send_counts, sdispls, MPI_INT,
                        recv_buf, recv_counts, rdispls, MPI_INT, comm);

    // data unpackaging, vertices order [local, remote]
    int *part_fill_offset = (int *)malloc(data_agg->np * sizeof(int));
    assert(part_fill_offset);
    for (int index_p = 0; index_p < data_agg->np; ++index_p)
    {
        part_fill_offset[index_p] = (data_agg->owner[index_p] == my_rank) ? data_agg->nlocal[index_p] : 0;
    }

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
            if (cnt > 0 &&
                data_agg->owner[index_p] == my_rank)
            {
                memcpy(data_agg->fine_gids[index_p] + part_fill_offset[index_p],
                       recv_buf + rdispls[index_r] + tmp_offset[index_r],
                       cnt * sizeof(int));
                part_fill_offset[index_p] += cnt;
                tmp_offset[index_r] += cnt;
            }
        }
    }

    // free memory
    free(part_fill_offset);
    free(recv_buf);
    free(send_buf);
    free(tmp_offset);
    free(send_counts);
    free(recv_counts);
    free(sdispls);
    free(rdispls);
    free(nlocal_all);

    return 0;
}

int SAMGPartitionOwnerRankConstructor(AggData **agg /*aggregation data*/,
                                      MeshData **mesh_f /*fine-level mesh data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    AggData *data_agg = *agg;

    int *flag_parts = (int *)calloc(data_mesh_f->np, sizeof(int));
    int *owner_rank_candidate_parts = (int *)malloc(data_mesh_f->np * sizeof(int));
    data_agg->owner = (int *)malloc(data_mesh_f->np * sizeof(int));
    assert(flag_parts && owner_rank_candidate_parts && data_agg->owner);

    for (int index = 0; index < data_mesh_f->data_vtx.local_nv; ++index)
    {
        flag_parts[data_mesh_f->parts[index]] = 1;
    }
    for (int index = 0; index < data_mesh_f->np; ++index)
    {
        owner_rank_candidate_parts[index] = flag_parts[index] ? my_rank : INT_MAX;
    }
    (void)MPI_Allreduce(owner_rank_candidate_parts,
                        data_agg->owner,
                        data_mesh_f->np,
                        MPI_INT,
                        MPI_MIN,
                        comm);

    // free memory
    free(owner_rank_candidate_parts);
    free(flag_parts);

    return 0;
}

int SAMGCoarseMeshConstructor(MeshData **mesh_f /*fine-level mesh data*/,
                              MeshData **mesh_c /*coarse-level mesh data*/,
                              AggData **agg /*aggregation data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    MeshData *data_mesh_f = *mesh_f;
    MeshData *data_mesh_c = *mesh_c;
    AggData *data_agg = *agg;

    data_agg->np = data_mesh_f->np;
    data_agg->nv = data_mesh_f->nv;

    // partition owner rank
    PetscCall(SAMGPartitionOwnerRankConstructor(agg, mesh_f));

    // coarse vtxdist
    int local_nv_c = 0;
    for (int index = 0; index < data_mesh_f->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            ++local_nv_c;
        }
    }

    data_mesh_c->data_vtx.local_nv = local_nv_c;
    data_mesh_c->data_adj.local_nv = local_nv_c;

    int *nv_rank_c = (int *)malloc(nprocs * sizeof(int));
    assert(nv_rank_c);
    (void)MPI_Allgather(&local_nv_c,
                        1,
                        MPI_INT,
                        nv_rank_c,
                        1,
                        MPI_INT,
                        comm);

    data_mesh_c->vtxdist = (int *)malloc((nprocs + 1) * sizeof(int));
    assert(data_mesh_c->vtxdist);

    data_mesh_c->vtxdist[0] = 0;
    for (int index = 0; index < nprocs; ++index)
    {
        data_mesh_c->vtxdist[index + 1] = data_mesh_c->vtxdist[index] + nv_rank_c[index];
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, owner rank:", index);
            printf("partition_id\t owner_rank\n");
            for (int index_i = 0; index_i < data_mesh_f->np; ++index_i)
            {
                if (data_agg->owner[index_i] == my_rank)
                {
                    printf("%d\t %d\n", index_i, my_rank);
                }
            }
            printf("vtxdist[]:\n");
            for (int index_i = 0; index_i < nprocs + 1; ++index_i)
            {
                printf("%d\t", data_mesh_c->vtxdist[index_i]);
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check owner rank and vtxdist

    // partition local data
    data_agg->nlocal = (int *)calloc(data_agg->np, sizeof(int));
    data_agg->local_gids = (int **)malloc(data_agg->np * sizeof(int *));
    assert(data_agg->nlocal && data_agg->local_gids);

    for (int index = 0; index < data_mesh_f->data_vtx.local_nv; ++index)
    {
        ++(data_agg->nlocal[data_mesh_f->parts[index]]);
    }
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->nlocal[index] > 0)
        {
            data_agg->local_gids[index] = (int *)malloc(data_agg->nlocal[index] * sizeof(int));
            assert(data_agg->local_gids[index]);
        }
    }

    int *cnt_local_gids = (int *)calloc(data_agg->np, sizeof(int));
    assert(cnt_local_gids);
    for (int index = 0; index < data_mesh_f->data_vtx.local_nv; ++index)
    {
        int p = data_mesh_f->parts[index];
        data_agg->local_gids[p][cnt_local_gids[p]] = data_mesh_f->data_vtx.idx[index];
        ++(cnt_local_gids[p]);
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, partition local data\n", index);
            printf("partition_id\t vertex_id\n");
            for (int index_i = 0; index_i < data_agg->np; ++index_i)
            {
                if (data_agg->nlocal[index_i] > 0)
                {
                    printf("%d\t ", index_i);
                    for (int index_j = 0; index_j < data_agg->nlocal[index_i]; ++index_j)
                    {
                        printf("%d\t ", data_agg->local_gids[index_i][index_j]);
                    }
                    printf("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check partition local data

    // partition owner vertex data
    PetscCall(SAMGPartitionOwnerVertexData(agg, mesh_f));
    PetscCall(SAMGPartitionRenumberID(agg, mesh_f, mesh_c));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, partition owner vertex data\n", index);
            printf("partition_id\t vertex_id\n");
            for (int index_i = 0; index_i < data_agg->np; ++index_i)
            {
                if (data_agg->owner[index_i] == my_rank)
                {
                    printf("%d\t ", index_i);
                    for (int index_j = 0; index_j < data_agg->n_fine[index_i]; ++index_j)
                    {
                        printf("%d\t ", data_agg->fine_gids[index_i][index_j]);
                    }
                    printf("\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check partition owner vertex data

    // mapping from fine-level vertex to partition id
    PetscCall(SAMGFineVertex2PartitionMap(agg));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (my_rank == index)
        {
            printf(">>>> in rank %d, fine-level vertex to partition id\n", index);
            printf("fine-level vertex_id\t partition_id\n");
            int index_fine_vtx_start = data_mesh_f->vtxdist[my_rank];
            int index_fine_vtx_end = data_mesh_f->vtxdist[my_rank + 1];
            for (int index_v = index_fine_vtx_start; index_v < index_fine_vtx_end; ++index_v)
            {
                printf("%d\t %d\n", index_v, data_agg->fgid2part[index_v]);
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check fine-level vertex mapping to partition id

    // construct ghost
    PetscCall(SAMGPartitionGhostDataMapping(agg));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, partition ghost mapping data\n", index);
            printf("partition_id\t vertex_id\n");
            printf("\t mapping vertex_id\n");
            for (int index_p = 0; index_p < data_agg->np; ++index_p)
            {
                if (data_agg->owner[index_p] == my_rank)
                {
                    printf("%d\t ", index_p);
                    for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
                    {
                        printf("%d\t ", data_agg->fine_gids[index_p][index_v]);
                    }
                    printf("\n\t ");
                    for (int index_v = 0; index_v < data_agg->n_fine[index_p]; ++index_v)
                    {
                        printf("%d\t ", data_agg->data_ghost_agg[index_p].fgid2flid[index_v]);
                    }
                    printf("\n\n");
                }
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check partition ghost mapping

    // construct coarse-level vertex coordinate
    PetscCall(SAMGCoarseVertexCoordinate(agg, mesh_f, mesh_c));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (my_rank == index)
        {
            printf(">>>> in rank %d, coarse-level vertex coordinate data\n", index);
            printf("vertex_id\t type\n");
            for (int index_v = 0; index_v < data_mesh_c->data_vtx.local_nv; ++index_v)
            {
                printf("%d\t %d\n", data_mesh_c->data_vtx.idx[index_v],
                       data_mesh_c->data_vtx.type[index_v]);
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check coarse-level vertex coordinate data

    // construct coarse-level adjacent list
    PetscCall(SAMGCoarseAdjacentListConstructor(agg, mesh_f, mesh_c));
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, coarse-level adjacent list data\n", index);
            printf("local_coarse_vtx\t global_coarse_vtx\t adjacent_vtxs\n");
            for (int index_v = 0; index_v < data_mesh_c->data_adj.local_nv; ++index_v)
            {
                printf("%d\t %d\t ", index_v, data_mesh_c->data_adj.idx[index_v]);
                int index_v_start = data_mesh_c->data_adj.xadj[index_v];
                int index_v_end = data_mesh_c->data_adj.xadj[index_v + 1];
                for (int index_v_adj = index_v_start; index_v_adj < index_v_end; ++index_v_adj)
                {
                    printf("%d\t ", data_mesh_c->data_adj.adjncy[index_v_adj]);
                }
                printf("\n");
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check coarse-level adjacent list data

    // free memory
    free(cnt_local_gids);
    free(nv_rank_c);

    return 0;
}

int SAMGLevel0Mesh(const CfgJson *data_cfg /*config data*/,
                   MeshData *data_f_mesh /*fine-level mesh data*/,
                   MeshData *data_c_mesh /*coarse-level mesh data*/,
                   AggData *data_agg /*aggregation data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    FileProcessMeshVtx(data_cfg->cfg_file.file_vtx,
                       &data_f_mesh->data_vtx);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, level 0 mesh vertex data, total nv = %d\n", index, data_f_mesh->data_vtx.nv);
            printf("local_vtx_id\t global_vtx_id\t type\t "
                   "x\t y\t z\t "
                   "nx\t ny\t nz\n");
            for (int index_v = 0; index_v < data_f_mesh->data_vtx.local_nv; ++index_v)
            {
                printf("%d\t %d\t %d\t "
                       "%021.16le\t %021.16le\t %021.16le\t "
                       "%021.16le\t %021.16le\t %021.16le\n",
                       index_v, data_f_mesh->data_vtx.idx[index_v], data_f_mesh->data_vtx.type[index_v],
                       data_f_mesh->data_vtx.data_coor[index_v].x,
                       data_f_mesh->data_vtx.data_coor[index_v].y,
                       data_f_mesh->data_vtx.data_coor[index_v].z,
                       data_f_mesh->data_vtx.data_coor[index_v].nx,
                       data_f_mesh->data_vtx.data_coor[index_v].ny,
                       data_f_mesh->data_vtx.data_coor[index_v].nz);
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check level 0 mesh vertex data

    FileProcessMeshAdj(data_cfg->cfg_file.file_adj,
                       &data_f_mesh->data_adj);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        (void)MPI_Barrier(comm);
        if (index == my_rank)
        {
            printf(">>>> in rank %d, level 0 mesh adjacent list data, total nv = %d, local nv = %d\n", index, data_f_mesh->data_adj.nv, data_f_mesh->data_adj.local_nv);
            printf("local_vtx_id\t global_vtx_id\t adj_vtxs\n");
            for (int index_v = 0; index_v < data_f_mesh->data_adj.local_nv; ++index_v)
            {
                printf("%d\t %d\t ", index_v, data_f_mesh->data_adj.idx[index_v]);
                int index_v_start = data_f_mesh->data_adj.xadj[index_v];
                int index_v_end = data_f_mesh->data_adj.xadj[index_v + 1];
                for (int index_v_adj = index_v_start; index_v_adj < index_v_end; ++index_v_adj)
                {
                    printf("%d\t ", data_f_mesh->data_adj.adjncy[index_v_adj]);
                }
                printf("\n");
            }
            puts("\n");
        }
        fflush(stdout);
    }
#endif // check level 0 mesh adjacent list data

    data_f_mesh->nv = data_f_mesh->data_vtx.nv;

    int cfg_mg_est_size_agg = data_cfg->cfg_mg.est_size_agg;
    data_f_mesh->np = data_f_mesh->nv / cfg_mg_est_size_agg;

    data_c_mesh->nv = data_f_mesh->np;
    data_c_mesh->np = data_c_mesh->nv / cfg_mg_est_size_agg;
    data_c_mesh->data_vtx.nv = data_c_mesh->nv;
    data_c_mesh->data_adj.nv = data_c_mesh->nv;

    data_f_mesh->vtxdist = (int *)malloc((nprocs + 1) * sizeof(int));
    assert(data_f_mesh->vtxdist);

    data_f_mesh->vtxdist[0] = 0;

    int *nv_rank = (int *)malloc(nprocs * sizeof(int));
    assert(nv_rank);

    // data_f_mesh->vtxdist
    for (int index = 0; index < nprocs; ++index)
    {
        nv_rank[index] = data_f_mesh->nv / nprocs;
    }
    for (int index = 0; index < data_f_mesh->nv % nprocs; ++index)
    {
        nv_rank[index] += 1;
    }
    for (int index = 0; index < nprocs; ++index)
    {
        data_f_mesh->vtxdist[index + 1] = data_f_mesh->vtxdist[index] + nv_rank[index];
    }
#if 0
    for (int rank = 0; rank < nprocs; ++rank)
    {
        if (rank == my_rank)
        {
            printf(">>>> in rank %d:\n", rank);
            for (int index = 0; index < nprocs + 1; ++index)
            {
                printf("in fine-level mesh, vtxdist[%d] = %d\n", index, data_f_mesh->vtxdist[index]);
            }
        }
    }
#endif // check vtxdist

    data_f_mesh->parts = (int *)malloc(data_f_mesh->data_vtx.local_nv * sizeof(int));
    assert(data_f_mesh->parts);

    // calling ParMetis
    int wgtflag = 0, numflag = 0, ncon = 1;
    int nparts = data_f_mesh->np;
    real_t *tpwgts = NULL, ubvec = 1.05; // double
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_SEED] = 42; // for any integer
    idx_t edgecut = 0;

    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(tpwgts);
    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / ncon / nparts;
    }

    int metis_status = ParMETIS_V3_PartKway(data_f_mesh->vtxdist,
                                            data_f_mesh->data_adj.xadj,
                                            data_f_mesh->data_adj.adjncy,
                                            NULL, NULL,
                                            &wgtflag, &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, data_f_mesh->parts, &comm);
#if 0
    for (int rank = 0; rank < nprocs; ++rank)
    {
        (void)MPI_Barrier(comm);
        if (rank == my_rank)
        {
            printf(">>>> in rank %d, partition result:\n", rank);
            printf("local vertex\t global vertex\t partition\n");
            for(int index = 0; index < data_f_mesh->data_vtx.local_nv; ++index)
            {
                printf("%d\t %d\t %d\n", index, data_f_mesh->vtxdist[rank] + index, data_f_mesh->parts[index]);
            }
            printf("\n");
        }
        fflush(stdout);
    }
#endif // check partition result

    // coarse-level mesh data construction
    PetscCall(SAMGCoarseMeshConstructor(&data_f_mesh, &data_c_mesh, &data_agg));

    // free memory
    free(tpwgts);
    free(nv_rank);

    return 0;
}

int SAMGLevelMesh(int cfg_mg_num_level /*config number of levels*/,
                  SAMGCtx **samg_ctx /*samg context data*/)
{
    SAMGCtx *data_samg_ctx = *samg_ctx;

    int cfg_mg_num_coarse_vtx = data_samg_ctx->data_cfg.cfg_mg.num_coarse_vtx;

    data_samg_ctx->levels = (MGLevel *)malloc(cfg_mg_num_level * sizeof(MGLevel));
    assert(data_samg_ctx->levels);

    int cnt_level = 0;
    PetscCall(SAMGLevel0Mesh(&data_samg_ctx->data_cfg,
                             &data_samg_ctx->levels[cnt_level].data_f_mesh,
                             &data_samg_ctx->levels[cnt_level].data_c_mesh,
                             &data_samg_ctx->levels[cnt_level].data_agg));

    while (cnt_level < cfg_mg_num_level &&
           data_samg_ctx->levels[cnt_level].data_c_mesh.nv > cfg_mg_num_coarse_vtx)
    {
        ++cnt_level;
    }

    data_samg_ctx->num_level = cnt_level + 1;

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

    // redo nlocal_all in fine_gids constructor
    int *nlocal_all = (int *)malloc(nprocs * data_agg->np * sizeof(int));
    assert(nlocal_all);
    (void)MPI_Allgather(data_agg->nlocal, data_agg->np, MPI_INT,
                        nlocal_all, data_agg->np, MPI_INT, comm);

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
            continue;

        for (int index_p = 0; index_p < data_agg->np; ++index_p)
        {
            int cnt = nlocal_all[index_r * data_agg->np + index_p];
            if (cnt > 0 && data_agg->owner[index_p] == my_rank)
            {
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
    free(nlocal_all);
    free(part_fill_offset);

    return 0;
}

static void QRFactorizationOpenBLAS(int m /*nrow*/, int n /*ncol*/,
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
#if 1
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
#endif // check size

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

    // set value for tentative prolongation operator
    for(int index_p = 0; index_p < data_level->data_agg.np; ++index_p)
    {
        if(data_level->data_agg.owner[index_p] == my_rank)
        {}
    }

    return 0;
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
    }

    return 0;
}

int SAMGSetupPhase(SAMGCtx *samg_ctx /*samg context data*/,
                   int sa_flag /*flag of sa*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    int cfg_mg_num_level = samg_ctx->data_cfg.cfg_mg.num_level;
    // int cfg_mg_pre_smooth = samg_ctx->data_cfg.cfg_mg.pre_smooth;
    // int cfg_mg_post_smooth = samg_ctx->data_cfg.cfg_mg.post_smooth;
    // int cfg_mg_num_coarse_vtx = samg_ctx->data_cfg.cfg_mg.num_coarse_vtx;
    // int cfg_mg_est_size_agg = samg_ctx->data_cfg.cfg_mg.est_size_agg;
    // int cfg_mg_ps_num_steps = samg_ctx->data_cfg.cfg_mg.ps_num_steps;
    // int cfg_mg_ps_type = samg_ctx->data_cfg.cfg_mg.ps_type;
    // double cfg_mg_ps_scale = samg_ctx->data_cfg.cfg_mg.ps_scale;

    // samg_ctx->levels = (MGLevel *)malloc(cfg_mg_num_level * sizeof(MGLevel));
    // assert(samg_ctx->levels);

    int cnt_level = 0;

    PetscCall(SAMGLevelMesh(cfg_mg_num_level, &samg_ctx)); // multilevel hierarchy
    PetscCall(MatDuplicate(samg_ctx->mysolver.solver_a, MAT_COPY_VALUES, &samg_ctx->levels[cnt_level].op_f));

    PetscCall(SAMGInitialNearNullSpace(&samg_ctx->levels[cnt_level].data_f_mesh,
                                       &samg_ctx->data_nullspace_level0)); // level 0 near null space
    PetscCall(SAMGLevelNearNullSpace(&samg_ctx));                          // multilevel near null space

    PetscCall(SAMGTentativeProlongationOperator(&samg_ctx)); // tentative prolongation operator constructor

    if (sa_flag == 1)
    {
        // SA
        PetscCall(PetscPrintf(comm, "==== Smoothed Aggregation-based Multigrid\n"));
    }
    else if (sa_flag == 0)
    {
        // UA
        PetscCall(PetscPrintf(comm, "==== Unsmoothed Aggregation-based Multigrid\n"));
    }

    // while (cnt_level < cfg_mg_num_level &&
    //        samg_ctx->levels[cnt_level].data_f_mesh.np > cfg_mg_num_coarse_vtx)
    // {
    //     ++cnt_level;
    // }

    return 0;
}
