#include "../include/main.h"

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
        g->fgid2flid = (int *)malloc(g->n_total * sizeof(int));
        assert(g->flid2fgid && g->fgid2flid);

        int lid_local = 0, lid_ghost = g->n_owned;
        for (int index_i = 0; index_i < n_fine; ++index_i)
        {
            int gid = fine_gids[index_i];
            if (CheckValueInArray(local_gids, nlocal, gid) == 1)
            {
                g->flid2fgid[lid_local] = gid;
                g->fgid2flid[lid_local] = lid_local;
                ++lid_local;
            }
            else
            {
                g->flid2fgid[lid_ghost] = gid;
                g->fgid2flid[lid_ghost] = lid_ghost;
                ++lid_ghost;
            }
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

    int *offset = (int *)calloc(data_agg->np, sizeof(int));
    assert(offset);

    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] == my_rank)
        {
            offset[index] = 0;

            if (data_agg->nlocal[index] > 0)
            {
                memcpy(data_agg->fine_gids[index],
                       data_agg->local_gids[index],
                       data_agg->nlocal[index] * sizeof(int));
                offset[index] += data_agg->nlocal[index];
            }
        }
    }

    // sending non-owner data
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->nlocal[index] > 0 &&
            data_agg->owner[index] != my_rank)
        {
            (void)MPI_Send(data_agg->local_gids[index],
                           data_agg->nlocal[index],
                           MPI_INT,
                           data_agg->owner[index],
                           1000 + index, // tag, partition_id + 1000
                           comm);
        }
    }

    // receiving non-owner data in owner rank
    MPI_Status status;
    for (int index = 0; index < data_agg->np; ++index)
    {
        if (data_agg->owner[index] != my_rank)
        {
            continue;
        }

        for (int index_r = 0; index_r < nprocs; ++index_r)
        {
            if (index_r == my_rank)
            {
                continue;
            }

            int tmp_cnt = nlocal_all[index_r * data_agg->np + index];
            if (tmp_cnt <= 0)
            {
                continue;
            }

            MPI_Recv(data_agg->fine_gids[index] + offset[index],
                     tmp_cnt,
                     MPI_INT,
                     index_r,
                     1000 + index, // tag, partition_id + 1000
                     comm,
                     &status);

            offset[index] += tmp_cnt;
        }
    }

    // free memory
    free(nlocal_all);
    free(offset);

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

    // construct ghost
    PetscCall(SAMGPartitionGhostDataMapping(agg));
#if 1
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
    FileProcessMeshAdj(data_cfg->cfg_file.file_adj,
                       &data_f_mesh->data_adj);

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

    data_samg_ctx->num_level = cnt_level;

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

int SAMGSetupPhase(SAMGCtx *samg_ctx /*samg context data*/)
{
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

    PetscCall(SAMGLevelMesh(cfg_mg_num_level, &samg_ctx));

    int cnt_level = 0;
    PetscCall(MatDuplicate(samg_ctx->mysolver.solver_a, MAT_COPY_VALUES, &samg_ctx->levels[cnt_level].op_f));

    // while (cnt_level < cfg_mg_num_level &&
    //        samg_ctx->levels[cnt_level].data_f_mesh.np > cfg_mg_num_coarse_vtx)
    // {
    //     ++cnt_level;
    // }

    return 0;
}
