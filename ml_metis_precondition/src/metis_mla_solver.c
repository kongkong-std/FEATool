#include "../include/main.h"

int TestMetisFunctionGmsh(DataGmsh data)
{
    int status = 0;

    idx_t ne = data.ne_in, nn = data.nn;
    idx_t nparts = data.nparts;
    idx_t options[METIS_NOPTIONS];
    idx_t objval;
    // idx_t *epart, *npart;

#if 0
    epart = (idx_t *)malloc(ne * sizeof(idx_t));
    npart = (idx_t *)malloc(nn * sizeof(idx_t));
    assert(epart && npart);
#endif

    METIS_SetDefaultOptions(options);

#if 0
    for (int index = 0; index < ne + 1; ++index)
    {
        printf("data.eptr_in[%d] = %ld\n", index, data.eptr_in[index]);
    }
    for (int index = 0; index < data.eptr_in[ne]; ++index)
    {
        printf("data.eind_in[%d] = %ld\n", index, data.eind_in[index]);
    }
#endif // mesh csr format information

    status = METIS_PartMeshNodal(&ne, &nn,
                                 data.eptr_in, data.eind_in,
                                 NULL, NULL,
                                 &nparts,
                                 NULL, options,
                                 &objval, data.epart_in, data.npart_in);

#if 1
    puts("\n==== metis function result ====");
    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    printf("partition objective value = %" PRIDX "\n", objval);
    printf("element partition:\n");
    for (int index = 0; index < ne; ++index)
    {
        printf("epart[%d] = %" PRIDX "\n", index, data.epart_in[index]);
    }
    printf("node partition:\n");
    for (int index = 0; index < nn; ++index)
    {
        // printf("npart[%d] = %ld\n", index, npart[index]);
        printf("npart[%d] = %" PRIDX "\n", index, data.npart_in[index]);
    }
#endif // metis function result

#if 0
    // free memory
    free(epart);
    free(npart);
#endif

    return status;
}

int DeepCopyKLevelMetisDataMesh(DataGmsh *dst, DataGmsh *src)
{
    dst->nn = src->nn;
    dst->ne_in = src->ne_in;
    dst->nparts = src->nparts;
    dst->nne_bd = src->nne_bd;
    dst->nne_in = src->nne_in;

    dst->coordinates = (double *)malloc(3 * dst->nn * sizeof(double));
    assert(dst->coordinates);
    memcpy(dst->coordinates, src->coordinates, 3 * dst->nn * sizeof(double));

    dst->eptr_in = (idx_t *)malloc((dst->ne_in + 1) * sizeof(idx_t));
    assert(dst->eptr_in);
    memcpy(dst->eptr_in, src->eptr_in, (dst->ne_in + 1) * sizeof(idx_t));

    dst->eind_in = (idx_t *)malloc(dst->eptr_in[dst->nn] * sizeof(idx_t));
    assert(dst->eind_in);
    memcpy(dst->eind_in, src->eind_in, dst->eptr_in[dst->nn] * sizeof(idx_t));

    dst->npart_in = (idx_t *)malloc(dst->nn * sizeof(idx_t));
    assert(dst->npart_in);
    memcpy(dst->npart_in, src->npart_in, dst->nn * sizeof(idx_t));

    return 0;
}

int DeepCopyMetisDataGmsh(DataGmsh *dst, DataGmsh *src)
{
    dst->nn = src->nn;
    dst->ne = src->ne;
    dst->ne_bd = src->ne_bd;
    dst->ne_in = src->ne_in;
    dst->nparts = src->nparts;
    dst->nne_bd = src->nne_bd;
    dst->nne_in = src->nne_in;
    dst->coordinates = (double *)malloc(3 * dst->nn * sizeof(double));
    dst->npart_in = (idx_t *)malloc(dst->nn * sizeof(idx_t));
    dst->epart_in = (idx_t *)malloc(dst->ne_in * sizeof(idx_t));
    dst->eptr_in = (idx_t *)malloc(dst->ne * sizeof(idx_t));
    dst->eind_in = (idx_t *)malloc(dst->ne * dst->nne_in * sizeof(idx_t));
    assert(dst->coordinates && dst->npart_in && dst->eptr_in && dst->eind_in);

    memcpy(dst->coordinates, src->coordinates, 3 * dst->nn * sizeof(double));
    memcpy(dst->eptr_in, src->eptr_in, dst->ne * sizeof(idx_t));
    memcpy(dst->eind_in, src->eind_in, dst->ne * dst->nne_in * sizeof(idx_t));

    return 0;
}

int MetisMLASolverSolvePhase(MLAContext *mla_ctx)
{
    return 0;
}

void MetisKLevelCoarseLevelGenerator(DataGmsh *coarse_data,
                                     DataGmsh *fine_data)
{
    coarse_data->nn = fine_data->nparts;
    coarse_data->ne_in = coarse_data->nn;
    coarse_data->nparts = coarse_data->nn / 4;
    coarse_data->nne_bd = NumNodeEleTypeMap(1); // segment line element
    coarse_data->nne_in = NumNodeEleTypeMap(1); // segment line element

    // coarse level node coordinates
    coarse_data->coordinates = (double *)malloc(coarse_data->nn * 3 * sizeof(double));
    assert(coarse_data->coordinates);
    memset(coarse_data->coordinates, 0, 3 * coarse_data->nn * sizeof(double));

    int *cnt_node_partition = NULL; // nodes in current partition
    cnt_node_partition = (int *)malloc(coarse_data->nn * sizeof(int));
    assert(cnt_node_partition);
    memset(cnt_node_partition, 0, coarse_data->nn * sizeof(int));

    for (int index = 0; index < fine_data->nn; ++index)
    {
        idx_t id_part = fine_data->npart_in[index]; // 0-base, node partition
        coarse_data->coordinates[id_part * 3] += fine_data->coordinates[index * 3];
        coarse_data->coordinates[id_part * 3 + 1] += fine_data->coordinates[index * 3 + 1];
        coarse_data->coordinates[id_part * 3 + 2] += fine_data->coordinates[index * 3 + 2];
        ++(cnt_node_partition[id_part]);
    }

    for (int index = 0; index < coarse_data->nn; ++index)
    {
        coarse_data->coordinates[3 * index] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 1] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 2] /= cnt_node_partition[index];
    }

    // coarse level adjacency list
    coarse_data->eptr_in = (idx_t *)malloc((coarse_data->ne_in + 1) * sizeof(idx_t));
    assert(coarse_data->eptr_in);
    memset(coarse_data->eptr_in, 0, (coarse_data->ne_in + 1) * sizeof(idx_t));

    // adjacency matrix
    bool **mat_adj = NULL;
    mat_adj = (bool **)malloc(coarse_data->nn * sizeof(bool *));
    assert(mat_adj);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_data->nn, sizeof(bool));
        assert(mat_adj[index]);
        // mat_adj[index][index] = true;
    }

    // adjacency list in fine level traversal
    for (int index = 0; index < fine_data->ne_in; ++index)
    {
        idx_t coarse_node_i = fine_data->npart_in[index];
        idx_t index_start = fine_data->eptr_in[index];
        idx_t index_end = fine_data->eptr_in[index + 1];

        for (idx_t index_j = index_start; index_j < index_end; ++index_j)
        {
            idx_t fine_node_j = fine_data->eind_in[index_j];
            idx_t coarse_node_j = fine_data->npart_in[fine_node_j];

            mat_adj[coarse_node_i][coarse_node_j] = true;
            mat_adj[coarse_node_j][coarse_node_i] = true;
        }
    }

    // diagonal = 0
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index][index] = 0;
    }

    // assign value to coarse_data->eptr_in
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        int cnt_tmp = 0;
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                ++cnt_tmp;
            }
        }
        coarse_data->eptr_in[index_i + 1] = coarse_data->eptr_in[index_i] + cnt_tmp;
    }

    // assign value to coarse_data->eind_in
    coarse_data->eind_in = (idx_t *)malloc(coarse_data->eptr_in[coarse_data->nn] * sizeof(idx_t));
    assert(coarse_data);
    int pos_tmp = 0;
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                coarse_data->eind_in[pos_tmp] = index_j;
                ++pos_tmp;
            }
        }
    }

    coarse_data->npart_in = (idx_t *)malloc(coarse_data->nn * sizeof(idx_t));
    assert(coarse_data->npart_in);

    // free memory
    free(cnt_node_partition);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
}

void GmshCoarseLevelGenerator(DataGmsh *coarse_data /*gmsh coarse level data pointer*/,
                              DataGmsh *fine_data /*gmsh fine level data pointer*/)
{
    coarse_data->nn = fine_data->nparts;
    coarse_data->ne_in = coarse_data->nn; // adjacency list size, number of elements equal to number of nodes
    coarse_data->nparts = coarse_data->nn / 4;
    coarse_data->nne_bd = NumNodeEleTypeMap(1); // segment line element
    coarse_data->nne_in = NumNodeEleTypeMap(1); // segment line element

    // coarse level node coordinates
    coarse_data->coordinates = (double *)malloc(coarse_data->nn * 3 * sizeof(double));
    assert(coarse_data->coordinates);
    memset(coarse_data->coordinates, 0, 3 * coarse_data->nn * sizeof(double));

    int *cnt_node_partition = NULL; // nodes in current partition
    cnt_node_partition = (int *)malloc(coarse_data->nn * sizeof(int));
    assert(cnt_node_partition);
    memset(cnt_node_partition, 0, coarse_data->nn * sizeof(int));

    for (int index = 0; index < fine_data->nn; ++index)
    {
        idx_t id_part = fine_data->npart_in[index]; // 0-base
        coarse_data->coordinates[id_part * 3] += fine_data->coordinates[index * 3];
        coarse_data->coordinates[id_part * 3 + 1] += fine_data->coordinates[index * 3 + 1];
        coarse_data->coordinates[id_part * 3 + 2] += fine_data->coordinates[index * 3 + 2];
        ++(cnt_node_partition[id_part]);
    }

    for (int index = 0; index < coarse_data->nn; ++index)
    {
        coarse_data->coordinates[3 * index] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 1] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 2] /= cnt_node_partition[index];
    }

    // coarse level adjacency list
    coarse_data->eptr_in = (idx_t *)malloc((coarse_data->ne_in + 1) * sizeof(idx_t));
    assert(coarse_data->eptr_in);
    memset(coarse_data->eptr_in, 0, (coarse_data->ne_in + 1) * sizeof(idx_t));

    // adjacency matrix
    bool **mat_adj = NULL;
    mat_adj = (bool **)malloc(coarse_data->nn * sizeof(bool *));
    assert(mat_adj);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_data->nn, sizeof(bool));
        assert(mat_adj[index]);
        // mat_adj[index][index] = true;
    }

    // elements in fine level traversal
    for (int index = 0; index < fine_data->ne_in; ++index)
    {
        idx_t index_start = fine_data->eptr_in[index];
        idx_t index_end = fine_data->eptr_in[index + 1];

        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            idx_t fine_node_i = fine_data->eind_in[index_i];
            idx_t coarse_node_i = fine_data->npart_in[fine_node_i];

            for (idx_t index_j = index_i + 1; index_j < index_end; ++index_j)
            {
                idx_t fine_node_j = fine_data->eind_in[index_j];
                idx_t coarse_node_j = fine_data->npart_in[fine_node_j];

                mat_adj[coarse_node_i][coarse_node_j] = true;
                mat_adj[coarse_node_j][coarse_node_i] = true;
            }
        }
    }

    // diagonal = 0
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index][index] = 0;
    }

    // assign value to coarse_data->eptr_in
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        int cnt_tmp = 0;
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                ++cnt_tmp;
            }
        }
        coarse_data->eptr_in[index_i + 1] = coarse_data->eptr_in[index_i] + cnt_tmp;
    }

    // assign value to coarse_data->eind_in
    coarse_data->eind_in = (idx_t *)malloc(coarse_data->eptr_in[coarse_data->nn] * sizeof(idx_t));
    assert(coarse_data);
    int pos_tmp = 0;
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                coarse_data->eind_in[pos_tmp] = index_j;
                ++pos_tmp;
            }
        }
    }

    coarse_data->npart_in = (idx_t *)malloc(coarse_data->nn * sizeof(idx_t));
    assert(coarse_data->npart_in);

    // free memory
    free(cnt_node_partition);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
}

int MetisMLASolverSetupPhase(MLAContext *mla_ctx)
{
    if (mla_ctx->setup == 1)
    {
        // setup phase has been done
        return 0;
    }

    mla_ctx->setup = 1;

    idx_t options[METIS_NOPTIONS];
    idx_t objval;
    idx_t nn = 0, ne = 0;
    int status_metis_part;

    int config_num_level = mla_ctx->config.mla_config.mla_level;
    int cnt_num_level = 0;

    mla_ctx->metis_mla[cnt_num_level].fine = (DataGmsh *)malloc(sizeof(DataGmsh));
    mla_ctx->metis_mla[cnt_num_level].coarse = (DataGmsh *)malloc(sizeof(DataGmsh));
    assert(mla_ctx->metis_mla[cnt_num_level].fine &&
           mla_ctx->metis_mla[cnt_num_level].coarse);

    DeepCopyMetisDataGmsh(mla_ctx->metis_mla[cnt_num_level].fine, mla_ctx->data_gmsh);
#if 1
    puts("==== data of fine level mesh ====");
    printf("%d level, number of nodes: %d\n", cnt_num_level,
           mla_ctx->metis_mla[cnt_num_level].fine->nn);
    printf("%d level, number of elements: %d\n", cnt_num_level,
           mla_ctx->metis_mla[cnt_num_level].fine->ne_in);
    printf("%d level, number of partitions: %" PRIDX "\n", cnt_num_level,
           mla_ctx->metis_mla[cnt_num_level].fine->nparts);
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->nn; ++index)
    {
        printf("node %d:\t%021.16le\t%021.16le\t%021.16le\n", index,
               mla_ctx->metis_mla[cnt_num_level].fine->coordinates[3 * index],
               mla_ctx->metis_mla[cnt_num_level].fine->coordinates[3 * index + 1],
               mla_ctx->metis_mla[cnt_num_level].fine->coordinates[3 * index + 2]);
    }
    printf("number of elements in inner: %d\n", mla_ctx->metis_mla[cnt_num_level].fine->ne_in);
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->ne_in + 1; ++index)
    {
        printf("fine.eptr_in[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].fine->eptr_in[index]);
    }
    printf("nodes in inner elements:\n");
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].fine->ne_in; ++index)
    {
        printf("element %d: ", index);
        for (int index_i = 0; index_i < mla_ctx->metis_mla[cnt_num_level].fine->nne_in; ++index_i)
        {
            printf("%" PRIDX "\t", mla_ctx->metis_mla[cnt_num_level].fine->eind_in[index * mla_ctx->metis_mla[cnt_num_level].fine->nne_in + index_i]);
        }
        putchar('\n');
    }
#endif // print fine level data

    METIS_SetDefaultOptions(options);
    nn = mla_ctx->metis_mla[cnt_num_level].fine->nn;
    ne = mla_ctx->metis_mla[cnt_num_level].fine->ne_in;
    status_metis_part = METIS_PartMeshNodal(&ne, &nn,
                                            mla_ctx->metis_mla[cnt_num_level].fine->eptr_in,
                                            mla_ctx->metis_mla[cnt_num_level].fine->eind_in,
                                            NULL, NULL,
                                            &(mla_ctx->metis_mla[cnt_num_level].fine->nparts),
                                            NULL,
                                            options, &objval,
                                            mla_ctx->metis_mla[cnt_num_level].fine->epart_in,
                                            mla_ctx->metis_mla[cnt_num_level].fine->npart_in);
#if 1
    printf("metis status_metis_part = %d\n", status_metis_part);
    printf("partition objective value = %" PRIDX "\n", objval);
    puts("fine level node partition:");
    for (int index = 0; index < nn; ++index)
    {
        printf("npart[%d] = %" PRIDX "\n", index, mla_ctx->metis_mla[cnt_num_level].fine->npart_in[index]);
    }
#endif

    // coarse level data gmsh
    GmshCoarseLevelGenerator(mla_ctx->metis_mla[cnt_num_level].coarse,
                             mla_ctx->metis_mla[cnt_num_level].fine);
#if 1
    puts("\ncoarse level data:");
    puts("nodes:");
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
    {
        printf("nodes[%d]\t%021.16le\t%021.16le\t%021.16le\n", index,
               mla_ctx->metis_mla[cnt_num_level].coarse->coordinates[3 * index],
               mla_ctx->metis_mla[cnt_num_level].coarse->coordinates[3 * index + 1],
               mla_ctx->metis_mla[cnt_num_level].coarse->coordinates[3 * index + 2]);
    }
    puts("adjacency:");
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
    {
        idx_t index_start = mla_ctx->metis_mla[cnt_num_level].coarse->eptr_in[index];
        idx_t index_end = mla_ctx->metis_mla[cnt_num_level].coarse->eptr_in[index + 1];
        // printf("node %d: ", index);
        printf("%d\t", index);
        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            printf("%" PRIDX "\t", mla_ctx->metis_mla[cnt_num_level].coarse->eind_in[index_i]);
        }
        putchar('\n');
    }
#endif // print coarse level data

    // coarse level partition
    idx_t nvtxs = mla_ctx->metis_mla[cnt_num_level].coarse->nn, ncon = 1;
    METIS_SetDefaultOptions(options);
    status_metis_part = METIS_PartGraphKway(&nvtxs, &ncon,
                                            mla_ctx->metis_mla[cnt_num_level].coarse->eptr_in,
                                            mla_ctx->metis_mla[cnt_num_level].coarse->eind_in,
                                            NULL, NULL, NULL,
                                            &(mla_ctx->metis_mla[cnt_num_level].coarse->nparts),
                                            NULL, NULL, options,
                                            &objval, mla_ctx->metis_mla[cnt_num_level].coarse->npart_in);
#if 1
    puts("==== coarse level partition");
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
    {
        printf("npart[%d] = %" PRIDX "\n", index,
               mla_ctx->metis_mla[cnt_num_level].coarse->npart_in[index]);
    }
#endif // print coarse level partition

    for (cnt_num_level = 1; cnt_num_level < config_num_level; ++cnt_num_level)
    {
        if (mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn < 100)
        {
            printf("nodes in level %d is %d, less than 100, setup done!\n", cnt_num_level,
                   mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn);
            break;
        }

        printf("in level %d:\n", cnt_num_level);
        mla_ctx->metis_mla[cnt_num_level].fine = (DataGmsh *)malloc(sizeof(DataGmsh));
        mla_ctx->metis_mla[cnt_num_level].coarse = (DataGmsh *)malloc(sizeof(DataGmsh));
        assert(mla_ctx->metis_mla[cnt_num_level].fine &&
               mla_ctx->metis_mla[cnt_num_level].coarse);

        DeepCopyKLevelMetisDataMesh(mla_ctx->metis_mla[cnt_num_level].fine,
                                    mla_ctx->metis_mla[cnt_num_level - 1].coarse);

        // coarse level data metis adjacency
        MetisKLevelCoarseLevelGenerator(mla_ctx->metis_mla[cnt_num_level].coarse,
                                        mla_ctx->metis_mla[cnt_num_level].fine);

        // coarse level partition
        nvtxs = mla_ctx->metis_mla[cnt_num_level].coarse->nn;
        ncon = 1;
        METIS_SetDefaultOptions(options);
        status_metis_part = METIS_PartGraphKway(&nvtxs, &ncon,
                                                mla_ctx->metis_mla[cnt_num_level].coarse->eptr_in,
                                                mla_ctx->metis_mla[cnt_num_level].coarse->eind_in,
                                                NULL, NULL, NULL,
                                                &(mla_ctx->metis_mla[cnt_num_level].coarse->nparts),
                                                NULL, NULL, options,
                                                &objval, mla_ctx->metis_mla[cnt_num_level].coarse->npart_in);
#if 1
        puts("==== coarse level partition");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
        {
            printf("npart[%d] = %" PRIDX "\n", index,
                   mla_ctx->metis_mla[cnt_num_level].coarse->npart_in[index]);
        }
#endif // print coarse level partition
    }
    mla_ctx->true_num_level = cnt_num_level;

    return 0;
}

int MetisMLASolver(MLAContext *mla_ctx,
                   int mla_phase)
{
    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        MetisMLASolverSetupPhase(mla_ctx);
    }

    // solve phase
    if (mla_phase == 1 || mla_phase == 2)
    {
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
            MetisMLASolverSolvePhase(mla_ctx);
            MLASolverRelativeResidual(&(mla_ctx->mysolver), &rela_resid);
            ++iter_cnt;
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));
        }
    }

#if 0

    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        /*
         * numer of levels
         * finest mesh
         * rbm order (prolongation operator constructor)
         */
        MLASolverSetupPhase(mysolver, graph, num_level, order_rbm, mla_ctx);
    }

    // solve phase
    if (mla_phase == 1 || mla_phase == 2)
    {
        int cnt = 0;
        double rela_resid = 0.;
        MLASolverRelativeResidual(mysolver, &rela_resid);

        while (cnt < config->mla_config.mla_max_it && rela_resid > config->mla_config.mla_rtol)
        {
#if 0
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));
#endif // print mla iteration information
            /*
             * solve phase recursive implementation
             *     linear system
             *     mla context, contains setup information
             *     config, pre- and post- smooth times
             *     a special case, rbm order is 2 and level is 1, coarse operator need shift
             *     finest level and then recursive implementation
             */
            MLASolverSolvePhase(config, mla_ctx, order_rbm, mysolver);

            MLASolverRelativeResidual(mysolver, &rela_resid);
            ++cnt;
        }
#if 0
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));
#endif // print mla iteration information
    }

#if 0
    puts("\n\n==== multilevel information ====");
    for (int index = 0; index < mla_ctx->num_level; ++index)
    {
        printf(">>>> in neighbouring level %d:\n", (mla_ctx->mla + index)->level);
        printf("there are %d vertices in fine mesh, %d vertices in coarse mesh!\n\n",
               (mla_ctx->mla + index)->fine->size,
               (mla_ctx->mla + index)->coarse->size);

        puts("\n\n======== prolongation operator ========");
        int m_pro = 0, n_pro = 0;
        int m_coar = 0, n_coar = 0;
        PetscCall(MatGetSize((mla_ctx->mla + index)->prolongation, &m_pro, &n_pro));
        PetscCall(MatGetSize((mla_ctx->mla + index)->operator_coarse, &m_coar, &n_coar));
        printf(">>>> level %d: prolongation size m = %d, n = %d, coarse size m = %d, n = %d\n",
               index, m_pro, n_pro, m_coar, n_coar);
        // PetscCall(MatView((mla_ctx->mla + index)->prolongation, PETSC_VIEWER_STDOUT_WORLD));

        puts("\n\n======== coarse operator ========");
        printf(">>>> level %d: prolongation size m = %d, n = %d, coarse size m = %d, n = %d\n",
               index, m_pro, n_pro, m_coar, n_coar);
        // PetscCall(MatView((mla_ctx->mla + index)->operator_coarse, PETSC_VIEWER_STDOUT_WORLD));
    }
#endif // basic information of multilevel
#endif

    return 0;
}
