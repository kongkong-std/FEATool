#include "../include/main.h"

int DeepCopyMetisDataGmsh(DataGmsh *dst, DataGmsh *src);

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

#if 0
    puts("\n==== metis function result ====");
    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    printf("partition objective value = %ld\n", objval);
    printf("element partition:\n");
    for (int index = 0; index < ne; ++index)
    {
        printf("epart[%d] = %ld\n", index, data.epart_in[index]);
    }
    printf("node partition:\n");
    for (int index = 0; index < nn; ++index)
    {
        // printf("npart[%d] = %ld\n", index, npart[index]);
        printf("%ld\n", data.npart_in[index]);
    }
#endif // metis function result

#if 0
    // free memory
    free(epart);
    free(npart);
#endif

    return status;
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

int MetisMLASolverSetupPhase(MLAContext *mla_ctx)
{
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
    printf("metis status_metis_part = %d\n", status_metis_part);

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
