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

int MetisMLANestedProcedurePreSmooth(int level /*current level*/,
                                     Vec *mla_recur_x /*solution in recursive*/,
                                     Vec *mla_recur_b /*rhs in recursive*/,
                                     MLAContext *mla_ctx /*mla context*/)
{
    int v_pre_smooth = mla_ctx->config.mla_config.pre_smooth_v;

    if (level == 0)
    {
        PetscCall(KSPSetType(mla_ctx->metis_mla[level].ksp_presmooth, KSPRICHARDSON));
        PetscCall(KSPGetPC(mla_ctx->metis_mla[level].ksp_presmooth,
                           &(mla_ctx->metis_mla[level].pc_presmooth)));
        PetscCall(PCSetType(mla_ctx->metis_mla[level].pc_presmooth, PCSOR));
        PetscCall(PCSORSetOmega(mla_ctx->metis_mla[level].pc_presmooth, 1.)); // gauss-seidel
        PetscCall(PCSORSetIterations(mla_ctx->metis_mla[level].pc_presmooth, v_pre_smooth, v_pre_smooth));
        PetscCall(PCSORSetSymmetric(mla_ctx->metis_mla[level].pc_presmooth, SOR_SYMMETRIC_SWEEP));
        PetscCall(KSPSetNormType(mla_ctx->metis_mla[level].ksp_presmooth, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetTolerances(mla_ctx->metis_mla[level].ksp_presmooth, 1e-10, 1e-10, PETSC_DEFAULT, 1));
        PetscCall(KSPSetInitialGuessNonzero(mla_ctx->metis_mla[level].ksp_presmooth, PETSC_TRUE));
        PetscCall(KSPSolve(mla_ctx->metis_mla[level].ksp_presmooth, mla_recur_b[level], mla_recur_x[level]));
    }
    else
    {
        PetscCall(VecDuplicate(mla_recur_b[level], mla_recur_x + level));
        PetscCall(KSPSetType(mla_ctx->metis_mla[level].ksp_presmooth, KSPRICHARDSON));
        PetscCall(KSPGetPC(mla_ctx->metis_mla[level].ksp_presmooth,
                           &(mla_ctx->metis_mla[level].pc_presmooth)));
        PetscCall(PCSetType(mla_ctx->metis_mla[level].pc_presmooth, PCSOR));
        PetscCall(PCSORSetOmega(mla_ctx->metis_mla[level].pc_presmooth, 1.)); // gauss-seidel
        PetscCall(PCSORSetIterations(mla_ctx->metis_mla[level].pc_presmooth, v_pre_smooth, v_pre_smooth));
        PetscCall(PCSORSetSymmetric(mla_ctx->metis_mla[level].pc_presmooth, SOR_SYMMETRIC_SWEEP));
        PetscCall(KSPSetNormType(mla_ctx->metis_mla[level].ksp_presmooth, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetTolerances(mla_ctx->metis_mla[level].ksp_presmooth, 1e-10, 1e-10, PETSC_DEFAULT, 1));
        PetscCall(KSPSetInitialGuessNonzero(mla_ctx->metis_mla[level].ksp_presmooth, PETSC_TRUE));
        PetscCall(KSPSolve(mla_ctx->metis_mla[level].ksp_presmooth, mla_recur_b[level], mla_recur_x[level]));
    }

    return 0;
}

int MetisMLANestedProcedurePostSmooth(int level /*current level*/,
                                      Vec *mla_recur_x /*solution in recursive*/,
                                      Vec *mla_recur_b /*rhs in recursive*/,
                                      MLAContext *mla_ctx /*mla context*/)
{
    int v_post_smooth = mla_ctx->config.mla_config.post_smooth_v;

    PetscCall(KSPSetType(mla_ctx->metis_mla[level].ksp_postsmooth, KSPRICHARDSON));
    PetscCall(KSPGetPC(mla_ctx->metis_mla[level].ksp_postsmooth,
                       &(mla_ctx->metis_mla[level].pc_postsmooth)));
    PetscCall(PCSetType(mla_ctx->metis_mla[level].pc_postsmooth, PCSOR));
    PetscCall(PCSORSetOmega(mla_ctx->metis_mla[level].pc_postsmooth, 1.)); // gauss-seidel
    PetscCall(PCSORSetIterations(mla_ctx->metis_mla[level].pc_postsmooth, v_post_smooth, v_post_smooth));
    PetscCall(PCSORSetSymmetric(mla_ctx->metis_mla[level].pc_postsmooth, SOR_SYMMETRIC_SWEEP));
    PetscCall(KSPSetNormType(mla_ctx->metis_mla[level].ksp_postsmooth, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetTolerances(mla_ctx->metis_mla[level].ksp_postsmooth, 1e-10, 1e-10, PETSC_DEFAULT, 1));
    PetscCall(KSPSetInitialGuessNonzero(mla_ctx->metis_mla[level].ksp_postsmooth, PETSC_TRUE));
    PetscCall(KSPSolve(mla_ctx->metis_mla[level].ksp_postsmooth, mla_recur_b[level], mla_recur_x[level]));

    return 0;
}

int MetisMLANestedProcedureCoarsestCorrection(int level /*current level*/,
                                              Vec *mla_recur_x /*solution in recursive*/,
                                              Vec *mla_recur_b /*rhs in recursive*/,
                                              MLAContext *mla_ctx /*mla context*/)
{
    PetscCall(VecDuplicate(mla_recur_b[level + 1], mla_recur_x + level + 1));
    int nrestart = mla_ctx->config.mla_config.coarse_restart;

    PetscCall(KSPSetType(mla_ctx->metis_mla[level].ksp_coarse, KSPGCR));
    PetscCall(KSPGCRSetRestart(mla_ctx->metis_mla[level].ksp_coarse, nrestart));
    PetscCall(KSPSetNormType(mla_ctx->metis_mla[level].ksp_coarse, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetTolerances(mla_ctx->metis_mla[level].ksp_coarse, 1e-10, 1e-10, PETSC_DEFAULT, 1000));
    PetscCall(KSPSolve(mla_ctx->metis_mla[level].ksp_coarse,
                       mla_recur_b[level + 1],
                       mla_recur_x[level + 1]));

    return 0;
}

int MetisMLANestedProcedure(Vec *mla_recur_x /*solution in recursive*/,
                            Vec *mla_recur_b /*rhs in recursive*/,
                            MLAContext *mla_ctx /*mla context*/)
{
    int num_level = mla_ctx->true_num_level;
    Vec *r_h = NULL, *tmp_e_h = NULL;

    r_h = (Vec *)malloc(num_level * sizeof(Vec));
    tmp_e_h = (Vec *)malloc(num_level * sizeof(Vec));
    assert(r_h && tmp_e_h);

    // loop implementation
    /*
     * downward sweep, from fine to coarse level
     */
    for (int level = 0; level < num_level; ++level)
    {
        PetscCall(VecDuplicate(mla_recur_b[level], r_h + level));

        // pre-smooth
        MetisMLANestedProcedurePreSmooth(level, mla_recur_x, mla_recur_b, mla_ctx);

        PetscCall(MatMult(mla_ctx->metis_mla[level].operator_fine,
                          mla_recur_x[level],
                          r_h[level]));
        PetscCall(VecAXPY(r_h[level], -1., mla_recur_b[level]));

        // restriction
        PetscInt m_prolongation = 0, n_prolongation = 0;
        PetscCall(MatGetSize(mla_ctx->metis_mla[level].prolongation,
                             &m_prolongation,
                             &n_prolongation));
        PetscCall(VecCreate(PETSC_COMM_WORLD, mla_recur_b + level + 1));
        PetscCall(VecSetSizes(mla_recur_b[level + 1], PETSC_DECIDE, n_prolongation));
        PetscCall(VecSetFromOptions(mla_recur_b[level + 1]));
        PetscCall(MatMultTranspose(mla_ctx->metis_mla[level].prolongation,
                                   r_h[level],
                                   mla_recur_b[level + 1]));
    }

    /*
     * coarset level
     */
    // PetscCall(VecSet(mla_recur_x[num_level], 0.));
#if 1
    MetisMLANestedProcedureCoarsestCorrection(num_level - 1,
                                              mla_recur_x,
                                              mla_recur_b,
                                              mla_ctx);
#endif

    /*
     * upward sweep, from coarse to fine level
     */
    for (int level = num_level - 1; level >= 0; --level)
    {
        PetscCall(VecDuplicate(mla_recur_x[level], tmp_e_h + level));

        // prolongation
        PetscCall(MatMult(mla_ctx->metis_mla[level].prolongation, mla_recur_x[level + 1], tmp_e_h[level]));
        PetscCall(VecAXPY(mla_recur_x[level], 1., tmp_e_h[level]));

        // post-smooth
        MetisMLANestedProcedurePostSmooth(level,
                                          mla_recur_x,
                                          mla_recur_b,
                                          mla_ctx);
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

int MetisMLASolverSolvePhase(MLAContext *mla_ctx)
{
    int num_level = mla_ctx->true_num_level;
    Vec *mla_recur_x, *mla_recur_b;

    mla_recur_x = (Vec *)malloc((num_level + 1) * sizeof(Vec));
    mla_recur_b = (Vec *)malloc((num_level + 1) * sizeof(Vec));
    assert(mla_recur_x && mla_recur_b);

    PetscCall(VecDuplicate(mla_ctx->mysolver.solver_x, mla_recur_x));
    PetscCall(VecDuplicate(mla_ctx->mysolver.solver_b, mla_recur_b));
    PetscCall(VecCopy(mla_ctx->mysolver.solver_x, mla_recur_x[0]));
    PetscCall(VecCopy(mla_ctx->mysolver.solver_b, mla_recur_b[0]));

    MetisMLANestedProcedure(mla_recur_x, mla_recur_b, mla_ctx);

    // updating solution
    PetscCall(VecCopy(mla_recur_x[0], mla_ctx->mysolver.solver_x));

    // free memory
    free(mla_recur_x);
    free(mla_recur_b);

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
#if 0
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
#if 0
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
#if 0
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
#if 0
    puts("==== coarse level partition");
    for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
    {
        printf("npart[%d] = %" PRIDX "\n", index,
               mla_ctx->metis_mla[cnt_num_level].coarse->npart_in[index]);
    }
#endif // print coarse level partition

    for (cnt_num_level = 1; cnt_num_level < config_num_level; ++cnt_num_level)
    {
        if (mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn < 1000)
        {
            printf("nodes in level %d (coarse) is %d, less than 1000, setup done!\n", cnt_num_level - 1,
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
#if 0
        puts("==== coarse level partition");
        for (int index = 0; index < mla_ctx->metis_mla[cnt_num_level].coarse->nn; ++index)
        {
            printf("npart[%d] = %" PRIDX "\n", index,
                   mla_ctx->metis_mla[cnt_num_level].coarse->npart_in[index]);
        }
#endif // print coarse level partition
    }
    mla_ctx->true_num_level = cnt_num_level;

    // operator assembling
    PetscInt mat_m = 0, mat_n = 0;
    PetscCall(MatGetSize(mla_ctx->mysolver.solver_a, &mat_m, &mat_n));
#if 0
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix size: m = %" PetscInt_FMT ", n = %" PetscInt_FMT "\n", mat_m, mat_n));
#endif // print matrix size information

    PetscCall(MatDuplicate(mla_ctx->mysolver.solver_a,
                           MAT_COPY_VALUES,
                           &(mla_ctx->metis_mla[0].operator_fine)));

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
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      PETSC_DECIDE, PETSC_DECIDE,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = 0; index_fine_node < mla_ctx->metis_mla[index_cnt_level].fine->nn; ++index_fine_node)
                {
                    double p_loc[6][6] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->npart_in[index_fine_node];
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
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      PETSC_DECIDE, PETSC_DECIDE,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = 0; index_fine_node < mla_ctx->metis_mla[index_cnt_level].fine->nn; ++index_fine_node)
                {
                    double p_loc[6][9] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->npart_in[index_fine_node];
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
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      PETSC_DECIDE, PETSC_DECIDE,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = 0; index_fine_node < mla_ctx->metis_mla[index_cnt_level].fine->nn; ++index_fine_node)
                {
                    double p_loc[6][6] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->npart_in[index_fine_node];
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
                PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[index_cnt_level].prolongation)));
                PetscCall(MatSetSizes(mla_ctx->metis_mla[index_cnt_level].prolongation,
                                      PETSC_DECIDE, PETSC_DECIDE,
                                      prolongation_m, prolongation_n));
                PetscCall(MatSetUp(mla_ctx->metis_mla[index_cnt_level].prolongation));

                for (int index_fine_node = 0; index_fine_node < mla_ctx->metis_mla[index_cnt_level].fine->nn; ++index_fine_node)
                {
                    double p_loc[6][9] = {0};
                    for (int index = 0; index < 6; ++index)
                    {
                        p_loc[index][index] = 1.;
                    }

                    idx_t index_coarse_node = mla_ctx->metis_mla[index_cnt_level].fine->npart_in[index_fine_node];
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

    for (int level = 0; level < cnt_num_level; ++level)
    {
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_presmooth)));
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_presmooth,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));

        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_postsmooth)));
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_postsmooth,
                                  mla_ctx->metis_mla[level].operator_fine,
                                  mla_ctx->metis_mla[level].operator_fine));

        PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx->metis_mla[level].ksp_coarse)));
        PetscCall(KSPSetOperators(mla_ctx->metis_mla[level].ksp_coarse,
                                  mla_ctx->metis_mla[level].operator_coarse,
                                  mla_ctx->metis_mla[level].operator_coarse));
    }

    return 0;
}

int MetisMLASolver(MLAContext *mla_ctx,
                   int mla_phase)
{
    PetscLogDouble time1, time2;

    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        PetscCall(PetscTime(&time1));
        MetisMLASolverSetupPhase(mla_ctx);
        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> setup time: %g (s)\n", time2 - time1););
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

    return 0;
}
