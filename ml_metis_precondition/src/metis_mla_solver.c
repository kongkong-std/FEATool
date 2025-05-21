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

int MetisMLANestedProcedurePreSmooth(KSP ksp, PC pc,
                                     int level,
                                     MLAContext *mla_ctx,
                                     Vec *mg_recur_x,
                                     Vec *mg_recur_b,
                                     int v_pre_smooth)
{
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
}

int MetisMLANestedProcedurePostSmooth(KSP ksp, PC pc,
                                      int level,
                                      MLAContext *mla_ctx,
                                      Vec *mg_recur_x,
                                      Vec *mg_recur_b,
                                      int v_post_smooth)
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
    PetscCall(PCSORSetIterations(pc, v_post_smooth, v_post_smooth));
    PetscCall(PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP));
    PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
    // PetscCall(KSPSetFromOptions(ksp));

    PetscCall(KSPSetTolerances(ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
    PetscCall(KSPSolve(ksp, mg_recur_b[level], mg_recur_x[level]));

    return 0;
}

int MetisMLASolverCoarsetCorrectionPhase(int order_rbm, KSP ksp, PC pc,
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

int MetisMLANestedProcedure(int level, int num_level,
                            MySolver *mysolver,
                            MLAContext *mla_ctx,
                            Vec *mg_recur_x,
                            Vec *mg_recur_b,
                            int v_pre_smooth,
                            int v_post_smooth,
                            int order_rbm)
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
        MetisMLANestedProcedurePreSmooth(mla_ctx->metis_mla[level].ksp_presmooth,
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
        PetscCall(VecCreate(PETSC_COMM_WORLD, mg_recur_b + level + 1));
        PetscCall(VecSetSizes(mg_recur_b[level + 1], PETSC_DECIDE, n_prolongation));
        PetscCall(VecSetFromOptions(mg_recur_b[level + 1]));
        PetscCall(MatMultTranspose(mla_ctx->metis_mla[level].prolongation, r_h[level], mg_recur_b[level + 1]));
    }

    // coarsest level
    MetisMLASolverCoarsetCorrectionPhase(order_rbm,
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
        MetisMLANestedProcedurePostSmooth(mla_ctx->metis_mla[level].ksp_postsmooth,
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

int MetisMLASolverSolvePhase(const ConfigJSON *config,
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

    MetisMLANestedProcedure(0, num_level,
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
    assert(coarse_data->eind_in);
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
        if (mla_ctx->metis_mla[cnt_num_level - 1].coarse->nn < 100)
        {
            printf("nodes in level %d (coarse) is %d, less than 100, setup done!\n", cnt_num_level - 1,
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
    mla_ctx->num_level = cnt_num_level;

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

                if (index_cnt_level == 0)
                {
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
                }
                else
                {
                    for (int index_fine_node = 0; index_fine_node < mla_ctx->metis_mla[index_cnt_level].fine->nn; ++index_fine_node)
                    {
                        double p_loc[9][9] = {0};
                        for (int index = 0; index < 9; ++index)
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
        // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));

        while (iter_cnt < mla_ctx->config.mla_config.mla_max_it &&
               rela_resid > mla_ctx->config.mla_config.mla_rtol)
        {
            MetisMLASolverSolvePhase(&(mla_ctx->config),
                                     mla_ctx,
                                     mla_ctx->order_rbm,
                                     &(mla_ctx->mysolver));
            MLASolverRelativeResidual(&(mla_ctx->mysolver), &rela_resid);
            ++iter_cnt;
            // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", iter_cnt, rela_resid));
        }

        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> solve time: %g (s)\n", time2 - time1));
    }
    return 0;
}
