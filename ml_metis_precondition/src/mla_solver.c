#include "../include/main.h"

int MLAMGNestedProcedurePreSmooth(KSP ksp, PC pc,
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
                                  (mla_ctx->mla + level)->operator_fine,
                                  (mla_ctx->mla + level)->operator_fine));
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
        double shift = 1e-10; // test metis aggregation implementation
        // double shift = 0.; // test metis aggregation implementation
        PetscCall(MatShift((mla_ctx->mla + level)->operator_fine, shift));
#endif
        PetscCall(KSPSetOperators(ksp,
                                  (mla_ctx->mla + level)->operator_fine,
                                  (mla_ctx->mla + level)->operator_fine)); // add shift to diagonal
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

int MLAMGNestedProcedurePostSmooth(KSP ksp, PC pc,
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
                              (mla_ctx->mla + level)->operator_fine,
                              (mla_ctx->mla + level)->operator_fine));
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

int MLASolverCoarsetCorrectionPhase(int order_rbm, KSP ksp, PC pc,
                                    int level,
                                    MLAContext *mla_ctx,
                                    Vec *mg_recur_x,
                                    Vec *mg_recur_b)
{
    PetscCall(VecDuplicate(mg_recur_b[level + 1], mg_recur_x + level + 1));
    int gcr_restart = 50;

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
                                      (mla_ctx->mla + level)->operator_coarse,
                                      (mla_ctx->mla + level)->operator_coarse));
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
            // double shift = 1e-12;     // test metis aggregation implementation
            double shift = 0.; // test metis aggregation implementation
            PetscCall(MatShift((mla_ctx->mla + level)->operator_coarse, shift));
#endif
            PetscCall(KSPSetOperators(ksp,
                                      (mla_ctx->mla + level)->operator_coarse,
                                      (mla_ctx->mla + level)->operator_coarse)); // add shift to diagonal
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
        double shift = 1e-12; // test metis aggregation implementation
        // double shift = 0.; // test metis aggregation implementation
        PetscCall(MatShift((mla_ctx->mla + level)->operator_coarse, shift));
#endif
        PetscCall(KSPSetOperators(ksp,
                                  (mla_ctx->mla + level)->operator_coarse,
                                  (mla_ctx->mla + level)->operator_coarse)); // add shift to diagonal
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

int MLAMGNestedProcedure(int level, int num_level,
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
        MLAMGNestedProcedurePreSmooth((mla_ctx->mla + level)->ksp_presmooth,
                                      (mla_ctx->mla + level)->pc_presmooth,
                                      level,
                                      mla_ctx,
                                      mg_recur_x,
                                      mg_recur_b,
                                      v_pre_smooth);
        PetscCall(MatMult((mla_ctx->mla + level)->operator_fine,
                          mg_recur_x[level],
                          r_h[level]));
        PetscCall(VecAYPX(r_h[level], -1., mg_recur_b[level]));

        // restriction
        PetscInt m_prolongation = 0, n_prolongation = 0; // size of prolongation operator
        PetscCall(MatGetSize((mla_ctx->mla + level)->prolongation,
                             &m_prolongation,
                             &n_prolongation));
        PetscCall(VecCreate(PETSC_COMM_WORLD, mg_recur_b + level + 1));
        PetscCall(VecSetSizes(mg_recur_b[level + 1], PETSC_DECIDE, n_prolongation));
        PetscCall(VecSetFromOptions(mg_recur_b[level + 1]));
        PetscCall(MatMultTranspose((mla_ctx->mla + level)->prolongation, r_h[level], mg_recur_b[level + 1]));
    }

    // coarsest level
    MLASolverCoarsetCorrectionPhase(order_rbm,
                                    (mla_ctx->mla + level - 1)->ksp_coarse,
                                    (mla_ctx->mla + level - 1)->pc_coarse,
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

        PetscCall(MatMult((mla_ctx->mla + level)->prolongation, mg_recur_x[level + 1], tmp_e_h[level]));
        PetscCall(VecAXPY(mg_recur_x[level], 1., tmp_e_h[level]));

        // post-smooth procedure
        MLAMGNestedProcedurePostSmooth((mla_ctx->mla + level)->ksp_postsmooth,
                                       (mla_ctx->mla + level)->pc_postsmooth,
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

int MLASolverSolvePhase(const ConfigJSON *config,
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

    MLAMGNestedProcedure(0, num_level,
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

int MLASolverSetupPhase(MySolver *mysolver,
                        const MeshGraph *graph,
                        int num_level,
                        int order_rbm,
                        MLAContext *mla_ctx)
{
    if (mla_ctx->setup == 1)
    {
        // prolongation operator has been constructed
        PetscFunctionReturn(PETSC_SUCCESS);
    }

#if 0
    MetisMLASolver(mla_ctx, 0);
#endif // test metis aggregtion implementation, added

#if 1
    mla_ctx->mla = (MLAGraph *)malloc(num_level * sizeof(MLAGraph));
    assert(mla_ctx->mla);
#endif

    int cnt_level = 0;

    // mesh information
    (mla_ctx->mla + cnt_level)->level = cnt_level;
    (mla_ctx->mla + cnt_level)->mesh_tmp = CreateMeshGraph(graph->size);
    (mla_ctx->mla + cnt_level)->fine = CreateMeshGraph(graph->size);

    CopyMeshGraph((mla_ctx->mla + cnt_level)->mesh_tmp, graph);
    CopyMeshGraph((mla_ctx->mla + cnt_level)->fine, graph); // neighbouring fine mesh
#if 0
    printf("\n\n==== level %d neighbouring fine mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->fine);
#endif // fine mesh

    (mla_ctx->mla + cnt_level)->aggregation = AggregationMeshGraph((mla_ctx->mla + cnt_level)->mesh_tmp); // aggregation mesh
#if 0
    printf("\n\n==== level %d neighbouring aggregation mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->aggregation);
#endif // aggregation mesh

    (mla_ctx->mla + cnt_level)->coarse_node = (MeshNode *)malloc((mla_ctx->mla + cnt_level)->aggregation->size * sizeof(MeshNode));
    assert((mla_ctx->mla + cnt_level)->coarse_node);
    AssembleCoarseMeshNode((mla_ctx->mla + cnt_level)->aggregation,
                           (mla_ctx->mla + cnt_level)->coarse_node);
    (mla_ctx->mla + cnt_level)->coarse = CreateMeshGraph((mla_ctx->mla + cnt_level)->aggregation->size);
    AssembleCoarseMeshGraph((mla_ctx->mla + cnt_level)->coarse,
                            (mla_ctx->mla + cnt_level)->aggregation,
                            (mla_ctx->mla + cnt_level)->fine,
                            (mla_ctx->mla + cnt_level)->coarse_node); // neighbouring coarse mesh
#if 0
    printf("\n\n==== level %d neighbouring coarse mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->coarse);
#endif

#if 1
    PetscCall(MatDuplicate(mysolver->solver_a,
                           MAT_COPY_VALUES,
                           &((mla_ctx->mla + cnt_level)->operator_fine)));
#endif // test metis aggregtion implementation

#if 0
    PetscCall(MatDuplicate(mla_ctx->metis_mla[cnt_level].operator_fine,
                           MAT_COPY_VALUES,
                           &((mla_ctx->mla + cnt_level)->operator_fine)));
#endif // test metis aggregation implementation, added

#if 0
    printf("\n\n==== level %d, neighbouring fine operator:\n", cnt_level);
    PetscCall(MatView((mla_ctx->mla + cnt_level)->operator_fine, PETSC_VIEWER_STDOUT_WORLD));
#endif // print neighbouring fine operator

#if 1
    int prolongation_block_row = (mla_ctx->mla + cnt_level)->fine->size;   // vertex in fine mesh
    int prolongation_block_col = (mla_ctx->mla + cnt_level)->coarse->size; // vertex in coarse mesh
    int m_prolongation = prolongation_block_row * 6;                       // row of prolongation operator
    int n_prolongation = 0;

    if (order_rbm == 1)
    {
        // 6 dof in coarse point
        n_prolongation = prolongation_block_col * 6;
    }
    else if (order_rbm == 2)
    {
        // 9 dof in coarse point
        n_prolongation = prolongation_block_col * 9;
    }

#if 1
    printf(">>>>>>>> setup phase:\n prolongation_block_row = %d, col = %d\n",
           prolongation_block_row,
           prolongation_block_col);
    printf("size of prolongation: m = %d, n = %d\n", m_prolongation, n_prolongation);
#endif

    // memory allocation for prolongation operator
    PetscCall(MatCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->prolongation)));
    PetscCall(MatSetSizes((mla_ctx->mla + cnt_level)->prolongation, PETSC_DECIDE, PETSC_DECIDE, m_prolongation, n_prolongation));
    PetscCall(MatSetUp((mla_ctx->mla + cnt_level)->prolongation));

    // assigning value to prolongation operator
    double **p_loc = (double **)malloc(6 * sizeof(double *));
    assert(p_loc);
    memset(p_loc, 0, 6 * sizeof(double *));

    for (int index = 0; index < 6; ++index)
    {
        if (order_rbm == 1)
        {
            p_loc[index] = (double *)malloc(6 * sizeof(double));
            memset(p_loc[index], 0, 6 * sizeof(double));
        }
        else if (order_rbm == 2)
        {
            p_loc[index] = (double *)malloc(9 * sizeof(double));
            memset(p_loc[index], 0, 9 * sizeof(double));
        }
        else
        {
            fprintf(stderr, "invalid order_rbm 1 or 2!!!\n");
            exit(EXIT_FAILURE);
        }
    }

    // lame constant
    double lame_e = 21e5, lame_nu = 0.28;
    double lame_lambda = 0., lame_mu = 0.;
    double lame_constant_modify = 0.;
    lame_lambda = lame_e * lame_nu / (1 + lame_nu) / (1 - 2 * lame_nu);
    lame_mu = lame_e / 2. / (1 + lame_nu);
    lame_constant_modify = lame_mu / (lame_lambda + lame_mu);
    lame_constant_modify = 1.;

    // aggregation graph traverse
    for (int index = 0; index < (mla_ctx->mla + cnt_level)->aggregation->size; ++index)
    {
        MeshGraphAdjNode *data_node_coarse = (mla_ctx->mla + cnt_level)->coarse->array[index].head;
        MeshGraphAdjNode *data_node_aggregation = (mla_ctx->mla + cnt_level)->aggregation->array[index].head;

        // coarse node coordinate
        double node_coarse_x = data_node_coarse->node->node_x;
        double node_coarse_y = data_node_coarse->node->node_y;
        double node_coarse_z = data_node_coarse->node->node_z;

        // aggregation graph adjacency list
        for (int index_i = 0; index_i < (mla_ctx->mla + cnt_level)->aggregation->array[index].size; ++index_i)
        {
            // 0-base
            int m_prolongation_block_tmp = data_node_aggregation->node->node_idx - 1;
            int n_prolongation_block_tmp = index;

            double node_aggregation_x = data_node_aggregation->node->node_x;
            double node_aggregation_y = data_node_aggregation->node->node_y;
            double node_aggregation_z = data_node_aggregation->node->node_z;

            if (order_rbm == 1)
            {
                for (int index_tmp = 0; index_tmp < 6; ++index_tmp)
                {
                    p_loc[index_tmp][index_tmp] = 1.;
                }
#if 0
                p_loc[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                p_loc[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                p_loc[1][3] = -p_loc[0][4];                       // (2, 4) -z
                p_loc[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                p_loc[2][3] = -p_loc[0][5];                       // (3, 4) y
                p_loc[2][4] = -p_loc[1][5];                       // (3, 5) -x
#endif
#if 1
                p_loc[0][3] = -(node_coarse_z - node_aggregation_z); // (1, 5) z, phi modified
                p_loc[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y, phi modified
                p_loc[1][4] = p_loc[0][4];                           // (2, 5) z, phi modified
                p_loc[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x, phi modified
                p_loc[2][3] = -p_loc[1][5];                          // (3, 4) -x, phi modified
                p_loc[2][4] = p_loc[0][5];                           // (3, 5) -y, phi modified
#endif

                // p_loc in prolongation operator location
                /*
                 * (m_tmp, n_tmp) <-> {(m_tmp + 1) * 6 - 1, (n_tmp + 1) * 6 - 1}
                 */
                for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                {
                    for (int index_tmp_j = 0; index_tmp_j < 6; ++index_tmp_j)
                    {
                        PetscCall(MatSetValue((mla_ctx->mla + cnt_level)->prolongation,
                                              m_prolongation_block_tmp * 6 + index_tmp_i,
                                              n_prolongation_block_tmp * 6 + index_tmp_j,
                                              p_loc[index_tmp_i][index_tmp_j],
                                              INSERT_VALUES));
                    }
                }
            }
            else if (order_rbm == 2)
            {
                for (int index_tmp = 0; index_tmp < 6; ++index_tmp)
                {
                    p_loc[index_tmp][index_tmp] = 1.;
                }
#if 0
                p_loc[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                p_loc[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                p_loc[1][3] = -p_loc[0][4];                       // (2, 4) -z
                p_loc[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                p_loc[2][3] = -p_loc[0][5];                       // (3, 4) y
                p_loc[2][4] = -p_loc[1][5];                       // (3, 5) -x
#endif
#if 0
                p_loc[0][4] = -(node_coarse_z - node_aggregation_z); // (1, 5) z
                p_loc[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y
                p_loc[1][3] = -p_loc[0][4];                          // (2, 4) -z
                p_loc[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x
                p_loc[2][3] = -p_loc[0][5];                          // (3, 4) y
                p_loc[2][4] = -p_loc[1][5];                          // (3, 5) -x
#endif
#if 1
                p_loc[0][3] = -(node_coarse_z - node_aggregation_z); // (1, 4) z, phi modified
                p_loc[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y, phi modified
                p_loc[1][4] = p_loc[0][3];                           // (2, 5) z, phi modified
                p_loc[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x, phi modified
                p_loc[2][3] = -p_loc[1][5];                          // (3, 4) -x, phi modified
                p_loc[2][4] = p_loc[0][5];                           // (3, 5) -y, phi modified
#endif

#if 0
                p_loc[0][6] = (node_aggregation_z - node_coarse_z) *
                              (node_coarse_x - node_aggregation_x); // (1, 7) -zx
#endif
#if 0
                p_loc[0][6] = (node_coarse_z - node_aggregation_z) *
                              (node_coarse_x - node_aggregation_x); // (1, 7) zx, modified
#endif
#if 0
                p_loc[0][6] = (node_coarse_z - node_aggregation_z) *
                              (node_coarse_x - node_aggregation_x); // (1, 7) zx, phi modified
#endif

#if 0 
                p_loc[0][8] = (node_aggregation_z - node_coarse_z) *
                              (node_coarse_y - node_aggregation_y); // (1, 9) -zy
#endif
#if 0
                p_loc[0][8] = -(node_aggregation_z - node_coarse_z) *
                              (node_coarse_y - node_aggregation_y) / 2.; // (1, 9) zy/2, modified
#endif
#if 0
                p_loc[0][8] = -(node_aggregation_z - node_coarse_z) *
                              (node_coarse_y - node_aggregation_y) / 2.; // (1, 9) zy/2, phi modified
#endif

#if 0
                p_loc[1][7] = p_loc[0][8]; // (2, 8) -zy
                p_loc[1][8] = p_loc[0][6]; // (2, 9) -zx
#endif
#if 0
                p_loc[1][7] = 2. * p_loc[0][8]; // (2, 8) zy, modified
                p_loc[1][8] = p_loc[0][6] / 2.; // (2, 9) zx/2, modified
#endif
#if 0
                p_loc[1][7] = 2. * p_loc[0][8]; // (2, 8) zy, phi modified
                p_loc[1][8] = p_loc[0][6] / 2.; // (2, 9) zx/2, phi modified
#endif

#if 0
                p_loc[2][6] = (node_coarse_x - node_aggregation_x) *
                              (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) xx/2
#endif
#if 0
                p_loc[2][6] = -(node_coarse_x - node_aggregation_x) *
                              (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) -xx/2, modified
                p_loc[2][7] = -(node_coarse_y - node_aggregation_y) *
                              (node_coarse_y - node_aggregation_y) / 2.; // (3, 8) -yy/2, modified
                p_loc[2][8] = -(node_coarse_x - node_aggregation_x) *
                              (node_coarse_y - node_aggregation_y) / 2.; // (3, 9) -xy/2, modified
#endif

#if 0 // lame constant modify
                p_loc[1][7] *= lame_constant_modify;
                p_loc[2][7] *= lame_constant_modify;
#endif

#if 0
                p_loc[3][7] = p_loc[0][5];       // (4, 8) -y, modified
                p_loc[3][8] = -p_loc[1][5] / 2.; // (4, 9) -x/2, modified
                p_loc[4][6] = -2. * p_loc[3][8]; // (5, 7) x, modified
                p_loc[4][8] = -p_loc[3][7] / 2.; // (5, 9) y/2, modified
#endif
#if 0
                p_loc[3][6] = p_loc[1][5];      // (4, 7) x, phi modified
                p_loc[3][8] = -p_loc[0][5] / 2.; // (4, 9) y/2, phi modified
                p_loc[4][7] = 2. * p_loc[3][8];  // (5, 8) y, phi modified
                p_loc[4][8] = p_loc[3][6] / 2.;  // (5, 9) x/2, phi modified
#endif

#if 1                                                                                                       // null space
                p_loc[0][7] = (node_coarse_x - node_aggregation_x) * (node_coarse_y - node_aggregation_y);  // (0, 7) xy
                p_loc[0][8] = (node_coarse_x - node_aggregation_x) * (node_coarse_z - node_aggregation_z);  // (0, 8) xz
                p_loc[1][6] = p_loc[0][7];                                                                  // (1, 6) xy
                p_loc[1][8] = -(node_coarse_y - node_aggregation_y) * (node_coarse_z - node_aggregation_z); // (1, 8) -yz
                p_loc[2][6] = -p_loc[0][8];                                                                 // (2, 6) -xz
                p_loc[2][7] = p_loc[1][8];                                                                  // (2, 7) -yz
#endif

                // p_loc in prolongation operator location
                /*
                 * (m_tmp, n_tmp) <-> {(m_tmp + 1) * 6 - 1, (n_tmp + 1) * 9 - 1}
                 */
                for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                {
                    for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                    {
                        PetscCall(MatSetValue((mla_ctx->mla + cnt_level)->prolongation,
                                              m_prolongation_block_tmp * 6 + index_tmp_i,
                                              n_prolongation_block_tmp * 9 + index_tmp_j,
                                              p_loc[index_tmp_i][index_tmp_j],
                                              INSERT_VALUES));
                    }
                }
            }

            data_node_aggregation = data_node_aggregation->next;
        }
    }
    PetscCall(MatAssemblyBegin((mla_ctx->mla + cnt_level)->prolongation, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((mla_ctx->mla + cnt_level)->prolongation, MAT_FINAL_ASSEMBLY));
#if 0
    printf("\n\n==== level %d, neighbouring prolongation operator:\n", cnt_level);
    PetscCall(MatView((mla_ctx->mla + cnt_level)->prolongation, PETSC_VIEWER_STDOUT_WORLD));
#endif // neighbouring prolongation operator

#endif // test metis aggregation implementation

#if 0
    PetscCall(MatDuplicate(mla_ctx->metis_mla[cnt_level].prolongation,
                           MAT_COPY_VALUES,
                           &((mla_ctx->mla + cnt_level)->prolongation)));
    puts("======== adjacency list setup, duplicate prolongation >>>>>>>>");
#endif // test metis aggregation implementation, added

    PetscCall(MatPtAP((mla_ctx->mla + cnt_level)->operator_fine,
                      (mla_ctx->mla + cnt_level)->prolongation,
                      MAT_INITIAL_MATRIX,
                      PETSC_DETERMINE,
                      &((mla_ctx->mla + cnt_level)->operator_coarse)));
    puts("======== adjacency list setup, duplicate coarse_operator >>>>>>>>");
#if 0
    printf("\n\n==== level %d, neighbouring coarse operator:\n", cnt_level);
    PetscCall(MatView((mla_ctx->mla + cnt_level)->operator_coarse, PETSC_VIEWER_STDOUT_WORLD));
#endif

    mla_ctx->setup = 1;

    // ksp and pc object
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_presmooth)));
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_postsmooth)));
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_coarse)));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_presmooth)));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_postsmooth)));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_coarse)));

    ++cnt_level;
#if 1
    while (cnt_level < num_level &&
           (mla_ctx->mla + cnt_level - 1)->coarse->size > 50)
#endif // test metis aggregation implementation
#if 0
    printf("==== before setup loop, cnt_level = %d, num_level = %d>>>>!!!!\n", cnt_level, num_level);
    while (cnt_level < num_level &&
           (mla_ctx->metis_mla[cnt_level - 1].coarse->nn > 100))
#endif // test metis aggregation implementation, added
    {
        printf("in initial loop, cnt_level = %d\n", cnt_level);
        // next level
        (mla_ctx->mla + cnt_level)->level = cnt_level;
        (mla_ctx->mla + cnt_level)->mesh_tmp = CreateMeshGraph((mla_ctx->mla + cnt_level - 1)->coarse->size);
        (mla_ctx->mla + cnt_level)->fine = CreateMeshGraph((mla_ctx->mla + cnt_level - 1)->coarse->size);

        CopyMeshGraph((mla_ctx->mla + cnt_level)->mesh_tmp, (mla_ctx->mla + cnt_level - 1)->coarse);
        CopyMeshGraph((mla_ctx->mla + cnt_level)->fine, (mla_ctx->mla + cnt_level - 1)->coarse); // neighbouring fine mesh
#if 0
    printf("\n\n==== level %d neighbouring fine mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->fine);
#endif // fine mesh

        (mla_ctx->mla + cnt_level)->aggregation = AggregationMeshGraph((mla_ctx->mla + cnt_level)->mesh_tmp); // aggregation mesh
#if 0
    printf("\n\n==== level %d neighbouring aggregation mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->aggregation);
#endif // aggregation mesh

        (mla_ctx->mla + cnt_level)->coarse_node = (MeshNode *)malloc((mla_ctx->mla + cnt_level)->aggregation->size * sizeof(MeshNode));
        assert((mla_ctx->mla + cnt_level)->coarse_node);
        AssembleCoarseMeshNode((mla_ctx->mla + cnt_level)->aggregation,
                               (mla_ctx->mla + cnt_level)->coarse_node);
        (mla_ctx->mla + cnt_level)->coarse = CreateMeshGraph((mla_ctx->mla + cnt_level)->aggregation->size);
        AssembleCoarseMeshGraph((mla_ctx->mla + cnt_level)->coarse,
                                (mla_ctx->mla + cnt_level)->aggregation,
                                (mla_ctx->mla + cnt_level)->fine,
                                (mla_ctx->mla + cnt_level)->coarse_node); // neighbouring coarse mesh
#if 0
    printf("\n\n==== level %d neighbouring coarse mesh:\n", cnt_level);
    PrintMeshGraph((mla_ctx->mla + cnt_level)->coarse);
#endif

#if 1
        PetscCall(MatDuplicate((mla_ctx->mla + cnt_level - 1)->operator_coarse,
                               MAT_COPY_VALUES,
                               &((mla_ctx->mla + cnt_level)->operator_fine)));
#endif // test metis aggregation implementation

#if 0
        PetscCall(MatDuplicate(mla_ctx->metis_mla[cnt_level].operator_fine,
                               MAT_COPY_VALUES,
                               &((mla_ctx->mla + cnt_level)->operator_fine)));
#endif // test metis aggregation implementation, added

#if 0
    printf("\n\n==== level %d, neighbouring fine operator:\n", cnt_level);
    PetscCall(MatView((mla_ctx->mla + cnt_level)->operator_fine, PETSC_VIEWER_STDOUT_WORLD));
#endif // print neighbouring fine operator

#if 1
        prolongation_block_row = (mla_ctx->mla + cnt_level)->fine->size;   // vertex in fine mesh
        prolongation_block_col = (mla_ctx->mla + cnt_level)->coarse->size; // vertex in coarse mesh

        if (order_rbm == 1)
        {
            // 6 dof in coarse point
            m_prolongation = prolongation_block_row * 6;
            n_prolongation = prolongation_block_col * 6;
        }
        else if (order_rbm == 2)
        {
            // 9 dof in coarse point
            m_prolongation = prolongation_block_row * 9;
            n_prolongation = prolongation_block_col * 9;
        }

#if 1
        printf(">>>>>>>> setup phase:\n prolongation_block_row = %d, col = %d\n",
               prolongation_block_row,
               prolongation_block_col);
        printf("size of prolongation: m = %d, n = %d\n", m_prolongation, n_prolongation);
#endif

        // memory allocation for prolongation operator
        PetscCall(MatCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->prolongation)));
        PetscCall(MatSetSizes((mla_ctx->mla + cnt_level)->prolongation, PETSC_DECIDE, PETSC_DECIDE, m_prolongation, n_prolongation));
        PetscCall(MatSetUp((mla_ctx->mla + cnt_level)->prolongation));

        // assigning value to prolongation operator
        double **p_loc_tmp = NULL;
        if (order_rbm == 1)
        {
            p_loc_tmp = (double **)malloc(6 * sizeof(double *));
            assert(p_loc_tmp);
            memset(p_loc_tmp, 0, 6 * sizeof(double *));

            for (int index = 0; index < 6; ++index)
            {
                p_loc_tmp[index] = (double *)malloc(6 * sizeof(double));
                memset(p_loc_tmp[index], 0, 6 * sizeof(double));
            }
        }
        else if (order_rbm == 2)
        {
            p_loc_tmp = (double **)malloc(9 * sizeof(double *));
            assert(p_loc_tmp);
            memset(p_loc_tmp, 0, 9 * sizeof(double *));

            for (int index = 0; index < 9; ++index)
            {
                p_loc_tmp[index] = (double *)malloc(9 * sizeof(double));
                memset(p_loc_tmp[index], 0, 9 * sizeof(double));
            }
        }
        else
        {
            fprintf(stderr, "invalid order_rbm 1 or 2!!!\n");
            exit(EXIT_FAILURE);
        }

        // aggregation graph traverse
        for (int index = 0; index < (mla_ctx->mla + cnt_level)->aggregation->size; ++index)
        {
            MeshGraphAdjNode *data_node_coarse = (mla_ctx->mla + cnt_level)->coarse->array[index].head;
            MeshGraphAdjNode *data_node_aggregation = (mla_ctx->mla + cnt_level)->aggregation->array[index].head;

            // coarse node coordinate
            double node_coarse_x = data_node_coarse->node->node_x;
            double node_coarse_y = data_node_coarse->node->node_y;
            double node_coarse_z = data_node_coarse->node->node_z;

            // aggregation graph adjacency list
            for (int index_i = 0; index_i < (mla_ctx->mla + cnt_level)->aggregation->array[index].size; ++index_i)
            {
                // 0-base
                int m_prolongation_block_tmp = data_node_aggregation->node->node_idx - 1;
                int n_prolongation_block_tmp = index;

                double node_aggregation_x = data_node_aggregation->node->node_x;
                double node_aggregation_y = data_node_aggregation->node->node_y;
                double node_aggregation_z = data_node_aggregation->node->node_z;

                if (order_rbm == 1)
                {
                    for (int index_tmp = 0; index_tmp < 6; ++index_tmp)
                    {
                        p_loc[index_tmp][index_tmp] = 1.;
                    }
#if 0
                    p_loc_tmp[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                    p_loc_tmp[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                    p_loc_tmp[1][3] = -p_loc_tmp[0][4];                   // (2, 4) -z
                    p_loc_tmp[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                    p_loc_tmp[2][3] = -p_loc_tmp[0][5];                   // (3, 4) y
                    p_loc_tmp[2][4] = -p_loc_tmp[1][5];                   // (3, 5) -x
#endif
#if 1
                    p_loc_tmp[0][3] = -(node_coarse_z - node_aggregation_z); // (1, 4) z, phi modified
                    p_loc_tmp[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y, phi modified
                    p_loc_tmp[1][4] = p_loc_tmp[0][3];                       // (2, 5) z, phi modified
                    p_loc_tmp[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x, phi modified
                    p_loc_tmp[2][3] = -p_loc_tmp[1][5];                      // (3, 4) -x, phi modified
                    p_loc_tmp[2][4] = p_loc_tmp[0][5];                       // (3, 5) -y, phi modified
#endif

                    // p_loc in prolongation operator location
                    /*
                     * (m_tmp, n_tmp) <-> {(m_tmp + 1) * 6 - 1, (n_tmp + 1) * 6 - 1}
                     */
                    for (int index_tmp_i = 0; index_tmp_i < 6; ++index_tmp_i)
                    {
                        for (int index_tmp_j = 0; index_tmp_j < 6; ++index_tmp_j)
                        {
                            PetscCall(MatSetValue((mla_ctx->mla + cnt_level)->prolongation,
                                                  m_prolongation_block_tmp * 6 + index_tmp_i,
                                                  n_prolongation_block_tmp * 6 + index_tmp_j,
                                                  p_loc_tmp[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }
                else if (order_rbm == 2)
                {
                    for (int index_tmp = 0; index_tmp < 9; ++index_tmp)
                    {
                        p_loc_tmp[index_tmp][index_tmp] = 1.;
                    }
#if 0
                    p_loc_tmp[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                    p_loc_tmp[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                    p_loc_tmp[1][3] = -p_loc_tmp[0][4];                   // (2, 4) -z
                    p_loc_tmp[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                    p_loc_tmp[2][3] = -p_loc_tmp[0][5];                   // (3, 4) y
                    p_loc_tmp[2][4] = -p_loc_tmp[1][5];                   // (3, 5) -x
#endif
#if 0
                    p_loc_tmp[0][4] = -(node_coarse_z - node_aggregation_z); // (1, 5) z
                    p_loc_tmp[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y
                    p_loc_tmp[1][3] = -p_loc_tmp[0][4];                      // (2, 4) -z
                    p_loc_tmp[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x
                    p_loc_tmp[2][3] = -p_loc_tmp[0][5];                      // (3, 4) y
                    p_loc_tmp[2][4] = -p_loc_tmp[1][5];                      // (3, 5) -x
#endif
#if 1
                    p_loc_tmp[0][3] = -(node_coarse_z - node_aggregation_z); // (1, 4) z, phi modified
                    p_loc_tmp[0][5] = -(node_aggregation_y - node_coarse_y); // (1, 6) -y, phi modified
                    p_loc_tmp[1][4] = p_loc_tmp[0][3];                       // (2, 5) z, phi modified
                    p_loc_tmp[1][5] = -(node_coarse_x - node_aggregation_x); // (2, 6) x, phi modified
                    p_loc_tmp[2][3] = -p_loc_tmp[1][5];                      // (3, 4) -x, phi modified
                    p_loc_tmp[2][4] = p_loc_tmp[0][5];                       // (3, 5) -y, phi modified
#endif

#if 0
                    p_loc_tmp[0][6] = (node_aggregation_z - node_coarse_z) *
                                      (node_coarse_x - node_aggregation_x); // (1, 7) -zx
#endif
#if 0
                    p_loc_tmp[0][6] = (node_coarse_z - node_aggregation_z) *
                                      (node_coarse_x - node_aggregation_x); // (1, 7) zx, modified
#endif
#if 1
                    p_loc_tmp[0][6] = (node_coarse_z - node_aggregation_z) *
                                      (node_coarse_x - node_aggregation_x); // (1, 7) zx, phi modified
#endif

#if 0 
                    p_loc_tmp[0][8] = (node_aggregation_z - node_coarse_z) *
                                      (node_coarse_y - node_aggregation_y); // (1, 9) -zy
#endif
#if 0
                    p_loc_tmp[0][8] = -(node_aggregation_z - node_coarse_z) *
                                      (node_coarse_y - node_aggregation_y) / 2.; // (1, 9) zy/2, modified
#endif
#if 1
                    p_loc_tmp[0][8] = -(node_aggregation_z - node_coarse_z) *
                                      (node_coarse_y - node_aggregation_y) / 2.; // (1, 9) zy/2, phi modified
#endif

#if 0
                    p_loc_tmp[1][7] = p_loc_tmp[0][8]; // (2, 8) -zy
                    p_loc_tmp[1][8] = p_loc_tmp[0][6]; // (2, 9) -zx
#endif
#if 0
                    p_loc_tmp[1][7] = 2. * p_loc_tmp[0][8]; // (2, 8) zy, modified
                    p_loc_tmp[1][8] = p_loc_tmp[0][6] / 2.; // (2, 9) zx/2, modified
#endif
#if 1
                    p_loc_tmp[1][7] = 2. * p_loc_tmp[0][8]; // (2, 8) zy, phi modified
                    p_loc_tmp[1][8] = p_loc_tmp[0][6] / 2.; // (2, 9) zx/2, phi modified
#endif

#if 0
                    p_loc_tmp[2][6] = (node_coarse_x - node_aggregation_x) *
                                      (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) xx/2
#endif
                    p_loc_tmp[2][6] = -(node_coarse_x - node_aggregation_x) *
                                      (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) -xx/2, modified
                    p_loc_tmp[2][7] = -(node_coarse_y - node_aggregation_y) *
                                      (node_coarse_y - node_aggregation_y) / 2.; // (3, 8) -yy/2, modified
                    p_loc_tmp[2][8] = -(node_coarse_x - node_aggregation_x) *
                                      (node_coarse_y - node_aggregation_y) / 2.; // (3, 9) -xy/2, modified

#if 1 // lame constant modify
                    p_loc_tmp[1][7] *= lame_constant_modify;
                    p_loc_tmp[2][7] *= lame_constant_modify;
#endif

#if 0
                    p_loc_tmp[3][7] = p_loc_tmp[0][5];       // (4, 8) -y, modified
                    p_loc_tmp[3][8] = -p_loc_tmp[1][5] / 2.; // (4, 9) -x/2, modified
                    p_loc_tmp[4][6] = -2. * p_loc_tmp[3][8]; // (5, 7) x, modified
                    p_loc_tmp[4][8] = -p_loc_tmp[3][7] / 2.; // (5, 9) y/2, modified
#endif
#if 0
                    p_loc_tmp[3][6] = p_loc_tmp[1][5];      // (4, 7) x, phi modified
                    p_loc_tmp[3][8] = -p_loc_tmp[0][5] / 2.; // (4, 9) y/2, phi modified
                    p_loc_tmp[4][7] = 2. * p_loc_tmp[3][8]; // (5, 8) y, phi modified
                    p_loc_tmp[4][8] = p_loc_tmp[3][6] / 2.; // (5, 9) x/2, phi modified
#endif

                    // p_loc in prolongation operator location
                    /*
                     * (m_tmp, n_tmp) <-> {(m_tmp + 1) * 9 - 1, (n_tmp + 1) * 9 - 1}
                     */
                    for (int index_tmp_i = 0; index_tmp_i < 9; ++index_tmp_i)
                    {
                        for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                        {
                            PetscCall(MatSetValue((mla_ctx->mla + cnt_level)->prolongation,
                                                  m_prolongation_block_tmp * 9 + index_tmp_i,
                                                  n_prolongation_block_tmp * 9 + index_tmp_j,
                                                  p_loc_tmp[index_tmp_i][index_tmp_j],
                                                  INSERT_VALUES));
                        }
                    }
                }

                data_node_aggregation = data_node_aggregation->next;
            }
        }
        PetscCall(MatAssemblyBegin((mla_ctx->mla + cnt_level)->prolongation, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd((mla_ctx->mla + cnt_level)->prolongation, MAT_FINAL_ASSEMBLY));
#if 0
        printf("\n\n==== level %d, neighbouring prolongation operator:\n", cnt_level);
        PetscCall(MatView((mla_ctx->mla + cnt_level)->prolongation, PETSC_VIEWER_STDOUT_WORLD));
#endif // neighbouring prolongation operator

#endif // test metis aggregation implementation

#if 0
        PetscCall(MatDuplicate(mla_ctx->metis_mla[cnt_level].prolongation,
                               MAT_COPY_VALUES,
                               &((mla_ctx->mla + cnt_level)->prolongation)));
#endif // test metis aggregation implementation, added

        PetscCall(MatPtAP((mla_ctx->mla + cnt_level)->operator_fine,
                          (mla_ctx->mla + cnt_level)->prolongation,
                          MAT_INITIAL_MATRIX,
                          PETSC_DETERMINE,
                          &((mla_ctx->mla + cnt_level)->operator_coarse)));
#if 0
        printf("\n\n==== level %d, neighbouring coarse operator:\n", cnt_level);
        PetscCall(MatView((mla_ctx->mla + cnt_level)->operator_coarse, PETSC_VIEWER_STDOUT_WORLD));
#endif

        // ksp and pc object
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_presmooth)));
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_postsmooth)));
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->ksp_coarse)));
        PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_presmooth)));
        PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_postsmooth)));
        PetscCall(PCCreate(PETSC_COMM_WORLD, &((mla_ctx->mla + cnt_level)->pc_coarse)));

        ++cnt_level;
        printf("in while loop, cnt_level = %d\n", cnt_level);

// memory free
#if 1
        if (order_rbm == 1)
        {
            for (int index = 0; index < 6; ++index)
            {
                free(p_loc_tmp[index]);
            }
            free(p_loc_tmp);
        }
        else if (order_rbm == 2)
        {
            for (int index = 0; index < 9; ++index)
            {
                free(p_loc_tmp[index]);
            }
            free(p_loc_tmp);
        }
#endif // test metis aggregation implementation
    }

    mla_ctx->num_level = (cnt_level < num_level) ? cnt_level : num_level;

// free memory
#if 1
    for (int index = 0; index < 6; ++index)
    {
        free(p_loc[index]);
    }
    free(p_loc);
#endif // test metis aggregation implementation

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

int MLASolver(const MeshGraph *graph,
              MySolver *mysolver,
              const ConfigJSON *config,
              int order_rbm,
              MLAContext *mla_ctx,
              int mla_phase)
{
    // int mla_phase = config->mla_config.mla_phase;
    int num_level = config->mla_config.mla_level;
    PetscLogDouble time1, time2;

    // setup phase
    if (mla_phase == 0 || mla_phase == 2)
    {
        /*
         * numer of levels
         * finest mesh
         * rbm order (prolongation operator constructor)
         */
        PetscCall(PetscTime(&time1));
        MLASolverSetupPhase(mysolver, graph, num_level, order_rbm, mla_ctx);

        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> setup time: %g (s)\n", time2 - time1));
    }

    // solve phase
    if (mla_phase == 1 || mla_phase == 2)
    {
        PetscCall(PetscTime(&time1));

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

        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> solve time: %g (s)\n", time2 - time1));
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
