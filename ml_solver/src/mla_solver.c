#include "../include/main.h"

void MLASolverSolvePhase(const ConfigJSON *config,
                         MLAContext *mla_ctx,
                         int order_rbm,
                         MySolver *mysolver)
{
}

void MLASolverSetupPhase(MySolver *mysolver,
                         const MeshGraph *graph,
                         int num_level,
                         int order_rbm,
                         MLAContext *mla_ctx)
{
    if (mla_ctx->setup == 1)
    {
        // setup has beed done
        return;
    }

    mla_ctx->setup = 1;

    mla_ctx->mla = (MLAGraph *)malloc(num_level * sizeof(MLAGraph));
    assert(mla_ctx->mla);

    /*
     * 1. mla_ctx->num_level = num_level
     * 2. in coarset mesh, number of points no less than 50
     */
    int cnt_level = 0; // level count

    // finest and 2nd finest mesh
    MeshGraph *graph_tmp = CreateMeshGraph(graph->size);
    CopyMeshGraph(graph_tmp, graph);

    (mla_ctx->mla + cnt_level)->fine = CreateMeshGraph(graph->size);
    CopyMeshGraph((mla_ctx->mla + cnt_level)->fine, graph);
#if 0
    puts("\n>>>> 1st fine mesh:");
    PrintMeshGraph((mla_ctx->mla + cnt_level)->fine);
#endif

    (mla_ctx->mla + cnt_level)->aggregation = AggregationMeshGraph(graph_tmp);
#if 0
    puts("\n>>>> 1st aggregation mesh:");
    PrintMeshGraph((mla_ctx->mla + cnt_level)->aggregation);
#endif

    (mla_ctx->mla + cnt_level)->coarse = CreateMeshGraph((mla_ctx->mla + cnt_level)->aggregation->size);

    (mla_ctx->mla + cnt_level)->coarse_node = (MeshNode *)malloc((mla_ctx->mla + cnt_level)->aggregation->size * sizeof(MeshNode));
    assert((mla_ctx->mla + cnt_level)->coarse_node);

    AssembleCoarseMeshNode((mla_ctx->mla + cnt_level)->aggregation,
                           (mla_ctx->mla + cnt_level)->coarse_node);

    AssembleCoarseMeshGraph((mla_ctx->mla + cnt_level)->coarse,
                            (mla_ctx->mla + cnt_level)->aggregation,
                            (mla_ctx->mla + cnt_level)->fine,
                            (mla_ctx->mla + cnt_level)->coarse_node);
    (mla_ctx->mla + cnt_level)->level = cnt_level;
#if 0
    puts("\n>>>> 1st coarse mesh:");
    PrintMeshGraph((mla_ctx->mla + cnt_level)->coarse);
#endif

    // free temporary memory
    // free(data_coarse_node);
    ClearMeshGraph(graph_tmp);

    // prolongation and coarse operator between finest mesh and 2nd finest mesh
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
    printf(">>>>>>>> setup phase:\n %d level: prolongation_block_row = %d, col = %d\n",
           cnt_level,
           prolongation_block_row,
           prolongation_block_col);
    printf("%d level, size of prolongation: m = %d, n = %d\n", cnt_level,
           m_prolongation,
           n_prolongation);
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
                p_loc[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                p_loc[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                p_loc[1][3] = -p_loc[0][4];                       // (2, 4) -z
                p_loc[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                p_loc[2][3] = -p_loc[0][5];                       // (3, 4) y
                p_loc[2][4] = -p_loc[1][5];                       // (3, 5) -x

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
                p_loc[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                p_loc[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                p_loc[1][3] = -p_loc[0][4];                       // (2, 4) -z
                p_loc[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                p_loc[2][3] = -p_loc[0][5];                       // (3, 4) y
                p_loc[2][4] = -p_loc[1][5];                       // (3, 5) -x
                p_loc[0][6] = (node_aggregation_z - node_coarse_z) *
                              (node_coarse_x - node_aggregation_x); // (1, 7) -zx
                p_loc[0][8] = (node_aggregation_z - node_coarse_z) *
                              (node_coarse_y - node_aggregation_y); // (1, 9) -zy
                p_loc[1][7] = p_loc[0][8];                          // (2, 8) -zy
                p_loc[1][8] = p_loc[0][6];                          // (2, 9) -zx
                p_loc[2][6] = (node_coarse_x - node_aggregation_x) *
                              (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) xx/2
                p_loc[2][7] = (node_coarse_y - node_aggregation_y) *
                              (node_coarse_y - node_aggregation_y) / 2.; // (3, 8) yy/2
                p_loc[2][8] = (node_coarse_x - node_aggregation_x) *
                              (node_coarse_y - node_aggregation_y); // (3, 9) xy
#if 0
                p_loc[3][7] = p_loc[0][4];                          // (4, 8) z
                p_loc[3][8] = p_loc[0][5];                          // (4, 9) -y
                p_loc[4][6] = p_loc[1][3];                          // (5, 7) -z
                p_loc[4][8] = p_loc[1][5];                          // (5, 9) x
                p_loc[5][6] = p_loc[2][3];                          // (6, 7) y
                p_loc[5][7] = p_loc[2][4];                          // (6, 8) -x
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

    // coarse operator
    PetscCall(MatPtAP(mysolver->solver_a, (mla_ctx->mla + cnt_level)->prolongation,
                      MAT_INITIAL_MATRIX, PETSC_DETERMINE,
                      &((mla_ctx->mla + cnt_level)->operator_coarse)));

#if 0
    puts(">>>> 1st level, coarse mesh:");
    PrintMeshGraph((mla_ctx->mla + cnt_level)->coarse);
#endif

// next level, after 1st level
/*
 * the coarset mesh contains at least 50 nodes
 */
// while (cnt_level < num_level && (mla_ctx->mla + cnt_level)->coarse->size >= 50)
// while ((mla_ctx->mla + cnt_level)->coarse->size >= 50)
#if 0
    printf(">>>> num_level = %d, cnt_level = %d, coarse num_vertex = %d\n", num_level,
           cnt_level,
           (mla_ctx->mla + cnt_level)->coarse->size);
#endif // 1st level
    while (cnt_level < num_level && (mla_ctx->mla + cnt_level)->coarse->size >= 50)
    {
        ++cnt_level;

        (mla_ctx->mla + cnt_level)->fine = CreateMeshGraph((mla_ctx->mla + cnt_level - 1)->coarse->size);
        CopyMeshGraph((mla_ctx->mla + cnt_level)->fine, (mla_ctx->mla + cnt_level - 1)->coarse);

        graph_tmp = CreateMeshGraph((mla_ctx->mla + cnt_level - 1)->coarse->size);
        CopyMeshGraph(graph_tmp, (mla_ctx->mla + cnt_level - 1)->coarse);

        (mla_ctx->mla + cnt_level)->aggregation = AggregationMeshGraph(graph_tmp);

        (mla_ctx->mla + cnt_level)->coarse = CreateMeshGraph((mla_ctx->mla + cnt_level)->aggregation->size);
        (mla_ctx->mla + cnt_level)->coarse_node = (MeshNode *)malloc((mla_ctx->mla + cnt_level)->aggregation->size * sizeof(MeshNode));
        assert((mla_ctx->mla + cnt_level)->coarse_node);

        AssembleCoarseMeshNode((mla_ctx->mla + cnt_level)->aggregation,
                               (mla_ctx->mla + cnt_level)->coarse_node);

        AssembleCoarseMeshGraph((mla_ctx->mla + cnt_level)->coarse,
                                (mla_ctx->mla + cnt_level)->aggregation,
                                (mla_ctx->mla + cnt_level)->fine,
                                (mla_ctx->mla + cnt_level)->coarse_node);
        (mla_ctx->mla + cnt_level)->level = cnt_level;
#if 0
        printf(">>>> num_level = %d, cnt_level = %d, coarse num_vertex = %d\n", num_level,
               cnt_level,
               (mla_ctx->mla + cnt_level)->coarse->size);
#endif // cnt_level-th level

        // prolongation and coarse operator
        prolongation_block_row = (mla_ctx->mla + cnt_level)->fine->size;   // vertex in fine mesh
        prolongation_block_col = (mla_ctx->mla + cnt_level)->coarse->size; // vertex in coarse mesh

        if (order_rbm == 1)
        {
            // 6 dof in coarse point
            m_prolongation = prolongation_block_row * 6; // row of prolongation operator
            n_prolongation = prolongation_block_col * 6; // col of prolongation operator
        }
        else if (order_rbm == 2)
        {
            // 9 dof in coarse point
            m_prolongation = prolongation_block_row * 9; // row of prolongation operator
            n_prolongation = prolongation_block_col * 9; // col of prolongation operator
        }

#if 1
        printf(">>>>>>>> setup phase:\n %d level: prolongation_block_row = %d, col = %d\n",
               cnt_level,
               prolongation_block_row,
               prolongation_block_col);
        printf("%d level, size of prolongation: m = %d, n = %d\n",
               cnt_level,
               m_prolongation,
               n_prolongation);
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
                        p_loc_tmp[index_tmp][index_tmp] = 1.;
                    }
                    p_loc_tmp[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                    p_loc_tmp[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                    p_loc_tmp[1][3] = -p_loc_tmp[0][4];                   // (2, 4) -z
                    p_loc_tmp[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                    p_loc_tmp[2][3] = -p_loc_tmp[0][5];                   // (3, 4) y
                    p_loc_tmp[2][4] = -p_loc_tmp[1][5];                   // (3, 5) -x

                    // p_loc_tmp in prolongation operator location
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
                    p_loc_tmp[0][4] = node_coarse_z - node_aggregation_z; // (1, 5) z
                    p_loc_tmp[0][5] = node_aggregation_y - node_coarse_y; // (1, 6) -y
                    p_loc_tmp[1][3] = -p_loc_tmp[0][4];                   // (2, 4) -z
                    p_loc_tmp[1][5] = node_coarse_x - node_aggregation_x; // (2, 6) x
                    p_loc_tmp[2][3] = -p_loc_tmp[0][5];                   // (3, 4) y
                    p_loc_tmp[2][4] = -p_loc_tmp[1][5];                   // (3, 5) -x
                    p_loc_tmp[0][6] = (node_aggregation_z - node_coarse_z) *
                                      (node_coarse_x - node_aggregation_x); // (1, 7) -zx
                    p_loc_tmp[0][8] = (node_aggregation_z - node_coarse_z) *
                                      (node_coarse_y - node_aggregation_y); // (1, 9) -zy
                    p_loc_tmp[1][7] = p_loc_tmp[0][8];                      // (2, 8) -zy
                    p_loc_tmp[1][8] = p_loc_tmp[0][6];                      // (2, 9) -zx
                    p_loc_tmp[2][6] = (node_coarse_x - node_aggregation_x) *
                                      (node_coarse_x - node_aggregation_x) / 2.; // (3, 7) xx/2
                    p_loc_tmp[2][7] = (node_coarse_y - node_aggregation_y) *
                                      (node_coarse_y - node_aggregation_y) / 2.; // (3, 8) yy/2
                    p_loc_tmp[2][8] = (node_coarse_x - node_aggregation_x) *
                                      (node_coarse_y - node_aggregation_y); // (3, 9) xy
#if 0
                p_loc_tmp[3][7] = p_loc_tmp[0][4];                          // (4, 8) z
                p_loc_tmp[3][8] = p_loc_tmp[0][5];                          // (4, 9) -y
                p_loc_tmp[4][6] = p_loc_tmp[1][3];                          // (5, 7) -z
                p_loc_tmp[4][8] = p_loc_tmp[1][5];                          // (5, 9) x
                p_loc_tmp[5][6] = p_loc_tmp[2][3];                          // (6, 7) y
                p_loc_tmp[5][7] = p_loc_tmp[2][4];                          // (6, 8) -x
#endif

                    // p_loc_tmp in prolongation operator location
                    /*
                     * (m_tmp, n_tmp) <-> {(m_tmp + 1) * 9 - 1, (n_tmp + 1) * 9 - 1}
                     */
                    for (int index_tmp_i = 0; index_tmp_i < 9; ++index_tmp_i)
                    {
                        for (int index_tmp_j = 0; index_tmp_j < 9; ++index_tmp_j)
                        {
                            PetscCall(MatSetValue((mla_ctx->mla + cnt_level)->prolongation,
                                                  m_prolongation_block_tmp * 6 + index_tmp_i,
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

        // coarse operator
        PetscCall(MatPtAP((mla_ctx->mla + cnt_level - 1)->operator_coarse,
                          (mla_ctx->mla + cnt_level)->prolongation,
                          MAT_INITIAL_MATRIX, PETSC_DETERMINE,
                          &((mla_ctx->mla + cnt_level)->operator_coarse)));

        // free temporary memory
        ClearMeshGraph(graph_tmp);
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
    }
#if 1
    printf(">>>> num_level = %d, cnt_level = %d, coarse num_vertex = %d\n", num_level,
           cnt_level,
           (mla_ctx->mla + cnt_level)->coarse->size);
#endif // final leve

    mla_ctx->num_level = (cnt_level < num_level) ? cnt_level : num_level; // number of level

// free p_loc
#if 1
    for (int index = 0; index < 6; ++index)
    {
        free(p_loc[index]);
    }
    free(p_loc);
#endif
}

void MLASolverRelativeResidual(MySolver *mysolver, double *value)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    *value = r_norm_2 / b_norm_2;
}

void MLASolver(const MeshGraph *graph,
               MySolver *mysolver,
               const ConfigJSON *config,
               int order_rbm,
               MLAContext *mla_ctx)
{
    int mla_phase = config->mla_config.mla_phase;
    int num_level = config->mla_config.mla_level;

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
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));

            /*
             * solve phase recursive implementation
             *     linear system
             *     mla context, contains setup information
             *     config, pre- and post- smooth times
             *     a special case, rbm order is 2 and level is 1, coarse operator need shift
             */
            MLASolverSolvePhase(config, mla_ctx, order_rbm, mysolver);

            MLASolverRelativeResidual(mysolver, &rela_resid);
            ++cnt;
        }
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));
    }

#if 0
    puts("\n\n==== multilevel information ====");
    for (int index = 0; index <= mla_ctx->num_level; ++index)
    {
        printf(">>>> in neighbouring level %d:\n", (mla_ctx->mla + index)->level);
        printf("there are %d vertices in fine mesh, %d vertices in coarse mesh!\n\n",
               (mla_ctx->mla + index)->fine->size,
               (mla_ctx->mla + index)->coarse->size);
    }
#endif // basic information of multilevel
}
