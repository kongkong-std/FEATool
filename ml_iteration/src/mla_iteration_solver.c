#include "../include/main.h"

void MLAPreSmoothPhase(MySolver *mysolver, int v_pre_smooth)
{
#if 0
    PetscViewer viewer;
    PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer));
    PetscCall(KSPMonitorSet(mysolver->ksp, KSPMonitorTrueResidual, &viewer, NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n==== pre-smoothing:\n"));
#endif

    // ==== mla step 1, gauss-seidel pre-smoothing
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));
    // PetscCall(KSPSetType(mysolver->ksp, KSPPREONLY));
    PetscCall(KSPSetType(mysolver->ksp, KSPRICHARDSON));
    PetscCall(KSPGetPC(mysolver->ksp, &(mysolver->pc)));
    PetscCall(PCSetType(mysolver->pc, PCSOR));
    PetscCall(PCSORSetOmega(mysolver->pc, 1.)); // gauss-seidel
    PetscCall(PCSORSetIterations(mysolver->pc, v_pre_smooth, v_pre_smooth));
    PetscCall(PCSORSetSymmetric(mysolver->pc, SOR_SYMMETRIC_SWEEP));
    PetscCall(KSPSetNormType(mysolver->ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetFromOptions(mysolver->ksp));

    PetscCall(KSPSetTolerances(mysolver->ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

    PetscCall(KSPSetInitialGuessNonzero(mysolver->ksp, PETSC_TRUE));
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
}

void MLAPostSmoothPhase(MySolver *mysolver, int v_post_smooth)
{
#if 0
    PetscViewer viewer;
    PetscCall(PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer));
    PetscCall(KSPMonitorSet(mysolver->ksp, KSPMonitorTrueResidual, &viewer, NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n==== post-smoothing:\n"));
#endif

    // ==== mla step , gauss-seidel post-smoothing
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));
    PetscCall(KSPSetType(mysolver->ksp, KSPRICHARDSON));
    PetscCall(KSPGetPC(mysolver->ksp, &(mysolver->pc)));
    PetscCall(PCSetType(mysolver->pc, PCSOR));
    PetscCall(PCSORSetOmega(mysolver->pc, 1.)); // gauss-seidel
    PetscCall(PCSORSetIterations(mysolver->pc, v_post_smooth, v_post_smooth));
    PetscCall(PCSORSetSymmetric(mysolver->pc, SOR_SYMMETRIC_SWEEP));
    PetscCall(KSPSetNormType(mysolver->ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetFromOptions(mysolver->ksp));

    PetscCall(KSPSetTolerances(mysolver->ksp, 1e-10, 1e-10, PETSC_DEFAULT, 1));

    PetscCall(KSPSetInitialGuessNonzero(mysolver->ksp, PETSC_TRUE));
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
}

void MLASolvePhase(MySolver *mysolver, MLAGraph *mla, int gcr_restart)
{
    // computing residual
    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAYPX(mysolver->solver_r, -1., mysolver->solver_b)); // b - Ax, ATTENTION!!!

    int m_solver_a = 0, n_solver_a = 0;
    int m_prolongation, n_prolongation = 0;
    PetscCall(MatGetSize(mysolver->solver_a, &m_solver_a, &n_solver_a));
    PetscCall(MatGetSize(mla->prolongation, &m_prolongation, &n_prolongation));

#if 0
    printf("\n>>>>>>>> m_solver_a = %d, n_solver_a = %d\n", m_solver_a, n_solver_a);
    printf(">>>>>>>> m_prolongation = %d, n_prolongation = %d\n", m_prolongation, n_prolongation);

    int n_solver_r = 0;
    PetscCall(VecGetSize(mysolver->solver_r, &n_solver_r));
    printf(">>>>>>>> n_solver_r = %d\n", n_solver_r);
#endif

    // coarse linear system and correction
    /*
     *     A_H size: n_prolongation x n_prolongation
     *     r_H size: n_prolongation x 1
     *     e_H size: n_prolongation x 1
     */
    // Mat a_H;
    Vec r_H, e_H;
    KSP ksp_H;
    PC pc;

    // PetscCall(MatPtAP(mysolver->solver_a, mla->prolongation, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &a_H));
#if 1
    int m_a_H = 0, n_a_H = 0;
    // PetscCall(MatGetSize(a_H, &m_a_H, &n_a_H));
    PetscCall(MatGetSize(mla->operator_coarse, &m_a_H, &n_a_H));
    printf("\n>>>>>>>> m_a_H = %d, n_a_H = %d\n", m_a_H, n_a_H);
#endif // coarse operator

#if 0
    int m_a_H = 0, n_a_H = 0;
    PetscCall(MatGetSize(a_H, &m_a_H, &n_a_H));
    printf(">>>>>>>> size of a_H: m = %d, n = %d\n", m_a_H, n_a_H);
    double prolongation_f_norm = 0., coarse_f_norm = 0., fine_f_norm = 0., mat_tmp_f_norm = 0.;
    PetscCall(MatNorm(mysolver->solver_a, NORM_FROBENIUS, &fine_f_norm));
    PetscCall(MatNorm(mla->prolongation, NORM_FROBENIUS, &prolongation_f_norm));
    PetscCall(MatNorm(a_H, NORM_FROBENIUS, &coarse_f_norm));

#if 0
    Mat mat_tmp;
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat_tmp));
    PetscCall(MatSetSizes(mat_tmp, PETSC_DECIDE, PETSC_DECIDE, m_a_H, n_a_H));
    PetscCall(MatSetUp(mat_tmp));
    PetscCall(MatSetFromOptions(mat_tmp));
    PetscCall(MatDuplicate(a_H, MAT_DO_NOT_COPY_VALUES, &mat_tmp));
    PetscCall(MatAXPY(mat_tmp, -1., a_H, UNKNOWN_NONZERO_PATTERN));
    PetscCall(MatAXPY(mat_tmp, 1., mysolver->solver_a, UNKNOWN_NONZERO_PATTERN));
    PetscCall(MatNorm(mat_tmp, NORM_FROBENIUS, &mat_tmp_f_norm));
#endif

    printf(">>>> ||*||_f\n");
    printf("prolongation = %021.16le\n"
           "fine_mat = %021.16le\n"
           "coarse_mat = %021.16le\n"
           "minus_fine_coarse = %021.16le\n",
           prolongation_f_norm, fine_f_norm, coarse_f_norm, mat_tmp_f_norm);

    PetscCall(MatView(a_H, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatView(mla->prolongation, PETSC_VIEWER_STDOUT_WORLD));
#endif

    PetscCall(VecCreate(PETSC_COMM_WORLD, &e_H));
    PetscCall(VecSetSizes(e_H, PETSC_DECIDE, n_prolongation));
    PetscCall(VecSetFromOptions(e_H));
    PetscCall(VecDuplicate(e_H, &r_H));
    // PetscCall(VecDuplicate(mysolver->solver_r, &e_H));
    // PetscCall(VecView(e_H, PETSC_VIEWER_STDOUT_WORLD));

#if 0
    int n_r_H = 0;
    PetscCall(VecGetSize(e_H, &n_r_H));
    printf(">>>>>>>> n_r_H = %d\n", n_r_H);
#endif

    // r_H = Pt r
    PetscCall(MatMultTranspose(mla->prolongation, mysolver->solver_r, r_H));
    // PetscCall(VecSet(r_H, 1.));

    if (mla->coarse->size * 6 == n_prolongation)
    {
        // rbm order 1
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp_H));
        PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
        PetscCall(KSPSetOperators(ksp_H, mla->operator_coarse, mla->operator_coarse));
        PetscCall(KSPSetType(ksp_H, KSPGCR));
        // PetscCall(KSPGetPC(ksp_H, &pc));
        // PetscCall(PCSetType(pc, PCLU));
        PetscCall(KSPGCRSetRestart(ksp_H, gcr_restart));
        PetscCall(KSPSetNormType(ksp_H, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetFromOptions(ksp_H));
        PetscCall(KSPSetTolerances(ksp_H, 1e-10, 1e-10, PETSC_DEFAULT, 1000));

        PetscCall(KSPSolve(ksp_H, r_H, e_H));
    }
    else if (mla->coarse->size * 9 == n_prolongation)
    {
// rbm order 2
#if 0
        printf("==== rbm order 2\n");
        int null_space_coarse_size = mla->coarse->size * 3;
        Vec null_space_coarse_basis[null_space_coarse_size];

        int cnt_mod_0 = 0, cnt_mod_1 = 0, cnt_mod_2 = 0;

        for (int index = 0; index < null_space_coarse_size; ++index)
        {
            PetscCall(VecCreate(PETSC_COMM_WORLD, null_space_coarse_basis + index));
            PetscCall(VecSetSizes(null_space_coarse_basis[index], PETSC_DECIDE, n_prolongation));
            PetscCall(VecSetFromOptions(null_space_coarse_basis[index]));

            double *values_tmp = (double *)calloc(n_prolongation, sizeof(double));
            assert(values_tmp);

            if (index % 3 == 0)
            {
                values_tmp[6 + 9 * cnt_mod_0] = 1.;

                for (int index_i = 0; index_i < n_prolongation; ++index_i)
                {
                    PetscCall(VecSetValues(null_space_coarse_basis[index],
                                           1, &index_i, values_tmp + index_i, INSERT_VALUES));
                }
                PetscCall(VecAssemblyBegin(null_space_coarse_basis[index]));
                PetscCall(VecAssemblyEnd(null_space_coarse_basis[index]));

                ++cnt_mod_0;
            }
            else if (index % 3 == 1)
            {
                values_tmp[7 + 9 * cnt_mod_1] = 1.;

                for (int index_i = 0; index_i < n_prolongation; ++index_i)
                {
                    PetscCall(VecSetValues(null_space_coarse_basis[index],
                                           1, &index_i, values_tmp + index_i, INSERT_VALUES));
                }
                PetscCall(VecAssemblyBegin(null_space_coarse_basis[index]));
                PetscCall(VecAssemblyEnd(null_space_coarse_basis[index]));

                ++cnt_mod_1;
            }
            else if (index % 3 == 2)
            {
                values_tmp[8 + 9 * cnt_mod_2] = 1.;

                for (int index_i = 0; index_i < n_prolongation; ++index_i)
                {
                    PetscCall(VecSetValues(null_space_coarse_basis[index],
                                           1, &index_i, values_tmp + index_i, INSERT_VALUES));
                }
                PetscCall(VecAssemblyBegin(null_space_coarse_basis[index]));
                PetscCall(VecAssemblyEnd(null_space_coarse_basis[index]));

                ++cnt_mod_2;
            }

            // free memory
            free(values_tmp);

            double null_spcae_norm_tmp = 0.;
            PetscCall(MatMult(a_H, null_space_coarse_basis[index], e_H));
            PetscCall(VecNorm(e_H, NORM_2, &null_spcae_norm_tmp));
            printf(">>>> null space vector %d, Ax = %021.16le\n", index, null_spcae_norm_tmp);
        }

        MatNullSpace null_coarse_space;
        PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, null_space_coarse_size,
                                     null_space_coarse_basis, &null_coarse_space));
        PetscCall(MatSetNullSpace(a_H, null_coarse_space));
#endif // null space

        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp_H));
        PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
        PetscCall(KSPSetOperators(ksp_H, mla->operator_coarse, mla->operator_coarse));
        PetscCall(KSPSetType(ksp_H, KSPGCR));
        // PetscCall(KSPSetType(ksp_H, KSPGMRES));
        PetscCall(KSPGetPC(ksp_H, &pc));
        PetscCall(PCSetType(pc, PCSVD));
        PetscCall(KSPGCRSetRestart(ksp_H, gcr_restart));
        PetscCall(KSPSetNormType(ksp_H, KSP_NORM_UNPRECONDITIONED));
        PetscCall(KSPSetFromOptions(ksp_H));
        PetscCall(KSPSetTolerances(ksp_H, 1e-4, 1e-4, PETSC_DEFAULT, 1));

        PetscCall(KSPSolve(ksp_H, r_H, e_H));
    }

#if 1
    double norm_r_H = 0.;
    PetscCall(VecNorm(r_H, NORM_2, &norm_r_H));
    Vec tmp_r_H;
    PetscCall(VecDuplicate(r_H, &tmp_r_H));
    PetscCall(MatMult(mla->operator_coarse, e_H, tmp_r_H));
    PetscCall(VecAXPY(tmp_r_H, -1., r_H));
    double norm_tmp_r_H = 0.;
    PetscCall(VecNorm(tmp_r_H, NORM_2, &norm_tmp_r_H));

    double norm_e_H = 0.;
    PetscCall(VecNorm(e_H, NORM_2, &norm_e_H));

    printf(">>>>>>>> coarse correction: norm_r_H = %021.16le\n", norm_r_H);
    printf(">>>>>>>> coarse correction: norm_e_H = %021.16le\n", norm_e_H);
    printf(">>>>>>>> coarse correction: relative = %021.16le\n\n", norm_tmp_r_H / norm_r_H);
#endif

    // updating solution
    /*
     * x = x + P eH
     */
    PetscCall(MatMult(mla->prolongation, e_H, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_x, 1., mysolver->solver_r));

    // free memory
    PetscCall(KSPDestroy(&ksp_H));
    PetscCall(VecDestroy(&r_H));
    PetscCall(VecDestroy(&e_H));
}

void MLASetupPhase(MySolver *mysolver, MLAGraph *mla, int order_rbm)
{
    if (mla->prolongation_set == 1)
    {
        // prolongation operator has been constructed
        return;
    }

    int prolongation_block_row = mla->fine->size;    // vertex in fine mesh
    int prolongation_block_col = mla->coarse->size;  // vertex in coarse mesh
    int m_prolongation = prolongation_block_row * 6; // row of prolongation operator
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
    PetscCall(MatCreate(PETSC_COMM_WORLD, &(mla->prolongation)));
    PetscCall(MatSetSizes(mla->prolongation, PETSC_DECIDE, PETSC_DECIDE, m_prolongation, n_prolongation));
    PetscCall(MatSetUp(mla->prolongation));

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
    for (int index = 0; index < mla->aggregation->size; ++index)
    {
        MeshGraphAdjNode *data_node_coarse = mla->coarse->array[index].head;
        MeshGraphAdjNode *data_node_aggregation = mla->aggregation->array[index].head;

        // coarse node coordinate
        double node_coarse_x = data_node_coarse->node->node_x;
        double node_coarse_y = data_node_coarse->node->node_y;
        double node_coarse_z = data_node_coarse->node->node_z;

        // aggregation graph adjacency list
        for (int index_i = 0; index_i < mla->aggregation->array[index].size; ++index_i)
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
                        PetscCall(MatSetValue(mla->prolongation,
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
                        PetscCall(MatSetValue(mla->prolongation,
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
    PetscCall(MatAssemblyBegin(mla->prolongation, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(mla->prolongation, MAT_FINAL_ASSEMBLY));

    PetscCall(MatPtAP(mysolver->solver_a, mla->prolongation,
                      MAT_INITIAL_MATRIX, PETSC_DETERMINE,
                      &(mla->operator_coarse)));

    mla->prolongation_set = 1;
}

void MLARelativeResidual(MySolver *mysolver, double *value)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    *value = r_norm_2 / b_norm_2;
}

void MLAIterationSolver(MySolver *mysolver, MLAGraph *mla,
                        int mla_pahse, int gcr_restart,
                        double rtol, int max_it,
                        int v_pre_smooth, int v_post_smooth,
                        int order_rbm)
{
    // setup phase
    if (mla_pahse == 0 || mla_pahse == 2)
    {
        MLASetupPhase(mysolver, mla, order_rbm);
    }

    // solve phase
    if (mla_pahse == 1 || mla_pahse == 2)
    {

        int cnt = 0;
        double rela_resid = 0.;
        MLARelativeResidual(mysolver, &rela_resid);

        while (cnt < max_it && rela_resid > rtol)
        {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));
            MLAPreSmoothPhase(mysolver, v_pre_smooth);
            MLASolvePhase(mysolver, mla, gcr_restart);
            MLAPostSmoothPhase(mysolver, v_post_smooth);
            MLARelativeResidual(mysolver, &rela_resid);
            ++cnt;
        }
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%d MLA ||r(i)||/||b|| %021.16le\n", cnt, rela_resid));
    }
}
