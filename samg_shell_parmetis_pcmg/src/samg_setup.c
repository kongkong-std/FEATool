#include "../include/main.h"

int SAMGLevel0Mesh(const CfgJson *data_cfg /*config data*/,
                   MeshData *data_f_mesh /*fine-level mesh data*/,
                   MeshData *data_c_mesh /*coarse-level mesh data*/)
{
    FileProcessMeshVtx(data_cfg->cfg_file.file_vtx,
                       &data_f_mesh->data_vtx);
    FileProcessMeshAdj(data_cfg->cfg_file.file_adj,
                       &data_f_mesh->data_adj);

    return 0;
}

int SAMGLevelMesh(int cfg_mg_num_level /*config number of levels*/,
                  SAMGCtx **samg_ctx /*samg context data*/)
{
    SAMGCtx *data_samg_ctx = *samg_ctx;

    int cfg_mg_num_level = data_samg_ctx->data_cfg.cfg_mg.num_level;

    data_samg_ctx->levels = (MGLevel *)malloc(cfg_mg_num_level * sizeof(MGLevel));
    assert(data_samg_ctx->levels);

    int cnt_level = 0;
    PetscCall(SAMGLevel0Mesh(&data_samg_ctx->data_cfg,
                             &data_samg_ctx->levels[cnt_level].data_f_mesh,
                             &data_samg_ctx->levels[cnt_level].data_c_mesh));

    while (cnt_level < cfg_mg_num_level)
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
    // int cfg_mg_pre_smooth = samg_ctx->data_cfg.cfg_mg.pre_smooth;
    // int cfg_mg_post_smooth = samg_ctx->data_cfg.cfg_mg.post_smooth;
    int cfg_mg_num_coarse_vtx = samg_ctx->data_cfg.cfg_mg.num_coarse_vtx;
    // int cfg_mg_est_size_agg = samg_ctx->data_cfg.cfg_mg.est_size_agg;
    int cfg_mg_ps_num_steps = samg_ctx->data_cfg.cfg_mg.ps_num_steps;
    int cfg_mg_ps_type = samg_ctx->data_cfg.cfg_mg.ps_type;
    double cfg_mg_ps_scale = samg_ctx->data_cfg.cfg_mg.ps_scale;

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
