#include "../include/main.h"

#define MAT_COL_MAJOR(r, c, m) ((c) * (m) + (r))

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

int SAMGSACoarseOperator(SAMGCtx **samg_ctx /*samg context data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    SAMGCtx *data_samg_ctx = *samg_ctx;

    int ps_num_steps = data_samg_ctx->data_cfg.cfg_mg.ps_num_steps;
    int ps_type = data_samg_ctx->data_cfg.cfg_mg.ps_type;
    double ps_scale = data_samg_ctx->data_cfg.cfg_mg.ps_scale;

    int num_level = data_samg_ctx->num_level;
    int cnt_level = 0;

    int ncol_p = 0;

    for (cnt_level = 0; cnt_level < num_level; ++cnt_level)
    {
        // setting for prolongaiton smoother op_s
        data_samg_ctx->levels[cnt_level].op_s.num_steps = ps_num_steps;
        data_samg_ctx->levels[cnt_level].op_s.smoother_type = ps_type;
        data_samg_ctx->levels[cnt_level].op_s.smoother_scale = ps_scale;

        // smoothing op_ua_p to op_sa_p
        PetscCall(SAMGSmoothedProlongation(data_samg_ctx->levels + cnt_level));

        // PtAP
        PetscCall(MatPtAP(data_samg_ctx->levels[cnt_level].op_f,
                          data_samg_ctx->levels[cnt_level].op_sa_p,
                          MAT_INITIAL_MATRIX,
                          PETSC_DETERMINE,
                          &(data_samg_ctx->levels[cnt_level].op_c)));

        // next level
        PetscCall(MatDuplicate(data_samg_ctx->levels[cnt_level].op_c,
                               MAT_COPY_VALUES,
                               &(data_samg_ctx->levels[cnt_level + 1].op_f)));
    }

    return 0;
}

int SAMGUACoarseOperator(SAMGCtx **samg_ctx /*samg context data*/)
{
    // int my_rank, nprocs;
    // MPI_Comm comm;
    // comm = PETSC_COMM_WORLD;
    // MPI_Comm_rank(comm, &my_rank);
    // MPI_Comm_size(comm, &nprocs);

    SAMGCtx *data_samg_ctx = *samg_ctx;

    int num_level = data_samg_ctx->num_level;
    int cnt_level = 0;

    for (cnt_level = 0; cnt_level < num_level; ++cnt_level)
    {
        // PtAP
        PetscCall(MatPtAP(data_samg_ctx->levels[cnt_level].op_f,
                          data_samg_ctx->levels[cnt_level].op_ua_p,
                          MAT_INITIAL_MATRIX,
                          PETSC_DETERMINE,
                          &(data_samg_ctx->levels[cnt_level].op_c)));

        // next level
        PetscCall(MatDuplicate(data_samg_ctx->levels[cnt_level].op_c,
                               MAT_COPY_VALUES,
                               &(data_samg_ctx->levels[cnt_level + 1].op_f)));
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

    PetscCall(SAMGLevelQMatrix(&samg_ctx)); // Q matrix, size and distribution same with near null space of fine-level

    PetscCall(SAMGTentativeProlongationOperator(&samg_ctx)); // tentative prolongation operator constructor

    if (sa_flag == 1)
    {
        // SA
        PetscCall(PetscPrintf(comm, "==== Smoothed Aggregation-based Multigrid\n"));
        PetscCall(SAMGSACoarseOperator(&samg_ctx));
    }
    else if (sa_flag == 0)
    {
        // UA
        PetscCall(PetscPrintf(comm, "==== Unsmoothed Aggregation-based Multigrid\n"));
        PetscCall(SAMGUACoarseOperator(&samg_ctx));
    }

    // while (cnt_level < cfg_mg_num_level &&
    //        samg_ctx->levels[cnt_level].data_f_mesh.np > cfg_mg_num_coarse_vtx)
    // {
    //     ++cnt_level;
    // }

    return 0;
}
