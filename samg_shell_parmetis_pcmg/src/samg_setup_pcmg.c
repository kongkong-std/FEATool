#include "../include/main.h"

int PCMGSetupFromSAMG(int sa_flag /*flag of sa*/,
                      SAMGCtx *samg_ctx /*samg context data*/,
                      MySolver *mysolver /*solver data*/)
{
    PetscCall(PCSetType(mysolver->pc, PCMG));
    PetscCall(PCMGSetLevels(mysolver->pc, samg_ctx->num_level + 1, NULL));
    PetscCall(PCMGSetType(mysolver->pc, PC_MG_MULTIPLICATIVE));
    PetscCall(PCMGSetCycleType(mysolver->pc, PC_MG_CYCLE_V));

    // prolongation oprator
    if (sa_flag == 1)
    {
        // SA
        for (int level = 1; level < samg_ctx->num_level + 1; ++level)
        {
            PetscCall(PCMGSetInterpolation(mysolver->pc, level, samg_ctx->levels[samg_ctx->num_level - level].op_sa_p));
        }
    }
    else if (sa_flag == 0)
    {
        // UA
        for (int level = 1; level < samg_ctx->num_level + 1; ++level)
        {
            PetscCall(PCMGSetInterpolation(mysolver->pc, level, samg_ctx->levels[samg_ctx->num_level - level].op_ua_p));
        }
    }

    // level operator
    for (int level = 1; level < samg_ctx->num_level + 1; ++level)
    {
        PetscCall(PCMGSetOperators(mysolver->pc, level,
                                   samg_ctx->levels[samg_ctx->num_level - level].op_f,
                                   samg_ctx->levels[samg_ctx->num_level - level].op_f));
    }

    PetscCall(PCMGSetOperators(mysolver->pc, 0,
                               samg_ctx->levels[samg_ctx->num_level - 1].op_c,
                               samg_ctx->levels[samg_ctx->num_level - 1].op_c));

    PetscCall(PCMGSetNumberSmooth(mysolver->pc, samg_ctx->data_cfg.cfg_mg.pre_smooth));

    // smoother
    for (int level = 0; level < samg_ctx->num_level + 1; ++level)
    {
        KSP smoother;
        PC pc_smoother;
        PetscCall(PCMGGetSmoother(mysolver->pc, level, &smoother));
        PetscCall(KSPGetPC(smoother, &pc_smoother));
        // PetscCall(PCSetType(pc_smoother, PCLU));
        // PetscCall(PCFactorSetMatSolverType(pc_smoother, MATSOLVERMUMPS));
        // PetscCall(PCHYPRESetType(pc_smoother, "pilut"));
        PetscCall(KSPView(smoother, PETSC_VIEWER_STDOUT_WORLD));
    }

    // coarsest level solver
    KSP coarse_ksp;
    PC coarse_pc;

    PetscCall(PCMGGetCoarseSolve(mysolver->pc, &coarse_ksp));
    PetscCall(KSPGetPC(coarse_ksp, &coarse_pc));
    PetscCall(PCSetType(coarse_pc, PCLU)); // Direct solver for coarse grid
    // PetscCall(PCFactorSetMatSolverType(coarse_pc, MATSOLVERMUMPS));
    PetscCall(KSPSetType(coarse_ksp, KSPPREONLY));
    PetscCall(KSPView(coarse_ksp, PETSC_VIEWER_STDOUT_WORLD));
    // MatSolverType coarse_stype;
    // PetscCall(PCFactorGetMatSolverType(coarse_pc, &coarse_stype));
    // PetscPrintf(PETSC_COMM_WORLD,"======== Coarse solver type: %s\n", coarse_stype);

    return 0;
}
