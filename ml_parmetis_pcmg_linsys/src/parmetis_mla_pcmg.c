#include "../include/main.h"

int PCMGSetupFromMLA(const MLAContext *mla_ctx /*mla context data*/,
                     MySolver *mysolver /*solver data*/)
{
    PetscCall(PCSetType(mysolver->pc, PCMG));
    // PetscCall(PCMGSetLevels(mysolver->pc, mla_ctx->num_level + 1, &comm));
    PetscCall(PCMGSetLevels(mysolver->pc, mla_ctx->num_level + 1, NULL));
    PetscCall(PCMGSetType(mysolver->pc, PC_MG_MULTIPLICATIVE));
    PetscCall(PCMGSetCycleType(mysolver->pc, PC_MG_CYCLE_V));

    // prolongation oprator
    for (int level = 1; level < mla_ctx->num_level + 1; ++level)
    {
#if 0
        PetscCall(PCMGSetInterpolation(mysolver->pc, level, mla_ctx->metis_mla[mla_ctx->num_level - level].usa_prolongation));
#endif // unsmoothed aggregation prolongation operator

#if 1
        PetscCall(PCMGSetInterpolation(mysolver->pc, level, mla_ctx->metis_mla[mla_ctx->num_level - level].sa_prolongation));
#endif // smoothed aggregation prolongation operator
    }

    // level operator
    for (int level = 1; level < mla_ctx->num_level + 1; ++level)
    {
#if 1
        double shift = 1e-12;
        PetscCall(MatShift(mla_ctx->metis_mla[mla_ctx->num_level - level].operator_fine, shift));
#endif // diagonal shift

        PetscCall(PCMGSetOperators(mysolver->pc, level,
                                   mla_ctx->metis_mla[mla_ctx->num_level - level].operator_fine,
                                   mla_ctx->metis_mla[mla_ctx->num_level - level].operator_fine));
    }
#if 1
    double shift = 1e-12;
    PetscCall(MatShift(mla_ctx->metis_mla[mla_ctx->num_level - 1].operator_coarse, shift));
#endif // diagonal shift
    PetscCall(PCMGSetOperators(mysolver->pc, 0,
                               mla_ctx->metis_mla[mla_ctx->num_level - 1].operator_coarse,
                               mla_ctx->metis_mla[mla_ctx->num_level - 1].operator_coarse)); // coarest level

    PetscCall(PCMGSetNumberSmooth(mysolver->pc, mla_ctx->config.mla_config.pre_smooth_v));

    // smoother
    for (int level = 0; level < mla_ctx->num_level + 1; ++level)
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
