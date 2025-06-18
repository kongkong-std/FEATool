#include "../include/main.h"

PetscErrorCode ParMetisMLAShellPCSetup(PC pc)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    int mla_phase = 0;
    ParMetisMLASolver(mla_ctx, mla_phase);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ParMetisMLAShellPCApply(PC pc, Vec x_in, Vec x_out)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    PetscCall(VecCopy(x_in, mla_ctx->mysolver.solver_b));

    int mla_phase = 1;
    ParMetisMLASolver(mla_ctx, mla_phase);

    //PetscCall(VecCopy(mla_ctx->mysolver.solver_x, x_out));
    PetscCall(VecCopy(x_in, x_out));

    PetscFunctionReturn(PETSC_SUCCESS);
}
