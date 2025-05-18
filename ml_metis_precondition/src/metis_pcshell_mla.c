#include "../include/main.h"

PetscErrorCode MetisMLAShellPCSetup(PC pc)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    int mla_phase = 0;
    MetisMLASolver(mla_ctx, mla_phase);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MetisMLAShellPCApply(PC pc, Vec x_in, Vec x_out)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    PetscCall(VecCopy(x_in, mla_ctx->mysolver.solver_b));

    int mla_phase = 1;
    MetisMLASolver(mla_ctx, mla_phase);

    PetscCall(VecCopy(mla_ctx->mysolver.solver_x, x_out));

    PetscFunctionReturn(PETSC_SUCCESS);
}
