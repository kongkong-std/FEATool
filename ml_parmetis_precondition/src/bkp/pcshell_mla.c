#include "../include/main.h"

PetscErrorCode MLAShellPCSetup(PC pc)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    int mla_phase = 0;
    MLASolver(mla_ctx->graph,
              &(mla_ctx->mysolver),
              &(mla_ctx->config),
              mla_ctx->order_rbm,
              mla_ctx,
              mla_phase);

    // puts("==== hhh, setup has been done ====");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MLAShellPCApply(PC pc, Vec x_in, Vec x_out)
{
    PetscFunctionBeginUser;

    MLAContext *mla_ctx;
    PetscCall(PCShellGetContext(pc, &mla_ctx));

    PetscCall(VecCopy(x_in, mla_ctx->mysolver.solver_b));

    int mla_phase = 1;
    MLASolver(mla_ctx->graph,
              &(mla_ctx->mysolver),
              &(mla_ctx->config),
              mla_ctx->order_rbm,
              mla_ctx,
              mla_phase);

    PetscCall(VecCopy(mla_ctx->mysolver.solver_x, x_out));

    PetscFunctionReturn(PETSC_SUCCESS);
}
