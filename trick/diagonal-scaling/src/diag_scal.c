#include "main.h"

int LinsysDiagScal(Mat *mat, Vec *rhs, Vec *diag)
{
    PetscInt local_m, global_m;
    PetscCall(MatGetSize(*mat, &global_m, NULL));
    PetscCall(MatGetLocalSize(*mat, &local_m, NULL));
    PetscInt rstart;
    PetscCall(MatGetOwnershipRange(*mat, &rstart, NULL));

    Vec diag_tmp;
    PetscCall(VecCreate(PetscObjectComm((PetscObject)*mat), &diag_tmp));
    PetscCall(VecSetSizes(diag_tmp, local_m, global_m));
    PetscCall(VecSetFromOptions(diag_tmp));
    PetscCall(MatGetDiagonal(*mat, diag_tmp));

    PetscScalar *diag_array;
    PetscCall(VecGetArray(diag_tmp, &diag_array));

    for (PetscInt i = 0; i < local_m; ++i)
    {
        PetscScalar dval = diag_array[i];
        if (dval > 0)
        {
            dval = PetscSqrtReal(dval);
        }
        else if (dval == 0.0)
        {
            dval = 1.0;
        }
        else
        {
            dval = PetscSqrtReal(-dval);
        }

        PetscInt gid = rstart + i;
        PetscCall(VecSetValues(*diag, 1, &gid, &dval, INSERT_VALUES));
    }
    PetscCall(VecRestoreArray(diag_tmp, &diag_array));
    PetscCall(VecDestroy(&diag_tmp));

    PetscCall(VecAssemblyBegin(*diag));
    PetscCall(VecAssemblyEnd(*diag));

    Vec diag_inv;
    PetscCall(VecDuplicate(*diag, &diag_inv));
    PetscCall(VecCopy(*diag, diag_inv));
    PetscCall(VecReciprocal(diag_inv));
    PetscCall(VecPointwiseMult(*rhs, *rhs, diag_inv));

    PetscCall(MatDiagonalScale(*mat, diag_inv, diag_inv));

    PetscCall(VecDestroy(&diag_inv));
    return 0;
}

int SolDiagScal(Vec *tilde_x /*tilde solution*/, Vec *diag /*processed diagonal*/, Vec *x /*solution*/)
{
    const PetscScalar *tilde_array, *diag_array;
    PetscScalar *x_array;

    PetscInt local_size;
    PetscCall(VecGetLocalSize(*x, &local_size));

    PetscCall(VecGetArrayRead(*tilde_x, &tilde_array));
    PetscCall(VecGetArrayRead(*diag, &diag_array));
    PetscCall(VecGetArray(*x, &x_array));

    for (PetscInt i = 0; i < local_size; ++i)
    {
        x_array[i] = tilde_array[i] / diag_array[i];
    }

    PetscCall(VecRestoreArrayRead(*tilde_x, &tilde_array));
    PetscCall(VecRestoreArrayRead(*diag, &diag_array));
    PetscCall(VecRestoreArray(*x, &x_array));

    return 0;
}
