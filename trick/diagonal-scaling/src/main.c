#include "main.h"

int SolverPetscResidualCheck(Mat solver_a, Vec solver_x, Vec solver_b, Vec solver_r);

int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    char *path_mat = NULL, *path_rhs = NULL;
    int flag_diag_scal = 0;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-file_mat", argv[index]))
        {
            path_mat = argv[index + 1];
        }
        if (strstr("-file_rhs", argv[index]))
        {
            path_rhs = argv[index + 1];
        }
        if (strstr("-flag_diag_scal", argv[index]))
        {
            flag_diag_scal = atoi(argv[index + 1]);
        }
    }

    Mat solver_a;
    Vec solver_b, solver_x, solver_r;

    // binary matrix load
    PetscViewer fd;

    // matrix file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_mat, FILE_MODE_READ, &fd));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &solver_a));
    PetscCall(MatLoad(solver_a, fd));
    PetscCall(PetscViewerDestroy(&fd));

    // rhs file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_rhs, FILE_MODE_READ, &fd));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &solver_b));
    PetscCall(VecLoad(solver_b, fd));
    PetscCall(PetscViewerDestroy(&fd));

    // sol vector and residual vector
    PetscCall(VecDuplicate(solver_b, &solver_x));
    PetscCall(VecDuplicate(solver_b, &solver_r));

    // size of matrix
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==== basic information of linear system ====\n"));
    PetscInt m_mat = 0, n_mat = 0; // nnz_mat = 0;
    PetscCall(MatGetSize(solver_a, &m_mat, &n_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix Row = %" PetscInt_FMT ", matrix Column = %" PetscInt_FMT "\n", m_mat, n_mat));
    MatInfo info_mat;
    PetscCall(MatGetInfo(solver_a, MAT_GLOBAL_SUM, &info_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix nz_allocated = %ld, matrix nz_used = %ld, matrix nz_unneeded = %ld\n",
                          (long)(info_mat.nz_allocated), (long)(info_mat.nz_used), (long)(info_mat.nz_unneeded)));
    PetscInt n_vec = 0;
    PetscCall(VecGetSize(solver_b, &n_vec));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "vector Row = %" PetscInt_FMT "\n", n_vec));

    PetscReal b_norm_2 = 0.;
    PetscCall(VecNorm(solver_b, NORM_2, &b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "initial     || b ||_2 = %021.16le\n", b_norm_2));

    // SolverPetscResidualCheck(solver_a, solver_b, solver_r);

    PetscLogDouble time1, time2;

    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));

    if (flag_diag_scal == 0)
    {
        PetscCall(KSPSetOperators(ksp, solver_a, solver_a));
        PetscCall(KSPSetFromOptions(ksp));

        PetscCall(PetscTime(&time1));
        PetscCall(KSPSolve(ksp, solver_b, solver_x));
        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> ksp solve time: %g (s)\n", time2 - time1));

        SolverPetscResidualCheck(solver_a, solver_x, solver_b, solver_r);
    }

    Mat tilde_solver_a;
    Vec tilde_solver_b, tilde_solver_x, tilde_solver_r, tilde_diag_solver_a;

    if (flag_diag_scal == 1)
    {

        // diagonal scaling
        PetscCall(VecDuplicate(solver_b, &tilde_solver_b));
        PetscCall(VecCopy(solver_b, tilde_solver_b));
        PetscCall(VecDuplicate(tilde_solver_b, &tilde_solver_x));
        PetscCall(VecDuplicate(tilde_solver_b, &tilde_solver_r));
        PetscCall(VecDuplicate(tilde_solver_b, &tilde_diag_solver_a));

        PetscCall(MatDuplicate(solver_a, MAT_COPY_VALUES, &tilde_solver_a));

        LinsysDiagScal(&tilde_solver_a, &tilde_solver_b, &tilde_diag_solver_a);

#if 0
        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            (void)MPI_Barrier(PETSC_COMM_WORLD);
            if (index_p == my_rank)
            {
                PetscInt global_row_start, global_row_end;
                PetscCall(MatGetOwnershipRange(tilde_solver_a, &global_row_start, &global_row_end));
                for (int index = global_row_start; index < global_row_end; ++index)
                {
                    PetscScalar val_tmp;
                    PetscCall(VecGetValues(tilde_diag_solver_a, 1, &index, &val_tmp));
                    printf("%d\t%021.16le\n", index, val_tmp);
                }
            }
        }
#endif

        PetscCall(KSPSetOperators(ksp, tilde_solver_a, tilde_solver_a));
        PetscCall(KSPSetFromOptions(ksp));

        PetscCall(PetscTime(&time1));
        PetscCall(KSPSolve(ksp, tilde_solver_b, tilde_solver_x));
        PetscCall(PetscTime(&time2));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> ksp solve time: %g (s)\n", time2 - time1));

        SolverPetscResidualCheck(tilde_solver_a, tilde_solver_x, tilde_solver_b, tilde_solver_r);

        SolDiagScal(&tilde_solver_x, &tilde_diag_solver_a, &solver_x);

        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "======== test original solution ========\n"));

        SolverPetscResidualCheck(solver_a, solver_x, solver_b, solver_r);
    }

    // free memory
    PetscCall(MatDestroy(&solver_a));
    PetscCall(VecDestroy(&solver_b));
    PetscCall(VecDestroy(&solver_x));
    PetscCall(VecDestroy(&solver_r));

    if (flag_diag_scal == 1)
    {
        PetscCall(MatDestroy(&tilde_solver_a));
        PetscCall(VecDestroy(&tilde_solver_b));
        PetscCall(VecDestroy(&tilde_solver_x));
        PetscCall(VecDestroy(&tilde_solver_r));
    }

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}

int SolverPetscResidualCheck(Mat solver_a, Vec solver_x, Vec solver_b, Vec solver_r)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(solver_a, solver_x, solver_r));
    PetscCall(VecAXPY(solver_r, -1., solver_b));
    PetscCall(VecNorm(solver_r, NORM_2, &r_norm_2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "            || b ||_2 = %021.16le\n", b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "|| r ||_2 / || b ||_2 = %021.16le\n", r_norm_2 / b_norm_2));

    return 0;
}

// command line
/*
 * mpirun -np <n> ./app_petsc_exe
 *     -file_mat </path/to/petscbin/mat>
 *     -file_rhs </path/to/petscbin/rhs>
 *     -flag_diag_scal <0/1>
 *     -ksp_type cg
 *     -ksp_max_it 1000
 *     -ksp_rtol 1e-8
 *     -ksp_norm_type unpreconditioned
 *     -ksp_monitor_true_residual
 */
