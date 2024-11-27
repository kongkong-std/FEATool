#include "../include/main.h"

int main(int argc, char **argv)
{
    int myrank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscBool path_flag;

    char path_config[PETSC_MAX_PATH_LEN];
    int order_rbm = 0;

    PetscCall(PetscOptionsGetString(NULL, NULL, "-config", path_config, sizeof(path_config), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Config file: %s\n", path_config));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-order_rbm", &order_rbm, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "order_rbm = %d\n", order_rbm));
    }

    ConfigFile file_config;
    ConfigMLA mla_config;
    if (ConfigParse(path_config, &file_config, &mla_config) != 0)
    {
        exit(EXIT_FAILURE);
    }
#if 0
    printf("File Config:\n");
    printf("file_mat: %s\n", file_config.file_mat);
    printf("file_rhs: %s\n", file_config.file_rhs);
    printf("file_mesh: %s\n", file_config.file_mesh);

    printf("MLA Config:\n");
    printf("pre_smooth_v: %d\n", mla_config.pre_smooth_v);
    printf("post_smooth_v: %d\n", mla_config.post_smooth_v);
    printf("mla_max_it: %d\n", mla_config.mla_max_it);
    printf("mla_rtol: %021.16le\n", mla_config.mla_rtol);
    printf("mla_level: %d\n", mla_config.mla_level);
    printf("mla_phase: %d\n", mla_config.mla_phase);
#endif // print config json

    /*
     * step 1, linear system assembling
     */
    MySolver mysolver;
    SolverPetscInitialize(file_config.file_mat, file_config.file_rhs, &mysolver);

    // computing residual
    SolverPetscResidualCheck(&mysolver);

    // free memory
    PetscCall(MatDestroy(&mysolver.solver_a));
    PetscCall(VecDestroy(&mysolver.solver_b));
    PetscCall(VecDestroy(&mysolver.solver_x));
    PetscCall(VecDestroy(&mysolver.solver_r));
    PetscCall(KSPDestroy(&mysolver.ksp));
    PetscCall(PCDestroy(&mysolver.pc));

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}