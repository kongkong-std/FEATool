#include "../include/main.h"

int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    // MPI_Init(&argc, &argv);

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    // MPI_Comm comm = MPI_COMM_WORLD;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    PetscBool path_flag;

    char path_config[PETSC_MAX_PATH_LEN];

    PetscCall(PetscOptionsGetString(NULL, NULL, "-config", path_config, sizeof(path_config), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(comm, "Config file: %s\n", path_config));
    }

    SAMGCtx mla_ctx;

    int status_config_json = ParseConfig(comm, path_config, &(mla_ctx.data_cfg));
    if (status_config_json != 0)
    {
        exit(EXIT_FAILURE);
    }
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        if (index == my_rank)
        {
            printf("in rank %d =====\n", my_rank);
            printf("file_mat: %s\n", mla_ctx.data_cfg.cfg_file.file_mat);
            printf("file_rhs: %s\n", mla_ctx.data_cfg.cfg_file.file_rhs);
            printf("file_vtx: %s\n", mla_ctx.data_cfg.cfg_file.file_vtx);
            printf("file_adj: %s\n", mla_ctx.data_cfg.cfg_file.file_adj);
            printf("pre_smooth: %d\n", mla_ctx.data_cfg.cfg_mg.pre_smooth);
            printf("post_smooth: %d\n", mla_ctx.data_cfg.cfg_mg.post_smooth);
            printf("num_level: %d\n", mla_ctx.data_cfg.cfg_mg.num_level);
            printf("num_coarse_vtx: %d\n", mla_ctx.data_cfg.cfg_mg.num_coarse_vtx);
            printf("est_size_agg: %d\n", mla_ctx.data_cfg.cfg_mg.est_size_agg);
            puts("\n");
        }
    }
#endif // print config file

    // linear system file process
    MySolver mysolver;
    PetscCall(SolverInitializeWithFile(mla_ctx.data_cfg.cfg_file.file_mat,
                                       mla_ctx.data_cfg.cfg_file.file_rhs,
                                       &mysolver));
    PetscCall(PetscPrintf(comm, "==== Initial L2 norm of residual\n"));
    SolverPetscResidualCheck(&mysolver);

    PetscCall(PetscFinalize());
    // MPI_Finalize();
    return 0;
}

// usage
/*
 * mpirun -np <nprocs> ./app_petsc_exe \
 *     -config </path/to/config/file>
 */
