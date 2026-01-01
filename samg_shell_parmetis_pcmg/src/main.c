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

    /*
     * SA flag: 1 represents using SA, 0 represents using UA
     *     SA: smoothed aggregation, UA: unsmoothed aggregation
     * default is 1
     */
    int sa_flag = 1;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-sa_flag", &sa_flag, &path_flag));
    if (path_flag)
    {
        // parameter sa_flag
        /*
         * sa_flag = 1, SA
         * sa_flag = 0, UA
         */
        PetscCall(PetscPrintf(comm, "sa_flag = %d\n", sa_flag));
    }

    SAMGCtx samg_ctx;

    int status_config_json = ParseConfig(comm, path_config, &(samg_ctx.data_cfg));
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
            printf("file_mat: %s\n", samg_ctx.data_cfg.cfg_file.file_mat);
            printf("file_rhs: %s\n", samg_ctx.data_cfg.cfg_file.file_rhs);
            printf("file_vtx: %s\n", samg_ctx.data_cfg.cfg_file.file_vtx);
            printf("file_adj: %s\n", samg_ctx.data_cfg.cfg_file.file_adj);
            printf("pre_smooth: %d\n", samg_ctx.data_cfg.cfg_mg.pre_smooth);
            printf("post_smooth: %d\n", samg_ctx.data_cfg.cfg_mg.post_smooth);
            printf("num_level: %d\n", samg_ctx.data_cfg.cfg_mg.num_level);
            printf("num_coarse_vtx: %d\n", samg_ctx.data_cfg.cfg_mg.num_coarse_vtx);
            printf("est_size_agg: %d\n", samg_ctx.data_cfg.cfg_mg.est_size_agg);
            puts("\n");
        }
    }
#endif // print config file

    // linear system file process
    MySolver mysolver;
    PetscCall(SolverInitializeWithFile(samg_ctx.data_cfg.cfg_file.file_mat,
                                       samg_ctx.data_cfg.cfg_file.file_rhs,
                                       &mysolver));
    PetscCall(PetscPrintf(comm, "==== Initial L2 norm of residual\n"));
    PetscCall(SolverPetscResidualCheck(&mysolver));

    // linear system copy to mg context
    PetscCall(DeepCopyMySolverData(&samg_ctx.mysolver, &mysolver));
#if 0
    PetscCall(PetscPrintf(comm, ">>>> SAMG context mysolver data:\n"));
    PetscCall(SolverPetscResidualCheck(&samg_ctx.mysolver));
#endif // check samg context mysolver data

    // setup + solve
    PetscCall(KSPSetFromOptions(mysolver.ksp));

    PetscLogDouble time1, time2;
    PetscCall(PetscTime(&time1));
    PetscCall(SAMGSetupPhase(&samg_ctx, sa_flag));
    PetscCall(PetscTime(&time2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> samg setup time: %g (s)\n", time2 - time1));

    PetscCall(KSPSetOperators(mysolver.ksp, mysolver.solver_a, mysolver.solver_a));
    PetscCall(KSPGetPC(mysolver.ksp, &(mysolver.pc)));
    PetscCall(PetscTime(&time1));
    PetscCall(KSPSolve(mysolver.ksp, mysolver.solver_b, mysolver.solver_x));
    PetscCall(PetscTime(&time2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> ksp solve time: %g (s)\n", time2 - time1));

    PetscCall(SolverPetscResidualCheck(&mysolver));

    // free memory
    for (int index = 0; index < samg_ctx.num_level; ++index)
    {
        free(samg_ctx.levels[index].data_f_mesh.data_vtx.idx);
        free(samg_ctx.levels[index].data_f_mesh.data_vtx.type);
        free(samg_ctx.levels[index].data_f_mesh.data_vtx.data_coor);
        free(samg_ctx.levels[index].data_f_mesh.data_adj.idx);
        free(samg_ctx.levels[index].data_f_mesh.data_adj.xadj);
        free(samg_ctx.levels[index].data_f_mesh.data_adj.adjncy);
        free(samg_ctx.levels[index].data_f_mesh.vtxdist);
        free(samg_ctx.levels[index].data_f_mesh.parts);
    }
    free(samg_ctx.levels);
    PetscCall(MatDestroy(&mysolver.solver_a));
    PetscCall(MatDestroy(&samg_ctx.mysolver.solver_a));
    PetscCall(VecDestroy(&mysolver.solver_b));
    PetscCall(VecDestroy(&mysolver.solver_x));
    PetscCall(VecDestroy(&mysolver.solver_r));
    PetscCall(VecDestroy(&samg_ctx.mysolver.solver_b));
    PetscCall(VecDestroy(&samg_ctx.mysolver.solver_x));
    PetscCall(VecDestroy(&samg_ctx.mysolver.solver_r));
    PetscCall(KSPDestroy(&mysolver.ksp));

    PetscCall(PetscFinalize());
    // MPI_Finalize();
    return 0;
}

// usage
/*
 * mpirun -np <nprocs> ./app_petsc_exe \
 *     -config </path/to/config/file>
 */
