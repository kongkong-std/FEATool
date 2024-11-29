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

    ConfigJSON config;
    if (ConfigParse(path_config, &config) != 0)
    {
        exit(EXIT_FAILURE);
    }
#if 0
    printf("File Config:\n");
    printf("file_mat: %s\n", config.file_config.file_mat);
    printf("file_rhs: %s\n", config.file_config.file_rhs);
    printf("file_mesh: %s\n", config.file_config.file_mesh);

    printf("MLA Config:\n");
    printf("pre_smooth_v: %d\n", config.mla_config.pre_smooth_v);
    printf("post_smooth_v: %d\n", config.mla_config.post_smooth_v);
    printf("mla_max_it: %d\n", config.mla_config.mla_max_it);
    printf("mla_rtol: %021.16le\n", config.mla_config.mla_rtol);
    printf("mla_level: %d\n", config.mla_config.mla_level);
    printf("mla_phase: %d\n", config.mla_config.mla_phase);
#endif // print config json

    /*
     * step 1, linear system assembling
     */
    MySolver mysolver;
    SolverPetscInitialize(config.file_config.file_mat, config.file_config.file_rhs, &mysolver);

    /*
     * step 2, mesh process
     */
    // mesh file process
    GenericList data_list_phy_tag, data_list_node;
    GenericList data_list_ele_bound, data_list_ele_omega;
    InitializeList(&data_list_phy_tag);
    InitializeList(&data_list_node);
    InitializeList(&data_list_ele_bound);
    InitializeList(&data_list_ele_omega);

    MeshFileProcess(config.file_config.file_mesh,
                    config.mesh_label_config.label_bound,
                    config.mesh_label_config.label_omega,
                    &data_list_phy_tag,
                    &data_list_node,
                    &data_list_ele_bound,
                    &data_list_ele_omega);

    MeshGraph *graph = CreateMeshGraph(data_list_node.size);
    AssembleMeshGraph(graph, &data_list_node, &data_list_ele_omega);
#if 0
    puts("\n>>>> initial mesh:");
    PrintMeshGraph(graph);
#endif

    /*
     * step 3, mla solver
     *     2 main parts:
     *         1. setup phase
     *         2. solve phase
     *     at least contains parameters list:
     *         maximum iterations
     *         tolerance
     *         number of levels
     *         pre-smooth times
     *         post-smooth times
     *         phase (setup/solve)
     *         rbm order (prolongation operator constructor)
     *         linear system
     *         graph data
     *         mla context
     */
    MLAContext mla_ctx;
    mla_ctx.setup = 0;
    MLASolver(graph, &mysolver, &config, order_rbm, &mla_ctx);

    // check, computing residual
    SolverPetscResidualCheck(&mysolver);

    // free memory
    for (int index = 0; index < mla_ctx.num_level; ++index)
    {
        free(mla_ctx.mla[index].coarse_node);
    }
    free(mla_ctx.mla);
    ClearMeshGraph(graph);
    ClearList(&data_list_ele_omega);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_node);
    ClearList(&data_list_phy_tag);
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

// command line
/*
 * ./app_petsc_exe -config </path/to/config/file>
 *     -order_rbm <int, 1st order, 2nd order>
 */
