#include "../include/main.h"

void DeepCopyMLAContextMySolver(MySolver *dst, MySolver *src);

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

    // user-defined precondition
    PetscBool def_pc_mla = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_mla", &def_pc_mla, NULL));

    // mla_ctx
    MLAContext mla_ctx;
    mla_ctx.setup = 0;
    mla_ctx.path_config = path_config;
    mla_ctx.order_rbm = order_rbm;

    if (ConfigParse(mla_ctx.path_config, &(mla_ctx.config)) != 0)
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
    SolverPetscInitialize(mla_ctx.config.file_config.file_mat,
                          mla_ctx.config.file_config.file_rhs,
                          &mysolver);

    // deep copy
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mla_ctx.mysolver.ksp)));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &(mla_ctx.mysolver.pc)));
    DeepCopyMLAContextMySolver(&(mla_ctx.mysolver), &mysolver);

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

    MeshFileProcess(mla_ctx.config.file_config.file_mesh,
                    mla_ctx.config.mesh_label_config.label_bound,
                    mla_ctx.config.mesh_label_config.label_omega,
                    &data_list_phy_tag,
                    &data_list_node,
                    &data_list_ele_bound,
                    &data_list_ele_omega);

    mla_ctx.graph = CreateMeshGraph(data_list_node.size);
    AssembleMeshGraph(mla_ctx.graph, &data_list_node, &data_list_ele_omega);
#if 0
    puts("\n>>>> initial mesh:");
    PrintMeshGraph(graph);
#endif

#if 1
    // pc shell
    PetscCall(KSPSetOperators(mysolver.ksp, mysolver.solver_a, mysolver.solver_a));
    PetscCall(KSPGetPC(mysolver.ksp, &(mysolver.pc)));
    if (def_pc_mla)
    {
        PetscCall(PCSetType(mysolver.pc, PCSHELL));
        PetscCall(PCShellSetContext(mysolver.pc, &mla_ctx));
        PetscCall(PCShellSetSetUp(mysolver.pc, MLAShellPCSetup));
        PetscCall(PCShellSetApply(mysolver.pc, MLAShellPCApply));
    }

    PetscCall(KSPSetFromOptions(mysolver.ksp));
    PetscCall(KSPSolve(mysolver.ksp, mysolver.solver_b, mysolver.solver_x));
#endif // pcshell

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

#if 0
    mla_ctx.mla = (MLAGraph *)malloc(config.mla_config.mla_level * sizeof(MLAGraph));
    assert(mla_ctx.mla);
#endif

#if 0
    // mla phase
    /*
     * 0 represents setup only
     * 1 represents solve only
     * 2 represents setup+solve
     */
    int mla_phase = 2;

    MLASolver(mla_ctx.graph,
              &(mla_ctx.mysolver),
              &(mla_ctx.config),
              mla_ctx.order_rbm,
              &mla_ctx,
              mla_phase);

    DeepCopyMLAContextMySolver(&mysolver, &(mla_ctx.mysolver));
#endif // mla solver

    // check, computing residual
    SolverPetscResidualCheck(&mysolver);

    // free memory
    for (int index = 0; index < mla_ctx.num_level; ++index)
    {
        free(mla_ctx.mla[index].coarse_node);
        PetscCall(MatDestroy(&(mla_ctx.mla[index].prolongation)));
        PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_coarse)));
        PetscCall(MatDestroy(&(mla_ctx.mla[index].operator_fine)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_coarse)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_presmooth)));
        PetscCall(KSPDestroy(&(mla_ctx.mla[index].ksp_postsmooth)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_coarse)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_presmooth)));
        PetscCall(PCDestroy(&(mla_ctx.mla[index].pc_postsmooth)));
        ClearMeshGraph(mla_ctx.mla[index].aggregation);
        ClearMeshGraph(mla_ctx.mla[index].fine);
        ClearMeshGraph(mla_ctx.mla[index].coarse);
        ClearMeshGraph(mla_ctx.mla[index].mesh_tmp);
    }
    free(mla_ctx.mla);
    ClearList(&data_list_ele_omega);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_node);
    ClearList(&data_list_phy_tag);
    ClearMeshGraph(mla_ctx.graph);
    PetscCall(MatDestroy(&mysolver.solver_a));
    PetscCall(VecDestroy(&mysolver.solver_b));
    PetscCall(VecDestroy(&mysolver.solver_x));
    PetscCall(VecDestroy(&mysolver.solver_r));
    PetscCall(KSPDestroy(&mysolver.ksp));
    // PetscCall(PCDestroy(&mysolver.pc));

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}

void DeepCopyMLAContextMySolver(MySolver *dst, MySolver *src)
{
    PetscCall(VecDuplicate(src->solver_b, &(dst->solver_b)));
    PetscCall(VecCopy(src->solver_b, dst->solver_b));
    PetscCall(VecDuplicate(src->solver_x, &(dst->solver_x)));
    PetscCall(VecCopy(src->solver_x, dst->solver_x));
    PetscCall(VecDuplicate(src->solver_r, &(dst->solver_r)));
    PetscCall(VecCopy(src->solver_r, dst->solver_r));

    PetscCall(MatDuplicate(src->solver_a, MAT_COPY_VALUES, &(dst->solver_a)));
}

// command line
/* >>>> for mla solver case
 * ./app_petsc_exe -config </path/to/config/file>
 *     -order_rbm <int, 1st order, 2nd order>
 */

/*
 * >>>> for mla preconditioning case
 * ./app_petsc_exe -config </path/to/config/file>
 *     -order_rbm <int, 1st order, 2nd order>
 *     -ksp_type fgmres
 *     -def_pc_mla
 *     -ksp_max_it 2000
 *     -ksp_rtol 1e-8
 *     -ksp_gmres_restart 50
 *     -ksp_monitor_true_residual
 *     -ksp_norm_type unpreconditioned
 */
