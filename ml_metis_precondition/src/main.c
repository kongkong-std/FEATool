#include "../include/main.h"

int DeepCopyMLAContextMySolver(MySolver *dst, MySolver *src);

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
    PetscInt order_rbm = 0, angle_type = 0;

    PetscCall(PetscOptionsGetString(NULL, NULL, "-config", path_config, sizeof(path_config), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Config file: %s\n", path_config));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-order_rbm", &order_rbm, &path_flag));
    if (path_flag)
    {
        // parameter order_rbm
        /*
         * order_rbm = 1, classical rbm type prolongation operator
         * order_rbm = 2, krichhoff type prolongation operator
         */
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "order_rbm = %" PetscInt_FMT "\n", order_rbm));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-angle_type", &angle_type, &path_flag));
    if (path_flag)
    {
        // parameter angle_type
        /*
         * angle_type = 0, rotational angle with axis
         * angle_type = 1, normal rotational angle with axis
         */
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "angle_type = %" PetscInt_FMT "\n", angle_type));
    }

    // user-defined precondition
    PetscBool def_pc_mla = PETSC_FALSE;
    PetscBool def_pc_metis_mla = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_mla", &def_pc_mla, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_metis_mla", &def_pc_metis_mla, NULL));

    // mla_ctx
    MLAContext mla_ctx;
    mla_ctx.setup = 0;
    mla_ctx.path_config = path_config;
    mla_ctx.order_rbm = order_rbm;
    mla_ctx.angle_type = angle_type;

    if (ConfigParse(mla_ctx.path_config, &(mla_ctx.config)) != 0)
    {
        exit(EXIT_FAILURE);
    }
    mla_ctx.metis_mla = (MetisMLAGraph *)malloc(mla_ctx.config.mla_config.mla_level *
                                                sizeof(MetisMLAGraph));
    assert(mla_ctx.metis_mla);
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
    PetscCall(DeepCopyMLAContextMySolver(&(mla_ctx.mysolver), &mysolver));

    /*
     * step 2, gmsh data process
     */
    mla_ctx.data_gmsh = (DataGmsh *)malloc(sizeof(DataGmsh));
    assert(mla_ctx.data_gmsh);
    MetisFileProcessGmsh(mla_ctx.config.file_config.file_mesh,
                         mla_ctx.data_gmsh);
#if 0
    puts("\n==== gmsh data information ====");
    printf("number of nodes: %d\n", mla_ctx.data_gmsh->nn);
    printf("number of elements: %d\n", mla_ctx.data_gmsh->ne);
    for (int index = 0; index < mla_ctx.data_gmsh->nn; ++index)
    {
        printf("node %d: %021.16le\t%021.16le\t%021.16le\n", index,
               mla_ctx.data_gmsh->coordinates[3 * index],
               mla_ctx.data_gmsh->coordinates[3 * index + 1],
               mla_ctx.data_gmsh->coordinates[3 * index + 2]);
    }
    printf("number of elements in boundary: %d\n", mla_ctx.data_gmsh->ne_bd);
    for (int index = 0; index < mla_ctx.data_gmsh->ne_bd + 1; ++index)
    {
        printf("data_gmsh.eptr_bd[%d] = %" PRIDX "\n", index, mla_ctx.data_gmsh->eptr_bd[index]);
    }
    printf("number of elements in inner: %d\n", mla_ctx.data_gmsh->ne_in);
    for (int index = 0; index < mla_ctx.data_gmsh->ne_in + 1; ++index)
    {
        printf("data_gmsh.eptr_in[%d] = %" PRIDX "\n", index, mla_ctx.data_gmsh->eptr_in[index]);
    }
    printf("nodes in inner elements:\n");
    for (int index = 0; index < mla_ctx.data_gmsh->ne_in; ++index)
    {
        printf("element %d: ", index);
        for (int index_i = 0; index_i < mla_ctx.data_gmsh->nne_in; ++index_i)
        {
            printf("%" PRIDX "\t", mla_ctx.data_gmsh->eind_in[index * mla_ctx.data_gmsh->nne_in + index_i]);
        }
        putchar('\n');
    }
#endif // print gmsh file data

#if 0
    puts("\n==== test metis function ====");
    TestMetisFunctionGmsh(*(mla_ctx.data_gmsh));
#endif // test metis function

    /*
     * step 3, calling MetisMLASolver()
     */
    // PetscCall(MetisMLASolver(&mla_ctx, 2));

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

#if 0
    MLASolver(mla_ctx.graph,
              &(mla_ctx.mysolver),
              &(mla_ctx.config),
              mla_ctx.order_rbm,
              &mla_ctx,
              2);
#endif // adjacency list implementation

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

    if (def_pc_metis_mla)
    {
        PetscCall(PCSetType(mysolver.pc, PCSHELL));
        PetscCall(PCShellSetContext(mysolver.pc, &mla_ctx));
        PetscCall(PCShellSetSetUp(mysolver.pc, MetisMLAShellPCSetup));
        PetscCall(PCShellSetApply(mysolver.pc, MetisMLAShellPCApply));
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
#if 0
    free(mla_ctx.data_gmsh->coordinates);
    free(mla_ctx.data_gmsh->eptr_bd);
    free(mla_ctx.data_gmsh->eind_bd);
    free(mla_ctx.data_gmsh->eptr_in);
    free(mla_ctx.data_gmsh->eind_in);
    free(mla_ctx.data_gmsh->epart_in);
    free(mla_ctx.data_gmsh->npart_in);
    free(mla_ctx.data_gmsh);
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
#endif // free memory

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}

int DeepCopyMLAContextMySolver(MySolver *dst, MySolver *src)
{
    PetscCall(VecDuplicate(src->solver_b, &(dst->solver_b)));
    PetscCall(VecCopy(src->solver_b, dst->solver_b));
    PetscCall(VecDuplicate(src->solver_x, &(dst->solver_x)));
    PetscCall(VecCopy(src->solver_x, dst->solver_x));
    PetscCall(VecDuplicate(src->solver_r, &(dst->solver_r)));
    PetscCall(VecCopy(src->solver_r, dst->solver_r));

    PetscCall(MatDuplicate(src->solver_a, MAT_COPY_VALUES, &(dst->solver_a)));

    return 0;
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
 *     -pc_type hypre
 *     -pc_hypre_type boomeramg
 *     -pc_hypre_boomeramg_numfunctions
 */
