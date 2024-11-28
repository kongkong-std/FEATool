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

    // ==== step 1, linear system assembling ==== //
    /*
     * assembling linear system with petsc binary file
     */
    char path_mat[PETSC_MAX_PATH_LEN];
    char path_rhs[PETSC_MAX_PATH_LEN];
    char path_mesh[PETSC_MAX_PATH_LEN];
    int label_bound = 0, label_omega = 0, order_rbm = 0;
    int mla_phase = 0, gcr_restart = 0;
    double mla_rtol = 0.;
    int mla_max_it = 0;
    int mla_v_pre_smooth = 0, mla_v_post_smooth = 0;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-file_mat", path_mat, sizeof(path_mat), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Matrix file: %s\n", path_mat));
    }

    PetscCall(PetscOptionsGetString(NULL, NULL, "-file_rhs", path_rhs, sizeof(path_rhs), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "RHS file: %s\n", path_rhs));
    }

    PetscCall(PetscOptionsGetString(NULL, NULL, "-mesh", path_mesh, sizeof(path_mesh), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Mesh file: %s\n", path_mesh));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-label_omega", &label_omega, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "label_omega = %d\n", label_omega));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-label_bound", &label_bound, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "label_bound = %d\n", label_bound));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-order_rbm", &order_rbm, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "order_rbm = %d\n", order_rbm));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-mla_phase", &mla_phase, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mla_phase = %d\n", mla_phase));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-gcr_restart", &gcr_restart, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "gcr_restart = %d\n", gcr_restart));
    }

    PetscCall(PetscOptionsGetReal(NULL, NULL, "-mla_rtol", &mla_rtol, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mla_rtol = %le\n", mla_rtol));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-mla_max_it", &mla_max_it, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mla_max_it = %d\n", mla_max_it));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-mla_v_pre_smooth", &mla_v_pre_smooth, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mla_v_pre_smooth time = %d\n", mla_v_pre_smooth));
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-mla_v_post_smooth", &mla_v_post_smooth, &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "mla_v_post_smooth time = %d\n", mla_v_post_smooth));
    }

    MySolver mysolver;
    SolverPetscInitialize(path_mat, path_rhs, &mysolver);

    // ==== step 2, aggregation-based technique implementation ==== //
    /*
     * fine mesh graph struct
     * aggregation mesh graph struct
     * coarse mesh graph struct
     */

    // mesh file process
    GenericList data_list_phy_tag, data_list_node;
    GenericList data_list_ele_bound, data_list_ele_omega;
    InitializeList(&data_list_phy_tag);
    InitializeList(&data_list_node);
    InitializeList(&data_list_ele_bound);
    InitializeList(&data_list_ele_omega);

    MeshFileProcess(path_mesh, label_bound, label_omega,
                    &data_list_phy_tag, &data_list_node,
                    &data_list_ele_bound, &data_list_ele_omega);

    MeshGraph *graph = CreateMeshGraph(data_list_node.size);
    AssembleMeshGraph(graph, &data_list_node, &data_list_ele_omega);
#if 1
    PrintMeshGraph(graph);
#endif

    // aggregation
    MeshGraph *graph_tmp = CreateMeshGraph(graph->size);
    CopyMeshGraph(graph_tmp, graph);

    // puts("\nCopied graph:");
    // PrintMeshGraph(graph_tmp);

    MeshGraph *graph_aggregation = AggregationMeshGraph(graph_tmp);
    // puts("\nAggregation graph:");
    // PrintMeshGraph(graph_tmp);

#if 1
    puts("\nUpdating aggregation graph:");
    PrintMeshGraph(graph_aggregation);
#endif

    // coarse level processing
    MeshGraph *coarse_graph = CreateMeshGraph(graph_aggregation->size);
    MeshNode *data_coarse_node = (MeshNode *)malloc(graph_aggregation->size * sizeof(MeshNode));
    assert(data_coarse_node);

    AssembleCoarseMeshNode(graph_aggregation, data_coarse_node);

    AssembleCoarseMeshGraph(coarse_graph, graph_aggregation, graph, data_coarse_node);
#if 1
    puts("\n\nCoarse Graph:");
    PrintMeshGraph(coarse_graph);
#endif

    // step 3, multilevel iteration
    /*
     * mesh processing transform to graph data struct
     *     1. fine graph
     *     2. aggregation graph
     *     3. coarse graph
     */
    MLAGraph mla; // multilevel aggregation variable
    mla.fine = graph;
#if 1
    mla.aggregation = graph_aggregation;
    mla.coarse = coarse_graph;
#endif // classic aggregation techinique
#if 0
    mla.aggregation = graph;
    mla.coarse = graph;
#endif // aggregation with 1 node
    mla.prolongation_set = 0;
    MLAIterationSolver(&mysolver, &mla, mla_phase, gcr_restart,
                       mla_rtol, mla_max_it,
                       mla_v_pre_smooth, mla_v_post_smooth,
                       order_rbm);

    // computing residual
    SolverPetscResidualCheck(&mysolver);

    // free memory
    ClearMeshGraph(coarse_graph);
    ClearMeshGraph(graph_tmp);
    ClearMeshGraph(graph_aggregation);
    ClearMeshGraph(graph);
    ClearList(&data_list_phy_tag);
    ClearList(&data_list_node);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_ele_omega);

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}

// command line
/*
 * ./app_petsc_exe -file_mat </path/to/petsc-bin-mat> -file_rhs </path/to/petsc-bin-rhs>
 *     -mesh </path/to/mesh/file>
 *     -label_bound 1
 *     -label_omega 2
 *     -order_rbm <1 or 2>
 *     -mla_phase <0 or 1 or 2>
 *     -gcr_restart <int num>
 *     -mla_rtol <double num>
 *     -mla_max_it <int num>
 *     -mla_v_pre_smooth <int num>
 *     -mla_v_post_smooth <int num>
 */
