#include "../include/main.h"

int main(int argc, char **argv)
{
    char *path_mesh = NULL;
    int label_bound = 0, label_omega = 0;
    int searchIdx = 0;
    int order_rbm = 0;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_mesh = argv[index + 1];
        }
        if (strstr("-label_bound", argv[index]))
        {
            label_bound = atoi(argv[index + 1]);
        }
        if (strstr("-label_omega", argv[index]))
        {
            label_omega = atoi(argv[index + 1]);
        }
        if (strstr("-index", argv[index]))
        {
            searchIdx = atoi(argv[index + 1]);
        }
        if (strstr("-order_rbm", argv[index]))
        {
            order_rbm = atoi(argv[index + 1]);
        }
    }

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

    // print data linked list
    printf("==== data linked list ====\n");
    printf("size physical tag: %d\n", data_list_phy_tag.size);
    printf("size node: %d\n", data_list_node.size);
    printf("size bound element: %d\n", data_list_ele_bound.size);
    printf("size omega element: %d\n", data_list_ele_omega.size);

    ListNode *current_mesh_node = data_list_node.head;
    while (current_mesh_node != NULL)
    {
        MeshNode *tmp_mesh_nodex = (MeshNode *)current_mesh_node->data;
        if (tmp_mesh_nodex->node_idx == searchIdx)
        {
            printf("node index %d: %021.16le\t%021.16le\t%021.16le\n",
                   tmp_mesh_nodex->node_idx,
                   tmp_mesh_nodex->node_x,
                   tmp_mesh_nodex->node_y,
                   tmp_mesh_nodex->node_z);
            break;
        }
        current_mesh_node = current_mesh_node->next;
    }

    MeshGraph *graph = CreateMeshGraph(data_list_node.size);
    AssembleMeshGraph(graph, &data_list_node, &data_list_ele_omega);
    PrintMeshGraph(graph);

    // aggregation
    MeshGraph *graph_tmp = CreateMeshGraph(graph->size);
    CopyMeshGraph(graph_tmp, graph);

    puts("\nCopied graph:");
    PrintMeshGraph(graph_tmp);

    puts("\nAggregation graph:");
    MeshGraph *graph_aggregation = AggregationMeshGraph(graph_tmp);
    PrintMeshGraph(graph_tmp);

    puts("\nUpdating aggregation graph:");
    PrintMeshGraph(graph_aggregation);

    // coarse level processing
    MeshGraph *coarse_graph = CreateMeshGraph(graph_aggregation->size);
    MeshNode *data_coarse_node = (MeshNode *)malloc(graph_aggregation->size * sizeof(MeshNode));
    assert(data_coarse_node);

    AssembleCoarseMeshNode(graph_aggregation, data_coarse_node);

    AssembleCoarseMeshGraph(coarse_graph, graph_aggregation, graph, data_coarse_node);
    puts("\n\nCoarse Graph:");
    PrintMeshGraph(coarse_graph);

    // prolongation operator
    double ***P_operator = (double ***)malloc(graph_aggregation->size * sizeof(double **));
    assert(P_operator);
    for (int index = 0; index < graph_aggregation->size; ++index)
    {
        *(P_operator + index) = (double **)malloc(graph_aggregation->array[index].size * sizeof(double *));
        assert(*(P_operator + index));
        for (int index_i = 0; index_i < graph_aggregation->array[index].size; ++index_i)
        {
            if (order_rbm == 1)
            {
                // local prolongation operator
                /*
                 * [I3  R1]
                 * [O3  I3]
                 */
                *(*(P_operator + index) + index_i) = (double *)malloc(15 * sizeof(double));
                assert(*(*(P_operator + index) + index_i));
            }
            else if (order_rbm == 2)
            {
                // local prolongation operator
                /*
                 * [I3  R1  R2]
                 * [O3  I3  O3]
                 */
                *(*(P_operator + index) + index_i) = (double *)malloc(24 * sizeof(double));
                assert(*(*(P_operator + index) + index_i));
            }
            else
            {
                fprintf(stderr, "invalid rbm order, 1 or 2!\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    ConstructProlongationOperator(P_operator, graph_aggregation, coarse_graph, order_rbm);
    puts("\nPrint prolongation operator:");
    PrintProlongationOperator(P_operator, graph_aggregation, order_rbm);

    int *map_row = (int *)malloc(graph->size * sizeof(int));
    int *map_column = (int *)malloc(graph->size * sizeof(int));
    assert(map_row && map_column);

    puts("\n==== maping aggregation to prolongation block row and column:");
    MapAggregationProlongationOperator(graph_aggregation, map_row, map_column, graph->size);
    for (int index = 0; index < graph->size; ++index)
    {
        printf("%d\t%d\n", map_row[index], map_column[index]);
    }

    // free memory
    free(map_row);
    free(map_column);
    for (int index = 0; index < graph_aggregation->size; ++index)
    {
        for (int index_i = 0; index_i < graph_aggregation->array[index].size; ++index_i)
        {
            free(P_operator[index][index_i]);
        }
        free(P_operator[index]);
    }
    free(P_operator);
    free(data_coarse_node);
    ClearMeshGraph(coarse_graph);
    ClearMeshGraph(graph_tmp);
    ClearMeshGraph(graph_aggregation);
    ClearMeshGraph(graph);
    ClearList(&data_list_phy_tag);
    ClearList(&data_list_node);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_ele_omega);

    return 0;
}

// command line
/*
 * ./app_exe -mesh <path/to/mesh/file>
 *     -label_bound 1
 *     -label_omega 2
 */
