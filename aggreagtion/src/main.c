#include "../include/main.h"

int main(int argc, char **argv)
{
    char *path_mesh = NULL;
    int label_bound = 0, label_omega = 0;
    int searchIdx = 0;
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

#if 0
    // 创建红黑树
    MeshNodeRBTree *rb_tree = CreateMeshNodeRBTree();
    ListNode *current_mesh_node = data_list_node.head;
    while (current_mesh_node != NULL)
    {
        MeshNode *tmp_mesh_node = (MeshNode *)current_mesh_node->data;
        InsertRBTree(rb_tree, *tmp_mesh_node);
        current_mesh_node = current_mesh_node->next;
    }

    // 查找操作
    MeshNodeRBNode *result = SearchRBTree(rb_tree, searchIdx);
    if (result != rb_tree->TNULL) // 查找成功，返回节点
    {
        printf("Node idx: %d\t(%.6f, %.6f, %.6f)\n",
               result->data.node_idx, result->data.node_x, result->data.node_y, result->data.node_z);
    }
    else
    {
        printf("Node with idx %d not found.\n", searchIdx);
    }
#endif

    MeshGraph *graph = CreateMeshGraph(data_list_node.size);
    CreateVertexMeshGraph(graph, &data_list_node);
    CreateEdgeMeshGraph(graph, &data_list_ele_omega);
    PrintMeshGraph(graph);

    // free memory
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
 *     -index <num>
 */
