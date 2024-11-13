#include "../include/main.h"

void ClearMeshAdjList(MeshGraphAdjList *list)
{
    MeshGraphAdjNode *curr = list->head;
    while (curr)
    {
        MeshGraphAdjNode *tmp = curr;
        free(tmp);
        curr = curr->next;
    }
    list->size = 0;
    list->head = NULL;
}

void ClearMeshGraph(MeshGraph *graph)
{
    for (int index = 0; index < graph->size; ++index)
    {
        ClearMeshAdjList(graph->array + index);
    }

    free(graph->array);
    free(graph);
}

void AssembleMeshGraph(MeshGraph *graph, GenericList *data_node, GenericList *data_ele)
{
    ListNode *curr_ele = data_ele->head;
    while (curr_ele)
    {
        MeshElement *curr_ele_data = (MeshElement *)curr_ele->data;
        int ele_num = NodeListSizeMapping(curr_ele_data->ele_type);
        for (int index = 0; index < ele_num; ++index)
        {
            AddEdgeMeshGraph(graph, data_node,
                             curr_ele_data->ele_node_lst[index],
                             curr_ele_data->ele_node_lst[(index + 1) % ele_num]);
        }

        curr_ele = curr_ele->next;
    }
}

void PrintMeshGraph(MeshGraph *graph)
{
    printf("number of vertex(s): %d\n", graph->size);
    for (int index = 0; index < graph->size; ++index)
    {
        PrintMeshGraphAdjList(graph->array + index);
    }
}

void AddEdgeMeshGraph(MeshGraph *graph, GenericList *data_list, int u, int v)
{
    // 1-base to 0-base
    --u;
    --v;

    ListNode *curr_node = data_list->head;

    for (int index = 0; index < u; ++index)
    {
        curr_node = curr_node->next;
    }
    MeshNode *node_u = (MeshNode *)curr_node->data;

    curr_node = data_list->head;
    for (int index = 0; index < v; ++index)
    {
        curr_node = curr_node->next;
    }
    MeshNode *node_v = (MeshNode *)curr_node->data;

    // add node to adjacency list, node u adjacency list
    AddNodeMeshGraphAdjList(graph->array + u, node_u);
    AddNodeMeshGraphAdjList(graph->array + u, node_v);

    // add node to adjacency list, node v adjacency list
    AddNodeMeshGraphAdjList(graph->array + v, node_v);
    AddNodeMeshGraphAdjList(graph->array + v, node_u);
}

void PrintMeshGraphAdjList(MeshGraphAdjList *list)
{
    printf("length of list: %d\n", list->size);
    MeshGraphAdjNode *curr = list->head;
    while (curr)
    {
        printf("%d->", curr->node->node_idx);
        curr = curr->next;
    }
    puts("NULL");
}

void AddNodeMeshGraphAdjList(MeshGraphAdjList *list, MeshNode *node)
{
#if 1
    // Check if the node is already in the adjacency list
    MeshGraphAdjNode *current = list->head;
    while (current != NULL)
    {
        if (current->node == node) // Node is already in the list
        {
            return; // Avoid adding the same node again
        }
        current = current->next;
    }
#endif

    MeshGraphAdjNode *new_node = (MeshGraphAdjNode *)malloc(sizeof(MeshGraphAdjNode));
    assert(new_node);

#if 1
    // new_node data
    new_node->node = node;
    new_node->next = NULL;

    if (list->head == NULL)
    {
        list->head = new_node;
        ++(list->size);
        return;
    }

    MeshGraphAdjNode *last_node = list->head;
    while (last_node->next != NULL)
    {
        last_node = last_node->next;
    }
    last_node->next = new_node;
    ++(list->size);
#endif // insert in tail of adjacency list

#if 0
    new_node->node = node;
    new_node->next = list->head;

    list->head = new_node;
    ++(list->size);
#endif // insert in head of adjacency list
}

void InitializeMeshGraphAdjList(MeshGraphAdjList *list)
{
    list->head = NULL;
    list->size = 0;
}

MeshGraph *CreateMeshGraph(int size)
{
    MeshGraph *graph = (MeshGraph *)malloc(sizeof(MeshGraph));
    assert(graph);

    graph->size = size;
    graph->array = (MeshGraphAdjList *)malloc(size * sizeof(MeshGraphAdjList));
    assert(graph->array);

    // initialize graph
    for (int index = 0; index < size; ++index)
    {
        InitializeMeshGraphAdjList(graph->array + index);
    }

    return graph;
}
