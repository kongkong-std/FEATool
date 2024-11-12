#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct mesh_node
{
    int idx;
    int add_flag;
} MeshNode;

typedef struct mesh_adj_node
{
    MeshNode *node;
    struct mesh_adj_node *next;
} MeshAdjNode;

typedef struct mesh_adj_list
{
    /* data */
    MeshAdjNode *head;
    int size;
} MeshAdjList;

typedef struct mesh_graph
{
    int size;           // number of vertices
    MeshAdjList *array; // adjacency list array
} MeshGraph;

MeshGraph *CreateMeshGraph(int size);
void ClearMeshGraph(MeshGraph *graph);
void AddEdgeMeshGraph(MeshGraph *graph, MeshNode *node, int u, int v);
void PrintMeshGraph(MeshGraph *graph);
void InitializeMeshAdjList(MeshAdjList *);
void AddNodeMeshAdjList(MeshAdjList *, MeshNode *);
void PrintMeshAdjList(MeshAdjList *);
void ClearMeshAdjList(MeshAdjList *);

int main(int argc, char **argv)
{
    int size = 0;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-size", argv[index]))
        {
            size = atoi(argv[index + 1]);
        }
    }

    // mesh node data
    MeshNode *data_node = (MeshNode *)malloc(size * sizeof(MeshNode));
    assert(data_node);
    for (int index = 0; index < size; ++index)
    {
        data_node[index].idx = index + 1; // 1-base
        data_node[index].add_flag = 0;    // unadded
    }

    MeshGraph *graph = CreateMeshGraph(size);

    // adding edge
    AddEdgeMeshGraph(graph, data_node, 1, 2);
    AddEdgeMeshGraph(graph, data_node, 1, 4);
    AddEdgeMeshGraph(graph, data_node, 2, 3);
    AddEdgeMeshGraph(graph, data_node, 2, 5);
    AddEdgeMeshGraph(graph, data_node, 3, 6);
    AddEdgeMeshGraph(graph, data_node, 4, 5);
    AddEdgeMeshGraph(graph, data_node, 4, 7);
    AddEdgeMeshGraph(graph, data_node, 5, 6);
    AddEdgeMeshGraph(graph, data_node, 6, 10);
    AddEdgeMeshGraph(graph, data_node, 6, 8);
    AddEdgeMeshGraph(graph, data_node, 7, 8);
    AddEdgeMeshGraph(graph, data_node, 7, 11);
    AddEdgeMeshGraph(graph, data_node, 8, 9);
    AddEdgeMeshGraph(graph, data_node, 8, 12);
    AddEdgeMeshGraph(graph, data_node, 9, 10);
    AddEdgeMeshGraph(graph, data_node, 9, 13);
    AddEdgeMeshGraph(graph, data_node, 11, 12);
    AddEdgeMeshGraph(graph, data_node, 11, 14);
    AddEdgeMeshGraph(graph, data_node, 12, 13);
    AddEdgeMeshGraph(graph, data_node, 13, 15);

    // print graph
    PrintMeshGraph(graph);

    // free memory
    ClearMeshGraph(graph);
    free(data_node);

    return 0;
}

void ClearMeshGraph(MeshGraph *graph)
{
    for(int index = 0; index < graph->size; ++index)
    {
        ClearMeshAdjList(graph->array + index);
    }
    free(graph->array);
    free(graph);
}

void ClearMeshAdjList(MeshAdjList *list)
{
    MeshAdjNode * curr = list->head;
    while(curr)
    {
        MeshAdjNode * tmp = curr;
        free(tmp);
        curr = curr->next;
    }
}

void PrintMeshAdjList(MeshAdjList *list)
{
    printf("length of list: %d\n", list->size);
    MeshAdjNode *curr = list->head;
    while (curr)
    {
        printf("%d->", curr->node->idx);
        curr = curr->next;
    }
    puts("NULL");
}

void PrintMeshGraph(MeshGraph *graph)
{
    printf("number of vertex(s): %d\n", graph->size);
    for (int index = 0; index < graph->size; ++index)
    {
        PrintMeshAdjList(graph->array + index);
    }
}

void AddEdgeMeshGraph(MeshGraph *graph, MeshNode *node, int u, int v)
{
    // 1-base to 0-base
    --u;
    --v;

    // add node to adjacency list, node u adjacency list
    AddNodeMeshAdjList(graph->array + u, node + u);
    AddNodeMeshAdjList(graph->array + u, node + v);

    // add node to adjacency list, node v adjacency list
    AddNodeMeshAdjList(graph->array + v, node + v);
    AddNodeMeshAdjList(graph->array + v, node + u);
}

void AddNodeMeshAdjList(MeshAdjList *list, MeshNode *node)
{
#if 1
    // Check if the node is already in the adjacency list
    MeshAdjNode *current = list->head;
    while (current != NULL)
    {
        if (current->node == node) // Node is already in the list
        {
            return; // Avoid adding the same node again
        }
        current = current->next;
    }
#endif

    MeshAdjNode *new_node = (MeshAdjNode *)malloc(sizeof(MeshAdjNode));
    assert(new_node);

    // new_node data
    new_node->node = node;
    new_node->next = NULL;

    if (list->head == NULL)
    {
        list->head = new_node;
        ++(list->size);
        return;
    }

    MeshAdjNode *last_node = list->head;
    while (last_node->next != NULL)
    {
        last_node = last_node->next;
    }
    last_node->next = new_node;
    ++(list->size);
}

void InitializeMeshAdjList(MeshAdjList *list)
{
    list->head = NULL;
    list->size = 0;
}

MeshGraph *CreateMeshGraph(int size)
{
    MeshGraph *graph = (MeshGraph *)malloc(sizeof(MeshGraph));
    assert(graph);

    graph->size = size;
    graph->array = (MeshAdjList *)malloc(size * sizeof(MeshAdjList));
    assert(graph->array);

    // initialize graph
    for (int index = 0; index < size; ++index)
    {
        InitializeMeshAdjList(graph->array + index);
    }

    return graph;
}
