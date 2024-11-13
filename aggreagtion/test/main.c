#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct mesh_node
{
    int idx;
    int add_flag;
    double x, y, z;
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
void CopyMeshGraph(MeshGraph *, MeshGraph *);
MeshGraph *AggregationMeshGraph(MeshGraph *);
void DeleteNodeMeshAdjList(MeshAdjList *, MeshAdjNode *);
void AssembleCoarseMeshNode(MeshGraph *, MeshNode *);
void AssembleCoarseMeshGraph(MeshGraph * /*coarse graph*/, MeshGraph * /*graph aggregation*/,
                             MeshGraph * /*fine graph*/, MeshNode * /*coarse node data*/);
int IsConnectAggregationMeshGraph(MeshGraph * /*graph aggregation*/, MeshGraph * /*fine graph*/,
                                  int /*adjacency list u*/, int /*adjacency list v*/);

int main(int argc, char **argv)
{
    int size = 0;
    char *path_node = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-size", argv[index]))
        {
            size = atoi(argv[index + 1]);
        }
        if (strstr("-node", argv[index]))
        {
            path_node = argv[index + 1];
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
    FILE *fp = fopen(path_node, "rb");
    assert(fp);
    for (int index = 0; index < size; ++index)
    {
        fscanf(fp, "%*d%lf%lf%lf", &(data_node[index].x),
               &(data_node[index].y),
               &(data_node[index].z));
    }
    fclose(fp);

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
#if 0
    for (int index = 0; index < graph_aggregation->size; ++index)
    {
        printf("node %d(flag = %d): %021.16le\t%021.16le\t%021.16le\n",
               data_coarse_node[index].idx,
               data_coarse_node[index].add_flag,
               data_coarse_node[index].x,
               data_coarse_node[index].y,
               data_coarse_node[index].z);
    }
#endif // check coarse node

    AssembleCoarseMeshGraph(coarse_graph, graph_aggregation, graph, data_coarse_node);
    puts("\n\nCoarse Graph:");
    PrintMeshGraph(coarse_graph);

    // free memory
    ClearMeshGraph(coarse_graph);
    ClearMeshGraph(graph_aggregation);
    ClearMeshGraph(graph_tmp);
    ClearMeshGraph(graph);
    free(data_coarse_node);
    free(data_node);

    return 0;
}

int IsConnectAggregationMeshGraph(MeshGraph *graph_aggregation, MeshGraph *graph, int u, int v)
{
    MeshAdjList *graph_aggregation_list_u = graph_aggregation->array + u;
    MeshAdjList *graph_aggregation_list_v = graph_aggregation->array + v;

    MeshAdjNode * curr_graph_aggregation_list_u = graph_aggregation_list_u->head;
    while(curr_graph_aggregation_list_u)
    {
        int tmp_u = curr_graph_aggregation_list_u->node->idx;    // 1-base
        MeshAdjList * graph_list_tmp_u = graph->array + tmp_u - 1;

        MeshAdjNode * curr_graph_aggregation_list_v = graph_aggregation_list_v->head;
        while(curr_graph_aggregation_list_v)
        {
            int tmp_v = curr_graph_aggregation_list_v->node->idx;    // 1-base
            MeshAdjNode * curr_graph_tmp_u = graph_list_tmp_u->head;
            while(curr_graph_tmp_u)
            {
                if(tmp_v == curr_graph_tmp_u->node->idx)
                {
                    return 1;
                }

                curr_graph_tmp_u = curr_graph_tmp_u->next;
            }

            curr_graph_aggregation_list_v = curr_graph_aggregation_list_v->next;
        }

        curr_graph_aggregation_list_u = curr_graph_aggregation_list_u->next;
    }

    return 0;
}

void AssembleCoarseMeshGraph(MeshGraph *coarse_graph, MeshGraph *graph_aggregation,
                             MeshGraph *fine_graph, MeshNode *coarse_node)
{
    for (int u = 0; u < graph_aggregation->size - 1; ++u)
    {
        for (int v = u + 1; v < graph_aggregation->size; ++v)
        {
            // check the connection of adjacency_list(u) and adjacency_list(v)
            /*
             * flag = 0 represents unconnected
             * flag = 1 represents connected
             */
            int flag = IsConnectAggregationMeshGraph(graph_aggregation, fine_graph, u, v);
            if (flag == 1)
            {
                // u, v from 0-base to 1-base
                AddEdgeMeshGraph(coarse_graph, coarse_node, u + 1, v + 1);
            }
        }
    }
}

void AssembleCoarseMeshNode(MeshGraph *graph, MeshNode *node)
{
    for (int index = 0; index < graph->size; ++index)
    {
        double tmp_x = 0., tmp_y = 0., tmp_z = 0.;
        MeshAdjNode *curr = graph->array[index].head;
        while (curr)
        {
            tmp_x += curr->node->x;
            tmp_y += curr->node->y;
            tmp_z += curr->node->z;

            curr = curr->next;
        }
        node[index].idx = index + 1;
        node[index].add_flag = 0;
        node[index].x = tmp_x / graph->array[index].size;
        node[index].y = tmp_y / graph->array[index].size;
        node[index].z = tmp_z / graph->array[index].size;
    }
}

void DeleteNodeMeshAdjList(MeshAdjList *list, MeshAdjNode *node)
{
    MeshAdjNode *to_delete = node;
    if (node == list->head)
    {
        list->head = node->next;
    }
    else
    {
        MeshAdjNode *prev = list->head;
        while (prev->next != node)
        {
            prev = prev->next;
        }

        prev->next = node->next;
    }

    free(to_delete);
    --(list->size);
}

MeshGraph *AggregationMeshGraph(MeshGraph *graph)
{
    int cnt_aggregation = 0;
    for (int index = 0; index < graph->size; ++index)
    {
        MeshAdjNode *curr_head_node = graph->array[index].head;

        // 如果当前顶点没有被访问过（add_flag == 0）
        if (curr_head_node->node->add_flag == 0)
        {
            ++cnt_aggregation; // 增加聚类计数
            // 遍历该节点的所有邻接节点
            while (curr_head_node)
            {
                if (curr_head_node->node->add_flag == 1)
                {
#if 0
                    // 删除已访问的节点
                    MeshAdjNode *to_delete = curr_head_node;
                    if (curr_head_node == graph->array[index].head)
                    {
                        // 如果是头节点，直接修改头指针
                        graph->array[index].head = curr_head_node->next;
                    }
                    else
                    {
                        // 找到前驱节点
                        MeshAdjNode *prev = graph->array[index].head;
                        while (prev->next != curr_head_node)
                        {
                            prev = prev->next;
                        }
                        // 删除当前节点
                        prev->next = curr_head_node->next;
                    }
                    // 释放内存
                    free(to_delete);

                    --(graph->array[index].size);
#endif

                    DeleteNodeMeshAdjList(graph->array + index, curr_head_node);

                    // 更新当前节点为下一个节点
                    curr_head_node = curr_head_node->next; // 正常跳到下一个节点
                    continue;                              // 继续处理下一个节点
                }
                else
                {
                    // 标记为已访问
                    curr_head_node->node->add_flag = 1;
                }
                // 更新当前节点
                curr_head_node = curr_head_node->next;
            }
        }
        else
        {
            // 如果当前顶点已经访问过，则清空邻接表
            ClearMeshAdjList(&graph->array[index]);
        }
    }
#if 0
    printf("\n\n==== cnt_aggregation: %d\n", cnt_aggregation);
    for(int index = 0; index < graph->size; ++index)
    {
        if(graph->array[index].head != NULL)
        {
            puts("1");
        }
    }
#endif

    MeshGraph *graph_aggregation = CreateMeshGraph(cnt_aggregation);
    int cnt = 0;
    for (int index = 0; index < graph->size; ++index)
    {
        if (graph->array[index].head != NULL)
        {
            MeshAdjNode *graph_node = graph->array[index].head;
            while (graph_node)
            {
                AddNodeMeshAdjList(graph_aggregation->array + cnt, graph_node->node);
                graph_node = graph_node->next;
            }
            ++cnt;
        }
    }

    return graph_aggregation;
}

void CopyMeshGraph(MeshGraph *graph_dst, MeshGraph *graph_src)
{
    for (int index = 0; index < graph_src->size; ++index)
    {
        MeshAdjNode *src_adj_node = graph_src->array[index].head;
        while (src_adj_node)
        {
            AddNodeMeshAdjList(graph_dst->array + index, src_adj_node->node);

            src_adj_node = src_adj_node->next;
        }
    }
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

void ClearMeshAdjList(MeshAdjList *list)
{
    list->size = 0;
    MeshAdjNode *curr = list->head;
    while (curr)
    {
        MeshAdjNode *tmp = curr;
        free(tmp);
        curr = curr->next;
    }
    list->head = NULL;
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
    puts("NULL\n");
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
