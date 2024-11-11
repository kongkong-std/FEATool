#include "../include/main.h"

void PrintMeshGraph(MeshGraph *graph)
{
    printf("size of graph :%d\n", graph->size);

    if (!graph)
    {
        printf("Graph is NULL\n");
        return;
    }

    // 遍历图中每个节点
    for (int i = 0; i < graph->size; ++i)
    {
        MeshNode *current_node = &graph->node[i];

        // 打印当前节点的信息
        printf("Node %d:\n",
               current_node->node_idx);

        // 打印该节点的所有邻接节点
        MeshGraphAdjList *adj_list = graph->adjlist[i];
        while (adj_list)
        {
            MeshNode *neighbor_node = adj_list->node;
            printf("\t->%d",
                   neighbor_node->node_idx);
            adj_list = adj_list->next;
        }
        putchar('\n');
    }
}

void CreateEdgeMeshGraph(MeshGraph *graph, GenericList *data_list)
{
    ListNode *curr = data_list->head;
    while (curr)
    {
        MeshElement *tmp = (MeshElement *)curr->data;
        int num_nodes = NodeListSizeMapping(tmp->ele_type);
        for (int index = 0; index < num_nodes; ++index)
        {
            int u = tmp->ele_node_lst[index];
            int v = tmp->ele_node_lst[(index + 1) % num_nodes];
            AddEdgeMeshGraph(graph, u, v);
        }
        curr = curr->next;
    }
}

int EdgeExistsMeshGraph(MeshGraph *graph, int u, int v)
{
    // vertex 1-base
    --u;
    --v;

    MeshGraphAdjList *curr = graph->adjlist[u];
    while (curr)
    {
        if (curr->node->node_idx == v + 1)
        {
            return 1; // 边已经存在
        }
        curr = curr->next;
    }
    return 0; // 边不存在
}

void AddEdgeMeshGraph(MeshGraph *graph, int u, int v)
{
    if (u > graph->size || v > graph->size || u < 0 || v < 0)
    {
        printf("Error: Invalid vertex index\n");
        return;
    }

    if (EdgeExistsMeshGraph(graph, u, v) || EdgeExistsMeshGraph(graph, v, u))
    {
        return; // 如果边已存在，跳过添加
    }

    MeshGraphAdjList *newEdgeU = (MeshGraphAdjList *)malloc(sizeof(MeshGraphAdjList));
    MeshGraphAdjList *newEdgeV = (MeshGraphAdjList *)malloc(sizeof(MeshGraphAdjList));
    assert(newEdgeU && newEdgeV);

    // vertex 1-base
    --u;
    --v;

    // 设置 u -> v 边
    newEdgeU->node = &graph->node[v];   // u 指向 v
    newEdgeU->next = graph->adjlist[u]; // 将当前邻接表挂到新节点的 next 上
    graph->adjlist[u] = newEdgeU;       // 更新 u 的邻接表

    // 设置 v -> u 边
    newEdgeV->node = &graph->node[u];   // v 指向 u
    newEdgeV->next = graph->adjlist[v]; // 将当前邻接表挂到新节点的 next 上
    graph->adjlist[v] = newEdgeV;       // 更新 v 的邻接表
}

void CreateVertexMeshGraph(MeshGraph *graph, GenericList *data_list)
{
    int cnt = 0;
    ListNode *curr = data_list->head;
    while (curr && cnt < graph->size)
    {
        MeshNode *tmp = (MeshNode *)curr->data;
        graph->node[cnt].node_idx = tmp->node_idx;
        graph->node[cnt].node_x = tmp->node_x;
        graph->node[cnt].node_y = tmp->node_y;
        graph->node[cnt].node_z = tmp->node_z;
        curr = curr->next;
        ++cnt;
    }
}

void ClearMeshGraph(MeshGraph *graph)
{
    if (!graph)
        return;

    // 释放邻接表中的每个链表
    for (int i = 0; i < graph->size; ++i)
    {
        MeshGraphAdjList *curr = graph->adjlist[i];
        while (curr)
        {
            MeshGraphAdjList *temp = curr;
            curr = curr->next;
            free(temp); // 释放链表节点
        }
    }

    // 释放节点数组
    free(graph->node);
    // 释放邻接表数组
    free(graph->adjlist);
    // 最后释放图结构体本身
    free(graph);
}

MeshGraph *CreateMeshGraph(int size)
{
    MeshGraph *graph = (MeshGraph *)malloc(sizeof(MeshGraph));
    assert(graph);

    graph->size = size;
    graph->node = (MeshNode *)malloc(size * sizeof(MeshNode));
    graph->adjlist = (MeshGraphAdjList **)malloc(size * sizeof(MeshGraphAdjList *));
    assert(graph->node && graph->adjlist);

    for (int index = 0; index < size; ++index)
    {
        graph->adjlist[index] = NULL;
    }

    // puts("\n==== mesh graph data structure has been constructed...\n");

    return graph;
}
