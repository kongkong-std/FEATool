#include "../include/main.h"

void AddEdgeMeshGraphWitNode(MeshGraph *graph, MeshNode *node, int u, int v)
{
    // 1-base to 0-base
    --u;
    --v;

    // add node to adjacency list, node u adjacency list
    AddNodeMeshGraphAdjList(graph->array + u, node + u);
    AddNodeMeshGraphAdjList(graph->array + u, node + v);

    // add node to adjacency list, node v adjacency list
    AddNodeMeshGraphAdjList(graph->array + v, node + v);
    AddNodeMeshGraphAdjList(graph->array + v, node + u);
}

int IsConnectAggregationMeshGraph(MeshGraph *graph_aggregation, MeshGraph *graph, int u, int v)
{
    MeshGraphAdjList *graph_aggregation_list_u = graph_aggregation->array + u;
    MeshGraphAdjList *graph_aggregation_list_v = graph_aggregation->array + v;

    MeshGraphAdjNode *curr_graph_aggregation_list_u = graph_aggregation_list_u->head;
    while (curr_graph_aggregation_list_u)
    {
        int tmp_u = curr_graph_aggregation_list_u->node->node_idx; // 1-base
        MeshGraphAdjList *graph_list_tmp_u = graph->array + tmp_u - 1;

        MeshGraphAdjNode *curr_graph_aggregation_list_v = graph_aggregation_list_v->head;
        while (curr_graph_aggregation_list_v)
        {
            int tmp_v = curr_graph_aggregation_list_v->node->node_idx; // 1-base
            MeshGraphAdjNode *curr_graph_tmp_u = graph_list_tmp_u->head;
            while (curr_graph_tmp_u)
            {
                if (tmp_v == curr_graph_tmp_u->node->node_idx)
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
                AddEdgeMeshGraphWitNode(coarse_graph, coarse_node, u + 1, v + 1);
            }
        }
    }
}

void AssembleCoarseMeshNode(MeshGraph *graph, MeshNode *node)
{
    for (int index = 0; index < graph->size; ++index)
    {
        double tmp_x = 0., tmp_y = 0., tmp_z = 0.;
        MeshGraphAdjNode *curr = graph->array[index].head;
        while (curr)
        {
            tmp_x += curr->node->node_x;
            tmp_y += curr->node->node_y;
            tmp_z += curr->node->node_z;

            curr = curr->next;
        }
        node[index].node_idx = index + 1;
        node[index].add_flag = 0;
        node[index].node_x = tmp_x / graph->array[index].size;
        node[index].node_y = tmp_y / graph->array[index].size;
        node[index].node_z = tmp_z / graph->array[index].size;
    }
}

void DeleteNodeMeshGraphAdjList(MeshGraphAdjList *list, MeshGraphAdjNode *node)
{
    MeshGraphAdjNode *to_delete = node;
    if (node == list->head)
    {
        list->head = node->next;
    }
    else
    {
        MeshGraphAdjNode *prev = list->head;
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
#if 1
    int cnt_aggregation = 0;
    for (int index = 0; index < graph->size; ++index)
    {
        MeshGraphAdjNode *curr_head_node = graph->array[index].head;

        // 如果当前顶点没有被访问过（add_flag == 0）
        if (curr_head_node->node->add_flag == 0)
        {
            ++cnt_aggregation; // 增加聚类计数
            // 遍历该节点的所有邻接节点
            while (curr_head_node)
            {
                if (curr_head_node->node->add_flag == 1)
                {
                    DeleteNodeMeshGraphAdjList(graph->array + index, curr_head_node);

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

    MeshGraph *graph_aggregation = CreateMeshGraph(cnt_aggregation);
    int cnt = 0;
    for (int index = 0; index < graph->size; ++index)
    {
        if (graph->array[index].head != NULL)
        {
            MeshGraphAdjNode *graph_node = graph->array[index].head;
            while (graph_node)
            {
                AddNodeMeshGraphAdjList(graph_aggregation->array + cnt, graph_node->node);
                graph_node = graph_node->next;
            }
            ++cnt;
        }
    }
#endif // aggregation technique

#if 0
    MeshGraph *graph_aggregation = CreateMeshGraph(graph->size);
    for (int index = 0; index < graph->size; ++index)
    {
        MeshGraphAdjNode *graph_node = graph->array[index].head;
        AddNodeMeshGraphAdjList(graph_aggregation->array + index, graph_node->node);
    }
#endif // every single node as a aggregation body

    return graph_aggregation;
}

void CopyMeshGraph(MeshGraph *graph_dst, MeshGraph *graph_src)
{
    for (int index = 0; index < graph_src->size; ++index)
    {
        MeshGraphAdjNode *src_adj_node = graph_src->array[index].head;
        while (src_adj_node)
        {
            AddNodeMeshGraphAdjList(graph_dst->array + index, src_adj_node->node);

            src_adj_node = src_adj_node->next;
        }
    }
}

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