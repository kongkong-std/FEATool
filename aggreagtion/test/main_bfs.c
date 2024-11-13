#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_VERTICES 100

// 邻接表结构体
typedef struct
{
    int adj[MAX_VERTICES][MAX_VERTICES]; // 邻接矩阵表示的图
    int V;                               // 顶点数量
} Graph;

// 队列结构体
typedef struct
{
    int items[MAX_VERTICES];
    int front;
    int rear;
} Queue;

// 队列操作函数
void initQueue(Queue *q)
{
    q->front = -1;
    q->rear = -1;
}

bool isEmpty(Queue *q)
{
    return q->front == -1;
}

void enqueue(Queue *q, int value)
{
    if (q->rear == MAX_VERTICES - 1)
        return; // 队列满
    if (q->front == -1)
        q->front = 0;
    q->items[++(q->rear)] = value;
}

int dequeue(Queue *q)
{
    if (isEmpty(q))
        return -1; // 队列为空
    int item = q->items[q->front];
    if (q->front == q->rear)
        q->front = q->rear = -1; // 队列变为空
    else
        q->front++;
    return item;
}

// 图的操作函数
void initGraph(Graph *g, int V)
{
    g->V = V;
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
        {
            g->adj[i][j] = 0; // 初始化为无边
        }
    }
}

void addEdge(Graph *g, int u, int v)
{
    g->adj[u][v] = 1;
    g->adj[v][u] = 1; // 无向边
}

// BFS实现聚类
void BFS(Graph *g, int start, bool visited[], int *cluster, int *cluster_size)
{
    Queue q;
    initQueue(&q);
    enqueue(&q, start);
    visited[start] = true;

    while (!isEmpty(&q))
    {
        int v = dequeue(&q);
        cluster[*cluster_size] = v; // 将当前顶点添加到聚类中
        (*cluster_size)++;

        // 遍历相邻节点
        for (int i = 0; i < g->V; i++)
        {
            if (g->adj[v][i] == 1 && !visited[i])
            { // 如果有边且未访问
                visited[i] = true;
                enqueue(&q, i);
            }
        }
    }
}

// 聚类函数
void findClusters(Graph *g)
{
    bool visited[g->V];
    for (int i = 0; i < g->V; i++)
    {
        visited[i] = false; // 初始化所有顶点未访问
    }

    int clusters[g->V][g->V]; // 存储聚类结果
    int cluster_count = 0;

    // 对每个未访问的顶点启动BFS
    for (int i = 0; i < g->V; i++)
    {
        if (!visited[i])
        {
            int cluster[g->V];
            int cluster_size = 0;
            BFS(g, i, visited, cluster, &cluster_size);

            // 存储当前聚类
            for (int j = 0; j < cluster_size; j++)
            {
                clusters[cluster_count][j] = cluster[j];
            }
            cluster_count++;

            // 打印当前聚类
            printf("Cluster %d: ", cluster_count);
            for (int j = 0; j < cluster_size; j++)
            {
                printf("%d ", cluster[j]);
            }
            printf("\n");
        }
    }
}

int main()
{
    Graph g;
    initGraph(&g, 7); // 假设有7个顶点

    // 添加一些边
    addEdge(&g, 0, 1);
    addEdge(&g, 1, 2);
    addEdge(&g, 3, 4);
    addEdge(&g, 4, 5);
    addEdge(&g, 0, 3);
    addEdge(&g, 0, 2);
    addEdge(&g, 0, 6);
    addEdge(&g, 6, 2);
    addEdge(&g, 6, 3);
    addEdge(&g, 6, 5);
    addEdge(&g, 5, 3);

    // 找到并输出图的聚类
    findClusters(&g);

    return 0;
}
