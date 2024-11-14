#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// macro
#define MAX_SIZE 1024
// #define RED 0   // RB tree
// #define BLACK 1 // RB tree

// struct, linked list
typedef enum
{
    TYPE_PHYSICAL_TAG,
    TYPE_NODE,
    TYPE_ELEMENT
} NodeType;

typedef struct
{
    /* data */
    int tag_dim; // tag dimensional
    int tag_idx; // tag index
} MeshPhysicalTag;

typedef struct list_node
{
    void *data;             // 指向数据的指针
    struct list_node *next; // 指向下一个节点的指针
    NodeType type;
} ListNode;

typedef struct generic_list
{
    ListNode *head; // 链表头
    ListNode *tail; // 链表尾
    int size;       // 链表大小
} GenericList;

typedef struct
{
    /* data */
    int node_idx;                  // node index
    int add_flag;                  // 0 represents unvisited, 1 represents visited
    double node_x, node_y, node_z; // node coordinate
} MeshNode;

typedef struct
{
    /* data */
    int ele_idx;       // element index
    int ele_type;      // element type
    int ele_num_tag;   // element tag number
    int *ele_tag_idx;  // element tag index, size ele_num_tag
    int *ele_node_lst; // element node list, size map(ele_type)
} MeshElement;

typedef struct mesh_graph_adj_node
{
    /* data */
    MeshNode *node;
    struct mesh_graph_adj_node *next;
} MeshGraphAdjNode;

typedef struct mesh_graph_adj_list
{
    /* data */
    MeshGraphAdjNode *head;
    int size;
} MeshGraphAdjList;

typedef struct mesh_graph
{
    /* data */
    int size; // number of vertices
    MeshGraphAdjList *array;
} MeshGraph;

// function prototype
/*
 * mapping of aggregation block to prolognation operator dof ordered by node indices,
 * to prolongation block row and clolumn
 */
void MapAggregationProlongationOperator(MeshGraph * /*graph aggreagtion*/,
                                        int * /*map_row*/,
                                        int * /*map_column*/,
                                        int /*size*/); // block row and column

/*
 * print prolongation operator
 */
void PrintProlongationOperator(double *** /*prolongation operator*/,
                               MeshGraph * /*graph aggregation*/, int);

/*
 * prolongation operator construction
 */
void ConstructProlongationOperator(double *** /*prolongation operator*/,
                                   MeshGraph * /*graph aggreation*/,
                                   MeshGraph * /*coarse graph*/,
                                   int /*rbm order*/);

/*
 * add node data node_u, node_v to graph
 */
void AddEdgeMeshGraphWitNode(MeshGraph *graph, MeshNode *node, int u, int v);

/*
 * check aggragation is connected or not
 */
int IsConnectAggregationMeshGraph(MeshGraph * /*graph aggregation*/, MeshGraph * /*fine graph*/,
                                  int /*adjacency list u*/, int /*adjacency list v*/);

/*
 * coarse graph adjacency list represents and node is coarse node
 */
void AssembleCoarseMeshGraph(MeshGraph * /*coarse graph*/, MeshGraph * /*graph aggregation*/,
                             MeshGraph * /*fine graph*/, MeshNode * /*coarse node data*/);

/*
 * center of aggregation is choosen as coarse node
 */
void AssembleCoarseMeshNode(MeshGraph * /*coarse graph*/, MeshNode * /*coarse node*/);

/*
 * delete specific node in adjacency list
 */
void DeleteNodeMeshAdjList(MeshGraphAdjList *, MeshGraphAdjNode *);

/*
 * fine graph aggregation, deleting adjacency list
 */
MeshGraph *AggregationMeshGraph(MeshGraph * /*fine graph*/);

/*
 * copying a graph to another graph
 */
void CopyMeshGraph(MeshGraph * /*destination graph*/, MeshGraph * /*source graph*/);

/*
 * clear adjacency list
 */
void ClearMeshAdjList(MeshGraphAdjList *);

/*
 * assembling graph with element data
 */
void AssembleMeshGraph(MeshGraph *, GenericList *, GenericList *);

/*
 * print adjacency list
 */
void PrintMeshGraphAdjList(MeshGraphAdjList *);

/*
 * add node to adjacency list
 */
void AddNodeMeshGraphAdjList(MeshGraphAdjList *, MeshNode *);

/*
 * initialize adjacency list
 */
void InitializeMeshGraphAdjList(MeshGraphAdjList *);

/*
 * print graph
 */
void PrintMeshGraph(MeshGraph *);

/*
 * add edge to a graph
 */
void AddEdgeMeshGraph(MeshGraph *, GenericList * /*node data*/, int, int);

/*
 * clean graph
 */
void ClearMeshGraph(MeshGraph *);

/*
 * create a graph with number of vertices
 */
MeshGraph *CreateMeshGraph(int);

/*
 * clear linked list
 */
void ClearList(GenericList * /*linked list*/);

/*
 * add node to the list
 */
void AddNodeToList(GenericList * /*linked list*/, void * /*data*/, NodeType /*node type*/);

/*
 * check whether linked list is empty or not
 *     0: not empty
 *     1: empty
 */
int IsListEmpty(GenericList * /*linked list*/);

/*
 * initialize linked list
 */
void InitializeList(GenericList * /*linked list*/);

/*
 * mapping of element type to node of element
 */
int NodeListSizeMapping(int);

/*
 * mesh file process
 */
void MeshFileProcess(const char * /*path to mesh file*/,
                     int /*label of bound*/, int /*label of omega*/,
                     GenericList * /*generic data list*/,
                     GenericList * /*mesh node data*/,
                     GenericList * /*mesh bound element data*/,
                     GenericList * /*mesh omega element data*/);

#endif // main.h