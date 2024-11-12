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

#if 0
/*
 * mesh node RB node
 */
typedef struct mesh_node_rb_node
{
    /* data */
    MeshNode data;
    int color;
    struct mesh_node_rb_node *left, *right, *parent;
} MeshNodeRBNode;

/*
 * mesh node RB tree
 */
typedef struct mesh_node_rb_tree
{
    /* data */
    MeshNodeRBNode *root;  // root node
    MeshNodeRBNode *TNULL; // nil node
} MeshNodeRBTree;
#endif // rb tree

// function prototype
#if 0
/*
 * search data in RB tree
 */
MeshNodeRBNode *SearchRBTree(MeshNodeRBTree *, int);

/*
 * insert data to RB tree
 */
void InsertRBTree(MeshNodeRBTree *, MeshNode);

/*
 * fix RB tree holds the property
 */
void FixInsert(MeshNodeRBTree *, MeshNodeRBNode *);

/*
 * RB tree left rotate
 */
void LeftRotate(MeshNodeRBTree *, MeshNodeRBNode *);

/*
 * RB tree right rotate
 */
void RightRotate(MeshNodeRBTree *, MeshNodeRBNode *);

/*
 * create RB node
 */
MeshNodeRBNode *CreateMeshNodeRBNode(MeshNode);

/*
 * create RB tree
 */
MeshNodeRBTree *CreateMeshNodeRBTree(void);
#endif // rb tree function

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