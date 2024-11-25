#ifndef MAIN_H_
#define MAIN_H_

// header
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <petscksp.h>

#define MAX_SIZE 1024

// data struct
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

typedef struct my_solver
{
    KSP ksp;
    PC pc;
    Mat solver_a;
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
} MySolver;

typedef struct mla_graph
{
    /* data */
    MeshGraph *fine;
    MeshGraph *aggregation;
    MeshGraph *coarse;
    Mat prolongation, operator_coarse;
    int prolongation_set; // 0 represents unset, 1 represents setting
} MLAGraph;

// function prototype
/*
 * mla solve phase
 */
void MLASolvePhase(MySolver *mysolver, MLAGraph *mla, int gcr_restart);

/*
 * mla post-smooth phase, gauss-seidel smoother
 */
void MLAPostSmoothPhase(MySolver *mysolver, int v_post_smooth);

/*
 * mla pre-smooth phase, gauss-seidel smoother
 */
void MLAPreSmoothPhase(MySolver *mysolver, int v_pre_smooth);

/*
 * mla setup phase constructs prolongation operator
 */
void MLASetupPhase(MySolver *, MLAGraph *, int /*rbm order*/);

/*
 * computing relative residual of mla
 */
void MLARelativeResidual(MySolver *, double *);

/*
 * multilevel iteration solver
 *     GS smoother
 *     mla_phase: 0 setup, 1 solve, 2 setup + solve
 */
void MLAIterationSolver(MySolver *, MLAGraph *, int /*mla phase*/, int /*gcr restart*/,
                        double /*tolerance*/, int /*max iteration*/,
                        int /*pre-smoothing time*/, int /*post-smoothing time*/,
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
 * center of aggregation is choosen as coarse node
 */
void AssembleCoarseMeshNode(MeshGraph * /*coarse graph*/, MeshNode * /*coarse node*/);

/*
 * coarse graph adjacency list represents and node is coarse node
 */
void AssembleCoarseMeshGraph(MeshGraph * /*coarse graph*/, MeshGraph * /*graph aggregation*/,
                             MeshGraph * /*fine graph*/, MeshNode * /*coarse node data*/);

/*
 * delete specific node in adjacency list
 */
void DeleteNodeMeshGraphAdjList(MeshGraphAdjList *, MeshGraphAdjNode *);

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
 * clean graph
 */
void ClearMeshGraph(MeshGraph *);

/*
 * assembling graph with element data
 */
void AssembleMeshGraph(MeshGraph *, GenericList *, GenericList *);

/*
 * print graph
 */
void PrintMeshGraph(MeshGraph *);

/*
 * add edge to a graph
 */
void AddEdgeMeshGraph(MeshGraph *, GenericList * /*node data*/, int, int);

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
 * create a graph with number of vertices
 */
MeshGraph *CreateMeshGraph(int);

/*
 * check whether linked list is empty or not
 *     0: not empty
 *     1: empty
 */
int IsListEmpty(GenericList * /*linked list*/);

/*
 * add node to the list
 */
void AddNodeToList(GenericList * /*linked list*/, void * /*data*/, NodeType /*node type*/);

/*
 * clear linked list
 */
void ClearList(GenericList * /*linked list*/);

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

void SolverPetscInitialize(char *path_mat, char *path_rhs, MySolver *mysolver);
void SolverPetscResidualCheck(MySolver *mysolver);

#endif // main.h