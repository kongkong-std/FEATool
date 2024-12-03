#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <cjson/cJSON.h>

#define MAX_SIZE PETSC_MAX_PATH_LEN

// data struct
/*
 * json data struct
 */
typedef struct config_file
{
    char file_mat[PETSC_MAX_PATH_LEN];
    char file_rhs[PETSC_MAX_PATH_LEN];
    char file_mesh[PETSC_MAX_PATH_LEN];
} ConfigFile;

typedef struct config_mla
{
    int pre_smooth_v;
    int post_smooth_v;
    int mla_max_it;
    double mla_rtol;
    int mla_level;
    int mla_phase;
} ConfigMLA;

typedef struct config_mesh_label
{
    /* data */
    int label_bound, label_omega;
} ConfigMeshLabel;

typedef struct config_json
{
    /* data */
    ConfigFile file_config;
    ConfigMLA mla_config;
    ConfigMeshLabel mesh_label_config;
} ConfigJSON;

/*
 * mesh data struct
 */
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

/*
 * solver data struct
 */
typedef struct my_solver
{
    KSP ksp;
    PC pc;
    Mat solver_a;
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
} MySolver;

/*
 * mla data struct
 */
typedef struct mla_graph
{
    /* data */
    MeshGraph *fine;        // fine mesh
    MeshGraph *aggregation; // aggregation mesh
    MeshGraph *coarse;      // coarse mesh
    MeshGraph *mesh_tmp;    // temporary mesh
    MeshNode *coarse_node;  // coarse mesh node data
    Mat prolongation;       // prolongation operator
    Mat operator_coarse;    // coarse operator
    Mat operator_fine;      // fine operator
    KSP ksp_presmooth;      // pre-smooth solver
    KSP ksp_postsmooth;     // post-smooth solver
    KSP ksp_coarse;         // coarse solver
    int level;              // current level
} MLAGraph;

typedef struct mla_context
{
    /* data */
    int num_level; // total level
    MLAGraph *mla; // neighbouring level, size: num_level - 1
    int setup;     // setup flag, 0 represents un-setup, 1 represents setup
} MLAContext;

// function prototype
/*
 * mla solver coarsest correction phase
 */
void MLASolverCoarsetCorrectionPhase(int order_rbm, KSP ksp, PC pc,
                                     int level,
                                     MLAContext *mla_ctx,
                                     Vec *mg_recur_x,
                                     Vec *mg_recur_b);

/*
 * mla solver nested mg procedure pre-smooth
 */
void MLAMGNestedProcedurePreSmooth(KSP ksp, PC pc,
                                   int level,
                                   MLAContext *mla_ctx,
                                   Vec *mg_recur_x,
                                   Vec *mg_recur_b,
                                   int v_pre_smooth);

/*
 * mla solver nested mg procedure post-smooth
 */
void MLAMGNestedProcedurePostSmooth(KSP ksp, PC pc,
                                    int level,
                                    MLAContext *mla_ctx,
                                    Vec *mg_recur_x,
                                    Vec *mg_recur_b,
                                    int v_post_smooth);

/*
 * mla solver nested mg procedure
 *     level and number of levels
 *     mla context
 *     solution
 *     rhs
 *     pre- and post- smooth times
 *     rbm order
 */
void MLAMGNestedProcedure(int /*level*/, int /*number of levels*/,
                          MySolver * /*solver data*/,
                          MLAContext * /*mla context*/,
                          Vec * /*x*/,
                          Vec * /*b*/,
                          int /*pre-smooth times*/,
                          int /*post-smooth times*/,
                          int /*rbm order*/);

/*
 * mla solver solve phase recursive implementation
 *     linear system
 *     mla context, contains setup information
 *     config, pre- and post- smooth times
 *     a special case, rbm order is 2 and level is 1, coarse operator need shift
 *     level to recursive implementation
 */
void MLASolverSolvePhase(const ConfigJSON * /*config json*/,
                         MLAContext * /*mla context*/,
                         int /*rbm order*/,
                         MySolver * /*linear system data*/);

/*
 * mla solver setup phase
 *     number of levels
 *     rbm order
 *     neighbouring fine and coarse mesh
 *     aggregation mesh
 *     prolongation operator
 *     coarse operator
 */
void MLASovlerSetupPhase(MySolver * /*linear system data*/,
                         const MeshGraph * /*finest mesh data*/,
                         int /*number of levels*/,
                         int /*rbm order*/,
                         MLAContext * /*mla context*/);

/*
 * mla solver relative residual computing
 *     linear system
 *     relative residual
 */
void MLASolverRelativeResidual(MySolver * /*linear system data*/, double * /*relative residual*/);

/*
 * mla solver
 */
void MLASolver(const MeshGraph * /*graph data*/,
               MySolver * /*linear system data*/,
               const ConfigJSON * /*config data*/,
               int /*rbm order*/,
               MLAContext * /*mla context*/);

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
void MeshFileProcess(const char *path,
                     int label_bound, int label_omega,
                     GenericList *data_list_phy_tag,
                     GenericList *data_list_node,
                     GenericList *data_list_ele_bound,
                     GenericList *data_list_ele_omega);

void SolverPetscInitialize(const char *path_mat, const char *path_rhs, MySolver *mysolver);
void SolverPetscResidualCheck(MySolver *mysolver);

/*
 * json config parse
 */
int ConfigParse(const char * /*path to config file*/,
                ConfigJSON * /*config data*/);

#endif // main.h