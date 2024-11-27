#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <cjson/cJSON.h>

// data struct
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
    MeshGraph *fine;                   // fine mesh
    MeshGraph *aggregation;            // aggregation mesh
    MeshGraph *coarse;                 // coarse mesh
    Mat prolongation, operator_coarse; // prolongation operator and coarse operator
    int level;                         // current level
} MLAGraph;

typedef struct mla_v_cycle
{
    /* data */
    int num_level; // total level
    MLAGraph *mla; // neighbouring level, size: num_level - 1
    int setup;     // setup flag, 0 represents un-setup, 1 represents setup
} MLAVCycle;

// function prototype
void SolverPetscInitialize(const char *path_mat, const char *path_rhs, MySolver *mysolver);
void SolverPetscResidualCheck(MySolver *mysolver);

/*
 * json config parse
 */
int ConfigParse(const char * /*path to config file*/,
                ConfigFile * /*config json of file*/,
                ConfigMLA * /*config json of aggregation-based multilevel*/);

#endif // main.h