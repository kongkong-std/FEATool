#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <cjson/cJSON.h>
#include <stdbool.h>
#include <parmetis.h>
#include <metis.h>
// #include <GKlib.h>

#define MAX_SIZE PETSC_MAX_PATH_LEN
#define BUF_MAX_SIZE PETSC_MAX_PATH_LEN

// data struct
typedef struct
{
    /* data */
    int ele_type;
    int num_ele_node; // number of nodes in current element
    int *ele_node;    // nodes list
} DataMeshEle;

typedef struct
{
    /* data */
    int nn, ne;                      // number of elements, nodes
    int dim;                         // dimension of coordinates, (1, 2 or 3)
    double *coordinates;             // coordinates of nodes
    int ne_solid, ne_shell, ne_beam; // number of solid, shell, beam element
    DataMeshEle *ele_solid;          // solid element data
    DataMeshEle *ele_shell;          // shell element data
    DataMeshEle *ele_beam;           // beam element data
} DataMesh;

typedef struct
{
    /* data */
    int nn;              // number of nodes
    int dim;             // dimensions
    double *coordinates; // coordinates of nodes
    idx_t *vtxdist;      // parmetis vtxdist parameter, global node indicies
    idx_t *xadj;         // csr row pointer
    idx_t *adjncy;       // adjacency nodes list
    idx_t nparts;        // number of super nodes (partitions)
    idx_t *part;         // partition value
    idx_t *local_xadj;   // local xadj
    idx_t *local_adjncy; // local adjncy
    idx_t *local_part;   // local part
} AdjDataMesh;

/*
 * distributed csr matrix
 */
typedef struct
{
    /* data */
    int nrows, ncols; // nrows, ncols
    int nnz;
    int *row_ptr;
    int *col_idx;
    double *val;

    int local_nrows;
    int local_nnz;
    int *local_row_ptr;
    int *local_col_idx;
    double *local_val;
} CSRMatrix;

/*
 * distributed rhs vector
 */
typedef struct
{
    /* data */
    int nrows;
    double *val;

    int local_nrows;
    double *local_val;
} CSRVector;

/*
 * mesh data block flag
 */
typedef enum
{
    NONE,           // 0
    MESH_FORMAT,    // 1
    PHYSICAL_NAMES, // 2
    NODES,          // 3
    ELEMENTS        // 4
} Flag_Data_Block;

typedef Flag_Data_Block FlagDataBlockGmsh;

/*
 * json data struct
 */
typedef struct
{
    char file_mat[PETSC_MAX_PATH_LEN];
    char file_rhs[PETSC_MAX_PATH_LEN];
    char file_mesh[PETSC_MAX_PATH_LEN];
} ConfigFile;

typedef struct
{
    int pre_smooth_v;
    int post_smooth_v;
    int mla_max_it;
    double mla_rtol;
    int mla_level;
    int mla_phase;
    int coarse_restart;
} ConfigMLA;

typedef struct
{
    /* data */
    int label_bound, label_omega;
} ConfigMeshLabel;

typedef struct
{
    /* data */
    ConfigFile file_config;
    ConfigMLA mla_config;
    ConfigMeshLabel mesh_label_config;
} ConfigJSON;

/*
 * solver data struct
 */
typedef struct
{
    KSP ksp;
    PC pc;
    Mat solver_a;
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
} MySolver;

/*
 * mla data struct
 */
typedef struct
{
    /* data */
    AdjDataMesh *fine;   // fine mesh
    AdjDataMesh *coarse; // coarse mesh
    Mat prolongation;    // prolongation operator
    Mat operator_coarse; // coarse operator
    Mat operator_fine;   // fine operator
    KSP ksp_presmooth;   // pre-smooth solver
    KSP ksp_postsmooth;  // post-smooth solver
    KSP ksp_coarse;      // coarse solver
    PC pc_presmooth;     // pre-smooth preconditioner
    PC pc_postsmooth;    // post-smooth preconditioner
    PC pc_coarse;        // coarse preconditioner
    int level;           // current level
} MetisMLAGraph;

typedef struct
{
    /* data */
    int num_level;            // total level
    int true_num_level;       // true level counts
    int setup;                // setup flag, 0 represents un-setup, 1 represents setup
    char *path_config;        // path to config file
    ConfigJSON config;        // json config
    int order_rbm;            // rbm order
    int angle_type;           // 0: theta (rotationan angle with axis), 1: phi (normal rotational angle with axis)
    MySolver mysolver;        // mysolver data
    MetisMLAGraph *metis_mla; // metis mla data
    AdjDataMesh *data_mesh;   // gmsh file data
} MLAContext;

// function prototype
/*
 * copy current level coarse graph to next level fine graph
 *     coarse_graph (I) current level coarse graph
 *     fine_graph (O) next level fine graph
 */
int DeepCopyCoarse2NextLevelFine(const AdjDataMesh *coarse_graph /*coarse level graph*/,
                                 AdjDataMesh *fine_graph /*next level fine graph*/);

/*
 * level k coarse mesh partition, from level k fine mesh
 *     fine_graph (I) level k fine graph
 *     coarse_graph (O) level k coarse graph
 */
int LevelKCoarsePartition(const AdjDataMesh *fine_graph /*level k fine data mesh*/,
                          AdjDataMesh *coarse_graph /*level k coarse data mesh*/);

/*
 * level 0 fine mesh partition, from initial mesh data
 *     data_mesh (I) initial data mesh
 *     fine_graph (O) level 0 fine graph data
 */
int Level0FinePartition(const AdjDataMesh *data_mesh /*initial data mesh*/,
                        AdjDataMesh *fine_graph /*level 0 fine graph data*/);

/*
 * coarse level graph data generated from fine level graph data
 *     fine_graph_data (I) fine level graph data
 *     coarse_graph_data (O) coarse level graph data
 */
int CoarseLevelGenerator(const AdjDataMesh *fine_graph_data /*fine level graph data*/,
                         AdjDataMesh *coarse_graph_data /*coarse level graph data*/);

/*
 * pcshell function solve
 *     pc (I) matrix
 *     x_in (I) rhs
 *     x_out (O) solution
 */
extern PetscErrorCode ParMetisMLAShellPCApply(PC pc, Vec x_in, Vec x_out);

/*
 * pcshell function setup
 *     pc (O) setup
 */
extern PetscErrorCode ParMetisMLAShellPCSetup(PC pc);

/*
 * parmetis mla solver, setup phase
 *     mla_ctx (IO) setup
 */
int ParMetisMLASolverSetupPhase(MLAContext *mla_ctx /*mla context data*/);

/*
 * parmetis mla solver
 *     mla_ctx (IO) setup, solve
 *     mla_pahse (I) 0: setup phase, 1: solve phase, 2 : setup + solve phase
 */
int ParMetisMLASolver(MLAContext *mla_ctx /*mla context data*/,
                      int mla_phase /*multilevel phase*/);

/*
 * csr graph data generator
 *     mesh_data (I) mesh data
 *     graph_data (O) csr type graph data
 */
int GlobalGraphCSRAdjGenerator(const DataMesh *mesh_data /*mesh data*/,
                               AdjDataMesh **graph_data /*csr graph data*/);

/*
 * mesh data file
 *     1. gmsh file
 *     2. comsol mesh file
 */
int FileProcessMesh(const char *path /*path to mesh file*/, DataMesh *mesh_data /*mesh data*/);

/*
 * computing residual
 *     mysolver (I) linear solver data
 */
int SolverPetscResidualCheck(MySolver *mysolver /*solver data*/);

/*
 * assembling petsc linear system
 *     path_mat (I) path to petsc binary matrix file
 *     path_rhs (I) path to petsc binary rhs file
 *     node_vtxdist (I) distribution of nodes
 *     mysolver (O) linear solver data
 */
int SolverPetscInitialize(const char *path_mat /*path to matrix file*/,
                          const char *path_rhs /*path to rhs file*/,
                          const int *node_vtxdist /*distribution of nodes*/
                          MySolver *mysolver /*solver data*/);

/*
 * json config parse
 *     path (I) path to json file
 *     config (O) json config data
 */
int ConfigParse(MPI_Comm comm /*communicator*/,
                const char *path /*path to json file*/,
                ConfigJSON *config /*json config data*/);

#endif // main.h