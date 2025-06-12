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
/*
 * gmsh struct
 */
typedef struct
{
    /* data */
    int nn, ne;               // number of nodes, elements
    int ne_bd, ne_in;         // number of elements of boundary, inner, ne = ne_bd + ne_in
    int nne_bd, nne_in;       // number of nodes in each element of boudary, inner
    idx_t nparts;             // number of partitions
    double *coordinates;      // coordinates of nodes [x1, y1, z1, x2, y2, z2, ...]
    idx_t *eptr_bd, *eind_bd; // csr mesh connectivity of boundary element
    idx_t *eptr_in, *eind_in; // csr mesh connectivity of inner element
    idx_t *npart_in;          // nodes partition of inner element
    idx_t *epart_in;          // elements partition of inner element
} DataGmsh;

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
    int coarse_restart;
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
typedef struct metis_mla_graph
{
    /* data */
    DataGmsh *fine;      // fine mesh
    DataGmsh *coarse;    // coarse mesh
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

typedef struct mla_context
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
    DataGmsh *data_gmsh;      // gmsh file data
} MLAContext;

// function prototype
/*
 * json config parse
 *     path (I) path to json file
 *     config (O) json config data
 */
int ConfigParse(MPI_Comm comm /*communicator*/,
                const char *path /*path to json file*/,
                ConfigJSON *config /*json config data*/);

#endif // main.h