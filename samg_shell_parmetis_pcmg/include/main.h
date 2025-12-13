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

#define MAX_SIZE PETSC_MAX_PATH_LEN
#define BUF_MAX_SIZE PETSC_MAX_PATH_LEN

// struct
/*
 * config file data
 */
typedef struct
{
    char file_mat[PETSC_MAX_PATH_LEN]; // base path to distributed matrix file
    char file_rhs[PETSC_MAX_PATH_LEN]; // base path to distributed vector file
    char file_vtx[PETSC_MAX_PATH_LEN]; // base path to distributed vertex file
    char file_adj[PETSC_MAX_PATH_LEN]; // base path to distributed adjacent list file
} CfgFile;

/*
 * config mg data
 */
typedef struct
{
    int pre_smooth, post_smooth; // pre/post-smooth times
    int num_level;               // number of levels
    int num_coarse_vtx;          // number of coarse level vertices
    int est_size_agg;            // estimation size of each aggregate
} CfgMG;

/*
 * config json data
 */
typedef struct
{
    CfgFile cfg_file; // config file data
    CfgMG cfg_mg;     // config mg data
} CfgJson;

/*
 * PETSc solver data
 */
typedef struct
{
    /* data */
    KSP ksp;
    PC pc;
    Mat solver_a;                     // A
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
} MySolver;

/*
 * CSR matrix data
 */
typedef struct
{
    int nrows, ncols, nnz;
    int *row_idx; // row indices, size: nrows
    int *row_ptr; // size: nrows + 1
    int *col_idx; // size: nnz
    double *val;  // size: nnz
} CSRMatrix;

/*
 * CSR vector data
 */
typedef struct
{
    int nrows;
    int *row_idx; // row indices, size: nrows
    double *val;  // size: nrows
} CSRVector;

/*
 * SAMG context data struct
 */
typedef struct
{
    /* data */
    CfgJson data_cfg; // config data
} SAMGCtx;

// function
/*
 * computing residual
 *     r = b - Ax
 */
int SolverPetscResidualCheck(MySolver *mysolver /*solver data*/);

/*
 * CSR vector file process
 */
int FileProcessCSRVector(const char *path_rhs /*path to rhs file*/,
                         CSRVector *data_rhs);

/*
 * CSR matrix file process
 */
int FileProcessCSRMatrix(const char *path_mat /*path to matrix file*/,
                         CSRMatrix *data_mat);

/*
 * solver initialize with file
 *     1. matrix
 *     2. rhs
 */
int SolverInitializeWithFile(const char *path_mat /*path to matrix file*/,
                             const char *path_rhs /*path to rhs file*/,
                             MySolver *mysolver /*solver data*/);

/*
 * parse config json file
 */
int ParseConfig(MPI_Comm comm /*communicator*/,
                const char *path /*path to json file*/,
                CfgJson *config /*json config data*/);

#endif // main.h