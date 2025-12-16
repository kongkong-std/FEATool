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
    int ps_num_steps;            // number of smoothed prolongation operators steps
    int ps_type;                 // prolongation operator smoother type
    double ps_scale;             // scaled prolongation operator smoother, default is 0.67
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
 * prolongation operator smoother
 *     P_sa = (I - P_smoother A)^{v} P_us
 */
typedef struct
{
    int num_steps;         // number of smoothing steps
    int smoother_type;     // 0: scaled weighted jacobi, 1: GS
    double smoother_scale; // scaled smoother, default 0.67
    PC pc;                 // smoother, PCJACOBI, ...
} PSmoother;

/*
 * multigrid level hierarchy data
 */
typedef struct
{
    Mat op_f, op_c; // fine-/coarse- level operators
    Mat op_ua_p;    // unsmoothed aggregation prolongation operator
    Mat op_sa_p;    // smoothed aggregation prolongation operator
    PSmoother op_s; // prolongation operator smoother
} MGLevel;

/*
 * SAMG context data struct
 */
typedef struct
{
    /* data */
    CfgJson data_cfg;  // config data
    int num_level;     // true number of levels
    MySolver mysolver; // setup phase, mysolver data
    MGLevel *levels;   // level hierarchy information
} SAMGCtx;

// function
/*
 * smoothed prolongation
 *     P_sa = (I - S A)^\nu P_ua
 */
int SAMGSmoothedProlongation(MGLevel *level /*level hierarchy data*/);

/*
 * apply one-time prolongation smoother
 *     P_sa = (I - S A) P_ua
 */
int SAMGApplyProlongationSmoother(int n /*column size of prolongation operator*/,
                                  double omega /*scaling weight parameter*/,
                                  PSmoother *p_s /*prolongation operator smoother*/,
                                  Mat *p_sa /*smoothed prolongation operator*/,
                                  Mat *p_ua /*unsmoothed prolongation operator*/,
                                  Mat *A /*fine-level operator*/);

/*
 * samg setup phase
 */
int SAMGSetupPhase(SAMGCtx *samg_ctx /*samg context data*/);

/*
 * deep copy
 */
int DeepCopyMySolverData(MySolver *dst /*destinated mysolver data*/,
                         MySolver *src /*source mysolver data*/);

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
