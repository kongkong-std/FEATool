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
#include <limits.h>
#include <lapacke.h>

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
 * coordinate data
 */
typedef struct
{
    /* data */
    double x, y, z;
    double nx, ny, nz;
} CoorData;

/*
 * mesh vertex information data
 */
typedef struct
{
    int nv; // global number of vertices
    int local_nv;
    int *idx;            // size: local_nv, local vertex index
    int *type;           // size: local_nv, 0: shell, 1: solid
    CoorData *data_coor; // size: local_nv, vertex coordinate and normal data
} MeshVtx;

/*
 * mesh adjacent list data
 */
typedef struct
{
    /* data */
    int nv; // global number of vertices
    int local_nv;
    int *idx;    // size: local_nv, local vertex index
    int *xadj;   // size: local_nv + 1, csr row ptr
    int *adjncy; // size: xadj[local_nv] - xadj[0], adjacent vertices
} MeshAdj;

/*
 * mesh data
 */
typedef struct
{
    /* data */
    int nv;           // number of vertices in mesh
    int np;           // number of partitions
    int *vtxdist;     // size: np + 1, vertex distribution
    MeshVtx data_vtx; // mesh vertex data
    MeshAdj data_adj; // mesh adjacent list data
    int *parts;       // size: local_nv, partition
} MeshData;

/*
 * near null space data of each vertex in level 0
 * LAPACK_COL_MAJOR remember!!!
 */
typedef struct
{
    /* data */
    int idx;        // vertex id
    int type;       // vertex type, 0: shell, 1: solid
    int nrow, ncol; // for shell, nrow = ncol = 6; for solid, nrow = 3, ncol = 6
    double val[36]; // size: for shell, 6 \times 6; for solid, 3 \times 6
} NearNullSpaceDataVertexLevel0;

/*
 * near null space data of each vertex in level k
 */
typedef struct
{
    /* data */
    int idx;
    int nrow, ncol; // in level k (except level 0), nrow = ncol = 6
    double val[36]; // size: nrow x ncol
} NearNullSpaceDataVertexLevelK;

/*
 * Q matrix data of each vertex in level k
 */
typedef struct
{
    /* data */
    int idx;
    int nrow, ncol; // same with near null space data
    double val[36]; // size: nrow x ncol
} QDataVertexLevelK;

/*
 * near null space of level 0
 */
typedef struct
{
    int size_global;                               // global number of vertices in file-level mesh
    int size_local;                                // local number of vertices in fine-level mesh
    NearNullSpaceDataVertexLevel0 *data_nullspace; // size: local number of vertices in fine-level mesh
} NearNullSpaceLevel0;

/*
 * near null space of level k, k > 0
 */
typedef struct
{
    /* data */
    int size_global;                               // global number of vertices in fine-level mesh
    int size_local;                                // local number of vertices in fine-level mesh
    NearNullSpaceDataVertexLevelK *data_nullspace; // size: local number of vertices in fine-level mesh
} NearNullSpaceLevelK;

/*
 * Q matrix from near null space data
 * corresponding to dof order or vertex id
 */
typedef struct
{
    /* data */
    int size_global;           // global number of vertices in fine-level mesh
    int size_local;            // local number of vertices in fine-level mesh
    QDataVertexLevelK *data_q; // size: local number of vertices in fine-level mesh
} QLevelK;

/*
 * ghost mapping on owner rank
 */
typedef struct
{
    int *flid2fgid; // owner-local fine-id mapping to global fine-id
    // int *fgid2flid; // global fine-id mapping to owner-local fine-id
    int n_owned; // number of owned fine vertices in this partition
    int nghost;  // number of ghost fine vertices in this partition
    int n_total; // = n_owned + nghost

    /*
     * mat_t = mat_q * mat_r
     *     size(mat_t) = m \times n
     *     size(mat_q) = m \times n
     *     size(mat_r) = n \times n
     * after calling OpenBLAS QR factorization API,
     *     LAPACKE_dgeqrf() implement QR factorization and mat_r in upper triangular part of mat_t
     *     LAPACKE_dorgqr() store mat_q explicitly in mat_t
     */
    int nrow, ncol; // size(mat_t) = nrow x ncol
    double *mat_t;  // copying mat_t to mat_q
    double *mat_q;  // after calling OpenBLAS, size: nrow x ncol
    double *mat_r;  // size: ncol x ncol
} GhostAggData;

/*
 * aggregation data
 */
typedef struct
{
    int np;     // number of partitions (global)
    int *owner; // size: np, owner rank of each partition

    /* local appearance */
    int *nlocal;      // size: np
    int **local_gids; // local fine global vertex IDs per partition

    int *nlocal_all; // size: nprocs * np. Stores vertex counts contributed by each rank to each partition.
    // char padding[1024];

    /* owner-only */
    int *n_fine;                                           // size: np
    int **fine_gids;                                       // full fine vertex list per partition
    NearNullSpaceDataVertexLevelK **fine_global_nullspace; // global near null space, global fine null space per partition

    /* mapping from fine-level vertex to partition */
    int nv;         // number of vertices in fine-level
    int *fgid2part; // size: nv

    /* ghost mapping on owner */
    GhostAggData *data_ghost_agg; // size: np
} AggData;

/*
 * multigrid level hierarchy data
 */
typedef struct
{
    MeshData data_f_mesh, data_c_mesh;         // fine-/coarse- level mesh data
    AggData data_agg;                          // aggregation data
    NearNullSpaceLevelK data_nullspace_levelk; // level k near null space data, in level k except level 0
    QLevelK data_q_levelk;                     // Q matrix data from near null space data, size same with near null space data
    Mat op_f, op_c;                            // fine-/coarse- level operators
    Mat op_ua_p;                               // unsmoothed aggregation prolongation operator
    Mat op_sa_p;                               // smoothed aggregation prolongation operator
    PSmoother op_s;                            // prolongation operator smoother
} MGLevel;

/*
 * SAMG context data struct
 */
typedef struct
{
    /* data */
    CfgJson data_cfg;                          // config data
    int num_level;                             // true number of levels
    MySolver mysolver;                         // setup phase, mysolver data
    NearNullSpaceLevel0 data_nullspace_level0; // level 0 near null space data, only in level 0
    MGLevel *levels;                           // level hierarchy information
} SAMGCtx;

// function
/*
 * free memory of samg level operator data
 */
int FreeSAMGMatData(MGLevel *level /*samg level data*/);

/*
 * free memory of Q matrix data
 */
int FreeQLevelK(QLevelK *q /*Q matrix data*/);

/*
 * free memory of aggregation data
 *     1. ghost data
 *     2. aggregation data
 */
int FreeGhostAggData(GhostAggData *g /*ghost data of aggregation data*/);
int FreeAggData(AggData *agg /*aggregation data*/, int my_rank /*current rank*/);

/*
 * free memory of near null space data
 *     1. level 0 near null space data
 *     2. level k near null space data
 */
int FreeNearNullSpaceLevel0(NearNullSpaceLevel0 *ns /*near null space level 0*/);
int FreeNearNullSpaceLevelK(NearNullSpaceLevelK *ns /*near null space level k*/);

/*
 * free memory of mesh data
 *     1. mesh adjacency data
 *     2. mesh vertex data
 */
int FreeMeshAdj(MeshAdj *data_adj /*mesh adjacency data*/);
int FreeMeshVtx(MeshVtx *data_vtx /*mesh vertex data*/);
int FreeMeshData(MeshData *mesh /*mesh data*/);

/*
 * PCMG setup from SAMG
 */
int PCMGSetupFromSAMG(int sa_flag /*flag of sa*/,
                      SAMGCtx *samg_ctx /*samg context data*/,
                      MySolver *mysolver /*solver data*/);

/*
 * level k mesh data generator
 */
int SAMGLevelKMesh(const CfgJson *data_cfg /*config data*/,
                   MeshData *data_f_mesh /*fine-level mesh data*/,
                   MeshData *data_c_mesh /*coarse-level mesh data*/,
                   AggData *data_agg /*aggregation data*/);

/*
 * copying coarse-level mesh data to next level fine-level mesh data
 */
int CoarseMesh2FineMesh(const MeshData *data_c_mesh /*coarse-level mesh*/,
                        MeshData *data_f_mesh /*next fine-level mesh*/);

/*
 * calling QR factorization from OpenBLAS
 */
void QRFactorizationOpenBLAS(int m /*nrow*/, int n /*ncol*/,
                             double *Q /*m x n*/, double *R /*n x n*/);

/*
 * smoothed aggregation coarse operator constructor
 */
int SAMGSACoarseOperator(SAMGCtx **samg_ctx /*samg context data*/);

/*
 * unsmoothed aggregation coarse operator constructor
 */
// int SAMGUACoarseOperator(SAMGCtx **samg_ctx /*samg context data*/);

/*
 * level k Q matrix
 */
int SAMGLevelKQMatrix(MeshData *data_mesh_f /*fine-level mesh data*/,
                      NearNullSpaceLevelK *data_nullspace_levelk /*fine-level near null space data*/,
                      AggData *data_agg /*aggregation data*/,
                      QLevelK *data_q_levelk /*level k Q matrix*/);

/*
 * Q matrix
 * size and parallel distribution same with near null space of fine-level
 */
int SAMGLevelQMatrix(SAMGCtx **samg_ctx /*samg context data*/);

/*
 * level k tentative prolongation operator constructor
 */
int SAMGLevelKTentativeProlongationOperator(int level /*current level*/,
                                            MGLevel *data_level /*level data*/);

/*
 * tentative prolongation operator constructor
 */
int SAMGTentativeProlongationOperator(SAMGCtx **samg_ctx /*samg context data*/);

/*
 * near null space ghost data of aggregation data
 */
int SAMGLevelKGhostDataNearNullSpace(const int *vtxdist_f /*fine-lelve mesh vtxdist array*/,
                                     NearNullSpaceLevelK *data_nullspace_f /*fine-level near null space*/,
                                     AggData *data_agg /*aggregation data*/);

/*
 * level k near null space
 */
int SAMGLevelKNearNullSpace(MeshData *data_mesh_f /*fine-level mesh data*/,
                            MeshData *data_mesh_c /*coarse-level mesh data*/,
                            NearNullSpaceLevelK *data_nullspace_f /*fine-level near null space*/,
                            AggData *data_agg /*aggregation data*/,
                            NearNullSpaceLevelK *data_nullspace_c /*coarse-level near null space*/);

/*
 * level 0 near null space
 */
int SAMGLevel0NearNullSpace(NearNullSpaceLevel0 *data_nullspace_level0 /*initial near null space*/,
                            NearNullSpaceLevelK *data_nullspace_levelk /*level 0 near null space*/);

/*
 * initial near null space depends on level 0 fine mesh
 */
int SAMGInitialNearNullSpace(MeshData *data_f_mesh /*fine level 0 mesh data*/,
                             NearNullSpaceLevel0 *data_nullspace_level0 /*level 0 near null space*/);

/*
 * construct multilevel near null space
 */
int SAMGLevelNearNullSpace(SAMGCtx **samg_ctx /*samg context data*/);

/*
 * mapping from fine-level vertex to partition id
 */
int SAMGFineVertex2PartitionMap(AggData **agg /*aggregation data*/);

/*
 * construct coarse-level adjacent list
 */
int SAMGCoarseAdjacentListConstructor(AggData **agg /*aggregation data*/,
                                      MeshData **mesh_f /*fine-level mesh data*/,
                                      MeshData **mesh_c /*coarse-level mesh data*/);

/*
 * coarse-level vertex coordinate data
 */
int SAMGCoarseVertexCoordinate(AggData **agg /*aggregation data*/,
                               MeshData **mesh_f /*fine-level mesh data*/,
                               MeshData **mesh_c /*coarse-level mesh data*/);

/*
 * construct partition ghost data mapping
 */
int SAMGPartitionGhostDataMapping(AggData **agg /*aggregation data*/);

/*
 * renumbering partition_id generated by ParMetis
 */
int SAMGPartitionRenumberID(AggData **agg /*aggregation data*/,
                            MeshData **mesh_f /*fine-level mesh data*/,
                            MeshData **mesh_c /*coarse-level mesh data*/);

/*
 * partition owner vertex data
 */
int SAMGPartitionOwnerVertexData(AggData **agg /*aggregation data*/,
                                 MeshData **mesh_f /*fine-level mesh data*/);

/*
 * partition owner rank constructor
 */
int SAMGPartitionOwnerRankConstructor(AggData **agg /*aggregation data*/,
                                      MeshData **mesh_f /*fine-level mesh data*/);

/*
 * coarse-level mesh constructor
 */
int SAMGCoarseMeshConstructor(MeshData **mesh_f /*fine-level mesh data*/,
                              MeshData **mesh_c /*coarse-level mesh data*/,
                              AggData **agg /*aggregation data*/);

/*
 * mesh adjacent list process
 */
int FileProcessMeshAdj(const char *path_file /*path to mesh adjacent list file*/,
                       MeshAdj *data_adj /*mesh adjacent list data*/);

/*
 * mesh vertex file process
 */
int FileProcessMeshVtx(const char *path_file /*path to mesh vertex file*/,
                       MeshVtx *data_vtx /*mesh vertex data*/);

/*
 * level 0 mesh generator
 */
int SAMGLevel0Mesh(const CfgJson *data_cfg /*config data*/,
                   MeshData *data_f_mesh /*fine-level mesh data*/,
                   MeshData *data_c_mesh /*coarse-level mesh data*/,
                   AggData *data_agg /*aggregation data*/);

/*
 * mesh hierarchy generator
 */
int SAMGLevelMesh(int cfg_mg_num_level /*config number of levels*/,
                  SAMGCtx **samg_ctx /*samg context data*/);

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
int SAMGSetupPhase(SAMGCtx *samg_ctx /*samg context data*/,
                   int sa_flag /*flag of sa*/);

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
