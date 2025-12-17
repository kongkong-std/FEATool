#include "../include/main.h"

int FileProcessMeshAdj(const char *path_file /*path to mesh adjacent list file*/,
                       MeshAdj *data_adj /*mesh adjacent list data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    char path[PETSC_MAX_PATH_LEN];
    sprintf(path, "%s/adjacent_list-%d.txt", path_file, my_rank);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        if (index == my_rank)
        {
            printf(">>>> in rank %d:\n", my_rank);
            printf("path to adjacent list file: %s\n", path);
        }
    }
#endif // test path

    FILE *fp = fopen(path, "rb");
    assert(fp);

    int local_nv = 0;
    fscanf(fp, " %d ", &local_nv);

    data_adj->idx = (int *)malloc(local_nv * sizeof(int));
    data_adj->xadj = (int *)malloc((local_nv + 1) * sizeof(int));
    assert(data_adj->idx && data_adj->xadj);

    for (int index = 0; index < local_nv; ++index)
    {
        fscanf(fp, " %d ", data_adj->idx + index);
    }

    for (int index = 0; index < local_nv + 1; ++index)
    {
        fscanf(fp, " %d ", data_adj->xadj + index);
    }

    int local_len_adjncy = data_adj->xadj[local_nv] - data_adj->xadj[0];
    data_adj->adjncy = (int *)malloc(local_len_adjncy * sizeof(int));
    assert(data_adj->adjncy);

    for (int index = 0; index < local_len_adjncy; ++index)
    {
        fscanf(fp, " %d ", data_adj->adjncy + index);
    }

    fclose(fp);

    (void)MPI_Allreduce(&local_nv, &data_adj->nv,
                        1, MPI_INT, MPI_SUM, comm);

    return 0;
}

int FileProcessMeshVtx(const char *path_file /*path to mesh vertex file*/,
                       MeshVtx *data_vtx /*mesh vertex data*/)
{
    int my_rank, nprocs;
    MPI_Comm comm;
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    char path[PETSC_MAX_PATH_LEN];
    sprintf(path, "%s/vertex_coor-%d.txt", path_file, my_rank);
#if 0
    for (int index = 0; index < nprocs; ++index)
    {
        if (index == my_rank)
        {
            printf(">>>> in rank %d:\n", my_rank);
            printf("path to vertex file: %s\n", path);
        }
    }
#endif // test path

    FILE *fp = fopen(path, "rb");
    assert(fp);

    int local_nv = 0;
    fscanf(fp, " %d ", &local_nv);
    data_vtx->local_nv = local_nv;

    data_vtx->idx = (int *)malloc(local_nv * sizeof(int));
    data_vtx->type = (int *)malloc(local_nv * sizeof(int));
    data_vtx->data_coor = (CoorData *)malloc(local_nv * sizeof(CoorData));
    assert(data_vtx->idx && data_vtx->type && data_vtx->data_coor);

    for (int index = 0; index < local_nv; ++index)
    {
        fscanf(fp, " %d %d %lf %lf %lf %lf %lf %lf ", data_vtx->idx + index,
               data_vtx->type + index,
               &data_vtx->data_coor[index].x,
               &data_vtx->data_coor[index].y,
               &data_vtx->data_coor[index].z,
               &data_vtx->data_coor[index].nx,
               &data_vtx->data_coor[index].ny,
               &data_vtx->data_coor[index].nz);
    }

    fclose(fp);

    (void)MPI_Allreduce(&local_nv, &data_vtx->nv,
                        1, MPI_INT, MPI_SUM, comm);

    return 0;
}
