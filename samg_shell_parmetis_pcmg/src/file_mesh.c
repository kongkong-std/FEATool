#include "../include/main.h"

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
            printf("path to matrix file: %s\n", path);
        }
    }
#endif // test path

    FILE *fp = fopen(path, "rb");
    assert(fp);

    fclose(fp);

    return 0;
}
