#include "../include/featool_main.h"

int main(int argc, char **argv)
{
    printf("==== test FEATool ====\n");

    int mesh_dim = 0;
    char *path_mesh_element = NULL, *path_mesh_node = NULL;

    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-dim", argv[index]))
        {
            // dimension of mesh
            mesh_dim = atoi(argv[index + 1]);
        }
        if(strstr("-mesh_element", argv[index]))
        {
            // path to mesh element file
            path_mesh_element = argv[index + 1];
        }
        if(strstr("-mesh_node", argv[index]))
        {
            // path to mesh node file
            path_mesh_node = argv[index + 1];
        }
    }

    printf("dimension of mesh: %d\n", mesh_dim);

    return 0;
}