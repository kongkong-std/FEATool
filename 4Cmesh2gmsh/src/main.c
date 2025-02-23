#include "../include/main.h"

int main(int argc, char **argv)
{
    // puts("testtest...");
    char *path_dat = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-4c_dat", argv[index]))
        {
            path_dat = argv[index + 1];
        }
    }

    FourCDatNodeData *node_data;
    FourCDatEleData *ele_data;

    FourCDatFileProcess(path_dat, &node_data, &ele_data);

    // free memory
    free(node_data);
    free(ele_data);

    return 0;
}

// command line
/*
 * -4c_dat <path/to/4c/dat/file>
 */
