#include "../include/main.h"

int main(int argc, char **argv)
{
    // puts("testtest...");
    char *path_dat = NULL, *path_msh = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-4c_dat", argv[index]))
        {
            path_dat = argv[index + 1];
        }
        if (strstr("-msh_file", argv[index]))
        {
            path_msh = argv[index + 1];
        }
    }

    FourCDatNodeData *node_data;
    FourCDatEleData *ele_data;
    FourCDatEdgeEleData *edge_ele_data;
    int num_node = 0, num_ele = 0, num_edge_ele = 0;

    FourCDatFileProcess(path_dat, &num_node, &node_data, &num_ele, &ele_data, &num_edge_ele, &edge_ele_data);
#if 0
    for (int index = 0; index < num_node; ++index)
    {
        printf("%d\t%d\t%021.16le\t%021.16le\t%021.16le\n", (node_data + index)->omega_flag,
               (node_data + index)->idx,
               (node_data + index)->x,
               (node_data + index)->y,
               (node_data + index)->z);
    }
    for (int index = 0; index < num_ele; ++index)
    {
        printf("%d\t%d\t", (ele_data + index)->idx, (ele_data + index)->msh_type);
        for (int index_ele = 0; index_ele < (ele_data + index)->num_cell; ++index_ele)
        {
            printf("%d\t", (ele_data + index)->node_list[index_ele]);
        }
        putchar('\n');
    }
    for (int index = 0; index < num_edge_ele; ++index)
    {
        printf("%d\t%d\t", (edge_ele_data + index)->idx, (edge_ele_data + index)->msh_type);
        for (int index_edge_ele = 0; index_edge_ele < (edge_ele_data + index)->num_cell; ++index_edge_ele)
        {
            printf("%d\t", (edge_ele_data + index)->node_list[index_edge_ele]);
        }
        putchar('\n');
    }
#endif

    FourCDat2MSHFile(path_msh, num_node, node_data,
                     num_ele, ele_data,
                     num_edge_ele, edge_ele_data);

    // puts("hhhhhhhh tettttt");

    // free memory
    free(edge_ele_data);
    free(node_data);
    for (int index = 0; index < num_ele; ++index)
    {
        free((ele_data + index)->node_list);
    }
    free(ele_data);

    return 0;
}

// command line
/*
 * -4c_dat <path/to/4c/dat/file>
 * -msh_file <path/to/msh/file>
 */
