#include "../include/main.h"

int main(int argc, char **argv)
{
    char *path_msh = NULL;
    char *dst_phy_tag = NULL, *dst_node = NULL;
    char *dst_ele_bound = NULL, *dst_ele_omega = NULL;
    int label_bound = 0, label_omega = 0;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_msh = argv[index + 1];
        }
        if (strstr("-physics_tag", argv[index]))
        {
            dst_phy_tag = argv[index + 1];
        }
        if (strstr("-node", argv[index]))
        {
            dst_node = argv[index + 1];
        }
        if (strstr("-element_bound", argv[index]))
        {
            dst_ele_bound = argv[index + 1];
        }
        if (strstr("-element_omega", argv[index]))
        {
            dst_ele_omega = argv[index + 1];
        }
        if (strstr("-label_bound", argv[index]))
        {
            label_bound = atoi(argv[index + 1]);
        }
        if (strstr("-label_omega", argv[index]))
        {
            label_omega = atoi(argv[index + 1]);
        }
    }

    GenericList data_list_phy_tag, data_list_node;
    GenericList data_list_ele_bound, data_list_ele_omega;
    InitializeList(&data_list_phy_tag);
    InitializeList(&data_list_node);
    InitializeList(&data_list_ele_bound);
    InitializeList(&data_list_ele_omega);

    MeshFileProcess(path_msh, label_bound, label_omega,
                    &data_list_phy_tag, &data_list_node,
                    &data_list_ele_bound, &data_list_ele_omega);

    // print data linked list
    printf("==== data linked list ====\n");
    printf("size physical tag: %d\n", data_list_phy_tag.size);
    printf("size node: %d\n", data_list_node.size);
    printf("size bound element: %d\n", data_list_ele_bound.size);
    printf("size omega element: %d\n", data_list_ele_omega.size);

#if 0
    PrintList(&data_list_phy_tag);
    PrintList(&data_list_node);
    PrintList(&data_list_ele_bound);
    PrintList(&data_list_ele_omega);
#endif // print data linked list

#if 1
    PrintListFile(dst_phy_tag, &data_list_phy_tag);
    PrintListFile(dst_node, &data_list_node);
    PrintListFile(dst_ele_bound, &data_list_ele_bound);
    PrintListFile(dst_ele_omega, &data_list_ele_omega);
#endif // output data linked list to file

    // free memory
    ClearList(&data_list_phy_tag);
    ClearList(&data_list_node);
    ClearList(&data_list_ele_bound);
    ClearList(&data_list_ele_omega);

    return 0;
}

// command line
/*
 * ./app_exe -mesh ../input/linear-elasticity-2d-22.msh
 *     -physics_tag physics-tag.txt
 *     -node node-data.txt
 *     -element_bound element-data-bound.txt
 *     -element_omega element-data-omega.txt
 *     -label_bound 1 -label_omega 2
 */
