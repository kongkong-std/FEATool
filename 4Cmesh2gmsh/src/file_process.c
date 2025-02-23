#include "../include/main.h"

void FourCDatFileProcess(const char *path, FourCDatNodeData **node_data, FourCDatEleData **ele_data)
{
    // printf("path to dat file: %s\n", path);
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[MAX_SIZE];
    int num_ele = 0, num_node = 0;
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "PROBLEM SIZE"))
        {
            fscanf(fp, "%*s%d", &num_ele);
            fscanf(fp, "%*s%d", &num_node);
            break;
        }
    }
#if 0
    printf("num_ele = %d, num_node = %d\n", num_ele, num_node);
#endif

    // FourCDatNodeData *node_data = NULL;
    // FourCDatEleData *ele_data = NULL;
    *node_data = (FourCDatNodeData *)malloc(num_node * sizeof(FourCDatNodeData));
    *ele_data = (FourCDatEleData *)malloc(num_ele * sizeof(FourCDatEleData));
    assert(*node_data && *ele_data);

    rewind(fp);

    // node data
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "NODE COORDS"))
        {
            for (int index = 0; index < num_node; ++index)
            {
                fscanf(fp, "%*s%d%*s%lf%lf%lf", &((*node_data + index)->idx),
                       &((*node_data + index)->x),
                       &((*node_data + index)->y),
                       &((*node_data + index)->z));
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < num_node; ++index)
    {
        printf("%d: %021.16le\t%021.16le\t%021.16le\n", (*node_data)[index].idx,
               (*node_data)[index].x,
               (*node_data)[index].y,
               (*node_data)[index].z);
    }
#endif

    rewind(fp);

    // element data
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "STRUCTURE ELEMENTS"))
        {
            // 1 WALL QUAD4 2 4 3 1 MAT 1 KINEM nonlinear EAS none THICK 1.0 STRESS_STRAIN plane_stress GP 2 2
            for (int index = 0; index < num_ele; ++index)
            {
                fgets(buffer, MAX_SIZE, fp);
                sscanf(buffer, "%d%*s%s", &((*ele_data + index)->idx),
                       (*ele_data + index)->cell_type);
                int num_node_cell_type = FourCNumNodeCellType((*ele_data)[index].cell_type, &((*ele_data + index)->msh_type));
                // printf(">>>> num node = %d, msh type = %d >>>>\n", num_node_cell_type, (*ele_data)[index].msh_type);
                (*ele_data + index)->node_list = (int *)malloc(num_node_cell_type * sizeof(int));
                assert((*ele_data + index)->node_list);
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < num_ele; ++index)
    {
        printf("%d: %s\n", (*ele_data)[index].idx, (*ele_data)[index].cell_type);
    }
#endif

    fclose(fp);
}

int FourCNumNodeCellType(const char *cell_type, int *msh_type)
{
    int value = 0;

    if (strcmp(cell_type, "POINT1") == 0)
    {
        value = 1;
    }
    else if (strcmp(cell_type, "LINE2") == 0)
    {
        value = 2;
        *msh_type = 1;
    }
    else if (strcmp(cell_type, "LINE3") == 0)
    {
        value = 3;
    }
    else if (strcmp(cell_type, "LINE4") == 0)
    {
        value = 4;
    }
    else if (strcmp(cell_type, "LINE5") == 0)
    {
        value = 5;
    }
    else if (strcmp(cell_type, "LINE6") == 0)
    {
        value = 6;
    }
    else if (strcmp(cell_type, "NURBS2") == 0)
    {
        value = 2;
    }
    else if (strcmp(cell_type, "NURBS3") == 0)
    {
        value = 3;
    }
    else if (strcmp(cell_type, "QUAD4") == 0)
    {
        value = 4;
        *msh_type = 3;
    }
    else if (strcmp(cell_type, "QUAD6") == 0)
    {
        value = 6;
    }
    else if (strcmp(cell_type, "QUAD8") == 0)
    {
        value = 8;
    }
    else if (strcmp(cell_type, "QUAD9") == 0)
    {
        value = 9;
    }
    else if (strcmp(cell_type, "TRI3") == 0)
    {
        value = 3;
        *msh_type = 2;
    }
    else if (strcmp(cell_type, "TRI6") == 0)
    {
        value = 6;
    }
    else if (strcmp(cell_type, "NURBS4") == 0)
    {
        value = 4;
    }
    else if (strcmp(cell_type, "NURBS9") == 0)
    {
        value = 9;
    }
    else if (strcmp(cell_type, "HEX8") == 0)
    {
        value = 8;
        *msh_type = 5;
    }
    else if (strcmp(cell_type, "HEX16") == 0)
    {
        value = 16;
    }
    else if (strcmp(cell_type, "HEX18") == 0)
    {
        value = 18;
    }
    else if (strcmp(cell_type, "HEX20") == 0)
    {
        value = 20;
    }
    else if (strcmp(cell_type, "HEX27") == 0)
    {
        value = 27;
    }
    else if (strcmp(cell_type, "TET4") == 0)
    {
        value = 4;
        *msh_type = 4;
    }
    else if (strcmp(cell_type, "TET10") == 0)
    {
        value = 10;
    }
    else if (strcmp(cell_type, "WEDGE6") == 0)
    {
        value = 6;
    }
    else if (strcmp(cell_type, "WEDGE15") == 0)
    {
        value = 15;
    }
    else if (strcmp(cell_type, "PYRAMID5") == 0)
    {
        value = 5;
        *msh_type = 7;
    }
    else if (strcmp(cell_type, "NURBS8") == 0)
    {
        value = 8;
    }
    else if (strcmp(cell_type, "NURBS27") == 0)
    {
        value = 27;
    }

    return value;
}
