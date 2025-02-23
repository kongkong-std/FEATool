#include "../include/main.h"

void FourCDat2MSHFile(const char *path, int num_node, const FourCDatNodeData *node_data,
                      int num_ele, const FourCDatEleData *ele_data,
                      int num_edge_ele, const FourCDatEdgeEleData *edge_ele_data)
{
    FILE *fp = fopen(path, "wb");
    assert(fp);

    fprintf(fp, "$MeshFormat\n"
                "2.2 0 8\n"
                "$EndMeshFormat\n");

    // nodes
    fprintf(fp, "$Nodes\n");
    fprintf(fp, "%d\n", num_node);
    for (int index = 0; index < num_node; ++index)
    {
        fprintf(fp, "%d %021.16le %021.16le %021.16le\n", (node_data + index)->idx,
                (node_data + index)->x,
                (node_data + index)->y,
                (node_data + index)->z);
    }
    fprintf(fp, "$EndNodes\n");

    // elements
    fprintf(fp, "$Elements\n");
    fprintf(fp, "%d\n", num_edge_ele + num_ele);
    for (int index = 0; index < num_edge_ele; ++index)
    {
        fprintf(fp, "%d " /*index*/
                    "%d " /*mesh type*/
                    "%d " /*tag num*/,
                (edge_ele_data + index)->idx,
                (edge_ele_data + index)->msh_type, 2);
        for (int index_i = 0; index_i < 2; ++index_i)
        {
            fprintf(fp, "%d ", (edge_ele_data + index)->edge_flag);
        }
        for (int index_i = 0; index_i < (edge_ele_data + index)->num_cell; ++index_i)
        {
            fprintf(fp, "%d ", (edge_ele_data + index)->node_list[index_i]);
        }
        fprintf(fp, "\n");
    }
    for (int index = 0; index < num_ele; ++index)
    {
        fprintf(fp, "%d "
                    "%d "
                    "%d ",
                (ele_data + index)->idx + num_edge_ele,
                (ele_data + index)->msh_type, 2);
        for (int index_i = 0; index_i < 2; ++index_i)
        {
            fprintf(fp, "%d ", (ele_data + index)->omega_flag);
        }
        for (int index_i = 0; index_i < (ele_data + index)->num_cell; ++index_i)
        {
            fprintf(fp, "%d ", (ele_data + index)->node_list[index_i]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "$EndElements");

    fclose(fp);
}

void FourCDatFileProcess(const char *path, int *ptr_num_node, FourCDatNodeData **node_data,
                         int *ptr_num_ele, FourCDatEleData **ele_data,
                         int *ptr_num_edge_ele, FourCDatEdgeEleData **edge_ele_data)
{
    // printf("path to dat file: %s\n", path);
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[MAX_SIZE];
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "PROBLEM SIZE"))
        {
            fscanf(fp, "%*s%d", ptr_num_ele);
            fscanf(fp, "%*s%d", ptr_num_node);
            break;
        }
    }
    int num_ele = *ptr_num_ele, num_node = *ptr_num_node;
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

    // oemga flag data
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "DSURF-NODE TOPOLOGY"))
        {
            // NODE 1 DSURFACE 1
            for (int index = 0; index < num_node; ++index)
            {
                int tmp_node_idx = 0, omega_flag = 0;
                fscanf(fp, "%*s%d%*s%d", &tmp_node_idx, &omega_flag);
                for (int index_node = 0; index_node < num_node; ++index_node)
                {
                    if (tmp_node_idx == (*node_data + index_node)->idx)
                    {
                        (*node_data + index_node)->omega_flag = omega_flag;
                        break;
                    }
                }
            }
            break;
        }
    }

    rewind(fp);

    // edge element data
    int cnt_edge = 0;
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "DNODE-NODE TOPOLOGY"))
        {
            while (1)
            {
                fgets(buffer, MAX_SIZE, fp);
                int tmp_node_1 = 0, tmp_node_2 = 0;
                if (sscanf(buffer, "%*s%d%*s%d", &tmp_node_1, &tmp_node_2) == 2)
                {
                    // printf("%d\t%d\n", tmp_node_1, tmp_node_2);
                    ++cnt_edge;
                }
                else
                {
                    break;
                }
            }
            break;
        }
    }

    rewind(fp);

    int *edge_flag = NULL;
    edge_flag = (int *)malloc(cnt_edge * sizeof(int));
    assert(edge_flag);

    int tmp_cnt_edge = 0;
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "DNODE-NODE TOPOLOGY"))
        {
            while (1)
            {
                fgets(buffer, MAX_SIZE, fp);
                int tmp_node_1 = 0;
                if (sscanf(buffer, "%*s%d%*s%d", &tmp_node_1, edge_flag + tmp_cnt_edge) == 2)
                {
                    // printf("%d\t%d\n", tmp_node_1, tmp_node_2);
                    ++tmp_cnt_edge;
                }
                else
                {
                    break;
                }
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < cnt_edge; ++index)
    {
        printf("%d\n", edge_flag[index]);
    }
#endif

    rewind(fp);

    int *num_edge_flag = NULL;
    num_edge_flag = (int *)malloc(cnt_edge * sizeof(int));
    assert(num_edge_flag);
    memset(num_edge_flag, 0, cnt_edge * sizeof(int));

    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "DLINE-NODE TOPOLOGY"))
        {
            while (1)
            {
                fgets(buffer, MAX_SIZE, fp);
                int tmp_edge_flag = 0;
                if (sscanf(buffer, "%*s%*s%*s%d", &tmp_edge_flag) == 1)
                {
                    for (int index_cnt_edge = 0; index_cnt_edge < cnt_edge; ++index_cnt_edge)
                    {
                        if (edge_flag[index_cnt_edge] == tmp_edge_flag)
                        {
                            ++num_edge_flag[index_cnt_edge];
                        }
                    }
                }
                else
                {
                    break;
                }
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < cnt_edge; ++index)
    {
        printf("%d\n", num_edge_flag[index]);
    }
#endif

    int **edge_node = NULL;
    edge_node = (int **)malloc(cnt_edge * sizeof(int *));
    for (int index = 0; index < cnt_edge; ++index)
    {
        *(edge_node + index) = (int *)malloc(num_edge_flag[index] * sizeof(int));
    }

    rewind(fp);
    while (!feof(fp))
    {
        fgets(buffer, MAX_SIZE, fp);
        if (strstr(buffer, "DLINE-NODE TOPOLOGY"))
        {
            for (int index = 0; index < cnt_edge; ++index)
            {
                for (int index_i = 0; index_i < num_edge_flag[index]; ++index_i)
                {
                    fscanf(fp, "%*s%d%*s%*d", edge_node[index] + index_i);
                }
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < cnt_edge; ++index)
    {
        for (int index_i = 0; index_i < num_edge_flag[index]; ++index_i)
        {
            printf("%d\t", edge_node[index][index_i]);
        }
        putchar('\n');
    }
#endif

    int num_edge_ele = 0;
    for (int index = 0; index < cnt_edge; ++index)
    {
        num_edge_ele += num_edge_flag[index] - 1;
    }
    // printf("number of edge element: %d\n", num_edge_ele);
    *ptr_num_edge_ele = num_edge_ele;

    *edge_ele_data = (FourCDatEdgeEleData *)malloc(num_edge_ele * sizeof(FourCDatEdgeEleData));
    assert(*edge_ele_data);

    int cnt_edge_ele = 0;
    for (int index = 0; index < cnt_edge; ++index)
    {
        for (int index_i = 0; index_i < num_edge_flag[index] - 1; ++index_i)
        {
            (*edge_ele_data + cnt_edge_ele)->idx = cnt_edge_ele + 1;
            (*edge_ele_data + cnt_edge_ele)->edge_flag = edge_flag[index];
            (*edge_ele_data + cnt_edge_ele)->msh_type = 1;
            (*edge_ele_data + cnt_edge_ele)->num_cell = 2;
            (*edge_ele_data + cnt_edge_ele)->node_list[0] = edge_node[index][index_i];
            (*edge_ele_data + cnt_edge_ele)->node_list[1] = edge_node[index][index_i + 1];
            ++cnt_edge_ele;
        }
    }

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
                const char *ptr_buffer = buffer;
                int buffer_offset;
                sscanf(buffer, "%d%*s%s%n", &((*ele_data + index)->idx),
                       (*ele_data + index)->cell_type, &buffer_offset);
                ptr_buffer += buffer_offset;

                int num_node_cell_type = FourCNumNodeCellType((*ele_data)[index].cell_type, &((*ele_data + index)->msh_type));
                (*ele_data + index)->num_cell = num_node_cell_type;
                // printf(">>>> num node = %d, msh type = %d >>>>\n", num_node_cell_type, (*ele_data)[index].msh_type);
                (*ele_data + index)->node_list = (int *)malloc(num_node_cell_type * sizeof(int));
                assert((*ele_data + index)->node_list);

                for (int index_buffer = 0; index_buffer < num_node_cell_type; ++index_buffer)
                {
                    sscanf(ptr_buffer, "%d%n", (*ele_data + index)->node_list + index_buffer, &buffer_offset);
                    ptr_buffer += buffer_offset;
                }
            }
            break;
        }
    }
#if 0
    for (int index = 0; index < num_ele; ++index)
    {
        printf("%d\t%d\t", (*ele_data + index)->idx, (*ele_data + index)->msh_type);
        for (int index_ele = 0; index_ele < (*ele_data + index)->num_cell; ++index_ele)
        {
            printf("%d\t", (*ele_data + index)->node_list[index_ele]);
        }
        putchar('\n');
    }
#endif

#if 1
    for (int index = 0; index < num_ele; ++index)
    {
        for (int index_i = 0; index_i < num_node; ++index_i)
        {
            if ((*node_data + index_i)->idx == (*ele_data + index)->node_list[0])
            {
                (*ele_data + index)->omega_flag = (*node_data + index_i)->omega_flag;
            }
        }
    }
#endif

    fclose(fp);

    // free memory
    for (int index = 0; index < cnt_edge; ++index)
    {
        free(edge_node[index]);
    }
    free(edge_node);
    free(num_edge_flag);
    free(edge_flag);
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
