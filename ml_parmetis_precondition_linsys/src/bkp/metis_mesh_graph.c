#include "../include/main.h"

int NumNodeEleTypeMap(int ele_type /*element type*/)
{
    int value = 0;

    switch (ele_type)
    {
    case 1:
        // 2-node line
        value = 2;
        break;

    case 2:
        // 3-node triangle
        value = 3;
        break;

    case 3:
        // 4-node quadrangle
        value = 4;
        break;

    case 4:
        // 4-node tetrahedron
        value = 4;
        break;

    case 5:
        // 8-node hexahedron
        value = 8;
        break;

    case 6:
        // 6-node prism
        value = 6;
        break;

    case 7:
        // 5-node pyramid
        value = 5;
        break;

    default:
        break;
    }

    return value;
}

void MetisFileProcessGmsh(const char *path /*path to gmsh file*/,
                          DataGmsh *data /*gmsh data pointer*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[BUF_MAX_SIZE];
    Flag_Data_Block current_status = NONE;
    int line_node = 0, line_ele_bd = 0, line_ele_in = 0;
    int nn_ele_bd = 0, nn_ele_in = 0;

    while (fgets(buffer, BUF_MAX_SIZE, fp))
    {
        char tmp_buffer[BUF_MAX_SIZE];
        if (buffer[0] == '$')
        {
            if (strstr(buffer, "$MeshFormat"))
            {
                current_status = MESH_FORMAT;
            }
            else if (strstr(buffer, "$PhysicalNames"))
            {
                current_status = PHYSICAL_NAMES;
            }
            else if (strstr(buffer, "$Nodes"))
            {
                current_status = NODES;
            }
            else if (strstr(buffer, "$Elements"))
            {
                current_status = ELEMENTS;
            }
            else
            {
                current_status = NONE;
            }
            continue;
        }

        switch (current_status)
        {
        case NODES:
        {
            // printf("%s", buffer);
            memcpy(tmp_buffer, buffer, BUF_MAX_SIZE);
            tmp_buffer[strcspn(tmp_buffer, "\n")] = '\0';
            int tmp_cnt = 0;
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &(data->nn));
                data->coordinates = (double *)malloc(data->nn * 3 * sizeof(double));
                assert(data->coordinates);
            }
            else
            {
                sscanf(buffer, "%*d%lf%lf%lf", data->coordinates + 3 * line_node,
                       data->coordinates + 3 * line_node + 1,
                       data->coordinates + 3 * line_node + 2);
                ++line_node;
            }

            break;
        }
        case ELEMENTS:
        {
            // printf("%s", buffer);
            memcpy(tmp_buffer, buffer, BUF_MAX_SIZE);
            tmp_buffer[strcspn(tmp_buffer, "\n")] = '\0';
            int tmp_cnt = 0;
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &(data->ne));
                data->eptr_bd = (idx_t *)malloc(data->ne * sizeof(idx_t));
                data->eptr_in = (idx_t *)malloc(data->ne * sizeof(idx_t));
                nn_ele_bd = NumNodeEleTypeMap(1); // segment element
                nn_ele_in = NumNodeEleTypeMap(2); // triangle element
                data->eind_bd = (idx_t *)malloc(data->ne * nn_ele_bd * sizeof(idx_t));
                data->eind_in = (idx_t *)malloc(data->ne * nn_ele_in * sizeof(idx_t));
                assert(data->eptr_bd &&
                       data->eptr_in &&
                       data->eind_bd &&
                       data->eind_in);
                data->eptr_bd[0] = 0;
                data->eptr_in[0] = 0;
            }
            else
            {
                int ele_type = 0;
                sscanf(buffer, "%*d%d", &ele_type);
                if (ele_type == 1)
                {
                    data->eptr_bd[line_ele_bd + 1] = data->eptr_bd[line_ele_bd] + nn_ele_bd;

                    // last nn_ele_bd number in buffer
                    char *token;
                    char *last2[2] = {NULL, NULL}; // 保存最后两个token

                    token = strtok(buffer, " \t");
                    while (token != NULL)
                    {
                        // 滚动保存最后两个
                        last2[0] = last2[1];
                        last2[1] = token;
                        token = strtok(NULL, " \t");
                    }

                    int second_last = atoi(last2[0]);
                    int last = atoi(last2[1]);

                    // change to 0-base
                    data->eind_bd[2 * line_ele_bd] = second_last - 1;
                    data->eind_bd[2 * line_ele_bd + 1] = last - 1;

                    ++line_ele_bd;
                }
                else if (ele_type == 2)
                {
                    data->eptr_in[line_ele_in + 1] = data->eptr_in[line_ele_in] + nn_ele_in;

                    // last nn_ele_bd number in buffer
                    char *token;
                    char *last3[3] = {NULL, NULL, NULL}; // 保存最后两个token

                    token = strtok(buffer, " \t");
                    while (token != NULL)
                    {
                        // 滚动保存最后两个
                        last3[0] = last3[1];
                        last3[1] = last3[2];
                        last3[2] = token;
                        token = strtok(NULL, " \t");
                    }

                    int third_last = atoi(last3[0]);
                    int second_last = atoi(last3[1]);
                    int last = atoi(last3[2]);

                    // change to 0-base
                    data->eind_in[3 * line_ele_in] = third_last - 1;
                    data->eind_in[3 * line_ele_in + 1] = second_last - 1;
                    data->eind_in[3 * line_ele_in + 2] = last - 1;

                    ++line_ele_in;
                }
            }
            break;
        }
        default:
        {
            break;
        }
        }
    }

    data->ne_bd = line_ele_bd;
    data->ne_in = line_ele_in;
    data->nne_bd = nn_ele_bd;
    data->nne_in = nn_ele_in;
    data->nparts = data->nn / 4;

    data->epart_in = (idx_t *)malloc(data->ne_in * sizeof(idx_t));
    data->npart_in = (idx_t *)malloc(data->nn * sizeof(idx_t));
    assert(data->epart_in && data->npart_in);

    fclose(fp);
}
