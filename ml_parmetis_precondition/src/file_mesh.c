#include "../include/main.h"

static int CountTrueMatAdj(const bool *a, int n)
{
    int value = 0;

    for (int index = 0; index < n; ++index)
    {
        if (a[index])
        {
            ++value;
        }
    }

    return value;
}

int GlobalGraphCSRAdjGenerator(const DataMesh *mesh_data /*mesh data*/,
                               AdjDataMesh **graph_data /*csr graph data*/)
{
    *graph_data = (AdjDataMesh *)malloc(sizeof(AdjDataMesh));
    assert(*graph_data);

    (*graph_data)->nn = mesh_data->nn;
    (*graph_data)->dim = mesh_data->dim;
    (*graph_data)->nparts = (*graph_data)->nn / 4;

    (*graph_data)->coordinates = (double *)malloc((*graph_data)->dim *
                                                  (*graph_data)->nn *
                                                  sizeof(double));
    assert((*graph_data)->coordinates);

    memcpy((*graph_data)->coordinates,
           mesh_data->coordinates,
           mesh_data->dim * mesh_data->nn * sizeof(double));

    int nn = mesh_data->nn;
    int ne = mesh_data->ne_shell;
    DataMeshEle *ele_data = mesh_data->ele_shell;

    (*graph_data)->xadj = (idx_t *)calloc(nn + 1, sizeof(idx_t));
    bool *mat_adj_tmp = (bool *)calloc(nn, sizeof(bool));
    int **adjncy_tmp = (int **)malloc(nn * sizeof(int *));
    int *adjncy_tmp_len = (int *)malloc(nn * sizeof(int)); // record length of adjncy_tmp[i]

    assert((*graph_data)->xadj && mat_adj_tmp && adjncy_tmp && adjncy_tmp_len);

    for (int index = 0; index < nn; ++index)
    {
        for (int index_e = 0; index_e < ne; ++index_e)
        {
            // triangle element
            if (ele_data[index_e].ele_node[0] == index ||
                ele_data[index_e].ele_node[1] == index ||
                ele_data[index_e].ele_node[2] == index)
            {
                mat_adj_tmp[ele_data[index_e].ele_node[0]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[1]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[2]] = true;
            }
        }
        mat_adj_tmp[index] = false;

        int cnt_tmp = CountTrueMatAdj(mat_adj_tmp, nn);

        (*graph_data)->xadj[index + 1] = (*graph_data)->xadj[index] + cnt_tmp;

        adjncy_tmp[index] = (int *)malloc(cnt_tmp * sizeof(int));
        adjncy_tmp_len[index] = cnt_tmp;
        int cnt_adjncy_tmp = 0;
        for (int index_mat_adj_tmp = 0; index_mat_adj_tmp < nn; ++index_mat_adj_tmp)
        {
            if (mat_adj_tmp[index_mat_adj_tmp])
            {
                adjncy_tmp[index][cnt_adjncy_tmp] = index_mat_adj_tmp;
                ++cnt_adjncy_tmp;
            }
        }

        memset(mat_adj_tmp, 0, nn * sizeof(bool));
    }

    (*graph_data)->adjncy = (idx_t *)malloc((*graph_data)->xadj[nn] * sizeof(idx_t));
    assert((*graph_data)->adjncy);

    int index_adjncy = 0;
    for (int index_i = 0; index_i < nn; ++index_i)
    {
        for (int index_j = 0; index_j < adjncy_tmp_len[index_i]; ++index_j)
        {
            (*graph_data)->adjncy[index_adjncy] = adjncy_tmp[index_i][index_j];
            ++index_adjncy;
        }
    }

    // free memory
    free(adjncy_tmp_len);
    for (int index = 0; index < nn; ++index)
    {
        free(adjncy_tmp[index]);
    }
    free(adjncy_tmp);
    free(mat_adj_tmp);

    return 0;
}

static int NumNodeEleTypeMap(int ele_type /*element type*/)
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

int FileProcessMesh(const char *path /*path to mesh file*/, DataMesh *mesh_data /*mesh data*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[BUF_MAX_SIZE];
    FlagDataBlockGmsh current_status = NONE;

    int line_node = 0, line_ele = 0;
    DataMeshEle *ele_tmp = NULL;

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
            /* code */
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
                sscanf(buffer, "%d", &mesh_data->nn);
                mesh_data->dim = 3;
                mesh_data->coordinates = (double *)malloc(mesh_data->dim * mesh_data->nn * sizeof(double));
                assert(mesh_data->coordinates);
            }
            else
            {
                sscanf(buffer, " %*d %lf %lf %lf ", mesh_data->coordinates + mesh_data->dim * line_node,
                       mesh_data->coordinates + mesh_data->dim * line_node + 1,
                       mesh_data->coordinates + mesh_data->dim * line_node + 2);
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
            int tmp_array[256];
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                tmp_array[tmp_cnt] = atoi(tmp_token);
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &mesh_data->ne);
                ele_tmp = (DataMeshEle *)malloc(mesh_data->ne * sizeof(DataMeshEle));
                assert(ele_tmp);
            }
            else
            {
                sscanf(buffer, " %*d %d ", &ele_tmp[line_ele].ele_type);
                ele_tmp[line_ele].num_ele_node = NumNodeEleTypeMap(ele_tmp[line_ele].ele_type);
                ele_tmp[line_ele].ele_node = (int *)malloc(ele_tmp[line_ele].num_ele_node * sizeof(int));
                assert(ele_tmp[line_ele].ele_node);

                for (int index = 0; index < ele_tmp[line_ele].num_ele_node; ++index)
                {
                    ele_tmp[line_ele].ele_node[index] = tmp_array[tmp_cnt - ele_tmp[line_ele].num_ele_node + index] - 1; // 0-base
                }

                ++line_ele;
            }

            break;
        }

        default:
            break;
        }
    }

    fclose(fp);

#if 0
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        printf("element %d:\t%d\t", index, ele_tmp[index].ele_type);
        for (int index_i = 0; index_i < ele_tmp[index].num_ele_node; ++index_i)
        {
            printf("%d\t", ele_tmp[index].ele_node[index_i]);
        }
        putchar('\n');
    }
#endif // element data

    int cnt_shell = 0;
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        if (ele_tmp[index].ele_type == 2)
        {
            // triangle mesh
            ++cnt_shell;
        }
    }

    mesh_data->ne_shell = cnt_shell;
    mesh_data->ele_shell = (DataMeshEle *)malloc(mesh_data->ne_shell * sizeof(DataMeshEle));
    assert(mesh_data->ele_shell);

    cnt_shell = 0;
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        if (ele_tmp[index].ele_type == 2)
        {
            mesh_data->ele_shell[cnt_shell].ele_type = ele_tmp[index].ele_type;
            mesh_data->ele_shell[cnt_shell].num_ele_node = ele_tmp[index].num_ele_node;

            mesh_data->ele_shell[cnt_shell].ele_node = (int *)malloc(mesh_data->ele_shell[cnt_shell].num_ele_node * sizeof(int));
            assert(mesh_data->ele_shell[cnt_shell].ele_node);

            memcpy(mesh_data->ele_shell[cnt_shell].ele_node, ele_tmp[index].ele_node, mesh_data->ele_shell[cnt_shell].num_ele_node * sizeof(int));

            ++cnt_shell;
        }
    }

    // free memory
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        free(ele_tmp[index].ele_node);
    }
    free(ele_tmp);

    return 0;
}
