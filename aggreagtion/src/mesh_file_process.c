#include "../include/main.h"

int NodeListSizeMapping(int ele_type)
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

void MeshFileProcess(const char *path,
                     int label_bound, int label_omega,
                     GenericList *data_list_phy_tag,
                     GenericList *data_list_node,
                     GenericList *data_list_ele_bound,
                     GenericList *data_list_ele_omega)
{
    int cnt_line = 0, flag_tag = 0, flag_node = 0, flag_ele = 0;
    FILE *fp = NULL;

    fp = fopen(path, "rb");
    assert(fp);
    char buffer[MAX_SIZE];

    while (fgets(buffer, MAX_SIZE, fp) != NULL)
    {
        // physical tag
        if (strstr(buffer, "$PhysicalNames"))
        {
            // printf("%s", buffer);
            cnt_line = 0;
            flag_tag = 1;
            continue;
        }
        if (flag_tag == 1)
        {
            if (strstr(buffer, "$EndPhysicalNames"))
            {
                // printf("%s", buffer);
                flag_tag = 0;
                continue;
            }
            // printf("%s", buffer);
            if (cnt_line != 0)
            {
                MeshPhysicalTag *data_phy_tag = (MeshPhysicalTag *)malloc(sizeof(MeshPhysicalTag));
                assert(data_phy_tag);
                sscanf(buffer, "%d%d", &data_phy_tag->tag_dim, &data_phy_tag->tag_idx);
                AddNodeToList(data_list_phy_tag, data_phy_tag, TYPE_PHYSICAL_TAG);
            }
            ++cnt_line;
        }

        // node tag
        if (strstr(buffer, "$Nodes"))
        {
            // printf("%s", buffer);
            cnt_line = 0;
            flag_node = 1;
            continue;
        }
        if (flag_node == 1)
        {
            if (strstr(buffer, "$EndNodes"))
            {
                // printf("%s", buffer);
                flag_node = 0;
                continue;
            }
            // printf("%s", buffer);
            if (cnt_line != 0)
            {
                MeshNode *data_node = (MeshNode *)malloc(sizeof(MeshNode));
                assert(data_node);
                sscanf(buffer, "%d%lf%lf%lf", &data_node->node_idx,
                       &data_node->node_x, &data_node->node_y, &data_node->node_z);
                AddNodeToList(data_list_node, data_node, TYPE_NODE);
            }
            ++cnt_line;
        }

        // element tag
        if (strstr(buffer, "$Elements"))
        {
            // printf("%s", buffer);
            cnt_line = 0;
            flag_ele = 1;
            continue;
        }
        if (flag_ele == 1)
        {
            if (strstr(buffer, "$EndElements"))
            {
                // printf("%s", buffer);
                flag_ele = 0;
                continue;
            }
            // printf("%s", buffer);
            if (cnt_line != 0)
            {
                MeshElement *data_ele = (MeshElement *)malloc(sizeof(MeshElement));
                assert(data_ele);
                int buffer_offset = sscanf(buffer, "%d%d%d", &data_ele->ele_idx,
                                           &data_ele->ele_type,
                                           &data_ele->ele_num_tag);
                int node_lst_size = NodeListSizeMapping(data_ele->ele_type);
                data_ele->ele_tag_idx = (int *)malloc(data_ele->ele_num_tag * sizeof(int));
                assert(data_ele->ele_tag_idx);

                char *ptr = buffer;
                for (int index = 0; index < buffer_offset; ++index)
                {
                    ptr = strchr(ptr, ' ') + 1; // 找到下一个空格并移动到下一个数字
                }

                for (int index = 0; index < data_ele->ele_num_tag; ++index)
                {
                    sscanf(ptr, "%d", data_ele->ele_tag_idx + index);
                    ptr = strchr(ptr, ' ') + 1;
                }

                data_ele->ele_node_lst = (int *)malloc(node_lst_size * sizeof(int));
                assert(data_ele->ele_node_lst);
                for (int index = 0; index < node_lst_size; ++index)
                {
                    sscanf(ptr, "%d", data_ele->ele_node_lst + index);
                    ptr = strchr(ptr, ' ') + 1;
                }
                if (data_ele->ele_type == label_bound)
                {
                    AddNodeToList(data_list_ele_bound, data_ele, TYPE_ELEMENT);
                }
                else if (data_ele->ele_type == label_omega)
                {
                    AddNodeToList(data_list_ele_omega, data_ele, TYPE_ELEMENT);
                }
            }
            ++cnt_line;
        }
    }

    fclose(fp);
}
