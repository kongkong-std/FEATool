#include "../include/main.h"

void MapAggregationProlongationOperator(MeshGraph *graph_agreagtion,
                                        int *map_row, int *map_column, int size)
{
    int cnt = 0;
    for (int index = 0; index < graph_agreagtion->size; ++index)
    {
        MeshGraphAdjNode *curr_node = graph_agreagtion->array[index].head;
        while (curr_node)
        {
            // 1-base
            *(map_row + cnt) = curr_node->node->node_idx;
            *(map_column + cnt) = index + 1;
            curr_node = curr_node->next;
            ++cnt;
        }
    }
}

void PrintProlongationOperator(double ***P_operator, MeshGraph *graph_aggregation, int order_rbm)
{
    for (int index = 0; index < graph_aggregation->size; ++index)
    {
        for (int index_i = 0; index_i < graph_aggregation->array[index].size; ++index_i)
        {
            if (order_rbm == 1)
            {
                printf("%lf\t-\t-\t%lf\t%lf\t%lf\n"
                       "-\t%lf\t-\t%lf\t%lf\t%lf\n"
                       "-\t-\t%lf\t%lf\t%lf\t%lf\n"
                       "-\t-\t-\t%lf\t-\t-\n"
                       "-\t-\t-\t-\t%lf\t-\n"
                       "-\t-\t-\t-\t-\t%lf\n",
                       P_operator[index][index_i][0],
                       P_operator[index][index_i][1],
                       P_operator[index][index_i][2],
                       P_operator[index][index_i][3],
                       P_operator[index][index_i][4],
                       P_operator[index][index_i][5],
                       P_operator[index][index_i][6],
                       P_operator[index][index_i][7],
                       P_operator[index][index_i][8],
                       P_operator[index][index_i][9],
                       P_operator[index][index_i][10],
                       P_operator[index][index_i][11],
                       P_operator[index][index_i][12],
                       P_operator[index][index_i][13],
                       P_operator[index][index_i][14]);
            }
            else if (order_rbm == 2)
            {
                printf("%lf\t-\t-\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
                       "-\t%lf\t-\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
                       "-\t-\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
                       "-\t-\t-\t%lf\t-\t-\t-\t-\t-\n"
                       "-\t-\t-\t-\t%lf\t-\t-\t-\t-\n"
                       "-\t-\t-\t-\t-\t%lf\t-\t-\t-\n",
                       P_operator[index][index_i][0],
                       P_operator[index][index_i][1],
                       P_operator[index][index_i][2],
                       P_operator[index][index_i][3],
                       P_operator[index][index_i][4],
                       P_operator[index][index_i][5],
                       P_operator[index][index_i][6],
                       P_operator[index][index_i][7],
                       P_operator[index][index_i][8],
                       P_operator[index][index_i][9],
                       P_operator[index][index_i][10],
                       P_operator[index][index_i][11],
                       P_operator[index][index_i][12],
                       P_operator[index][index_i][13],
                       P_operator[index][index_i][14],
                       P_operator[index][index_i][15],
                       P_operator[index][index_i][16],
                       P_operator[index][index_i][17],
                       P_operator[index][index_i][18],
                       P_operator[index][index_i][19],
                       P_operator[index][index_i][20],
                       P_operator[index][index_i][21],
                       P_operator[index][index_i][22],
                       P_operator[index][index_i][23]);
            }
            putchar('\n');
        }
    }
}

void ConstructProlongationOperator(double ***P_operator,
                                   MeshGraph *graph_aggregation,
                                   MeshGraph *coarse_graph,
                                   int rbm_order)
{
    for (int index = 0; index < graph_aggregation->size; ++index)
    {
        MeshGraphAdjNode *data_node_coarse = coarse_graph->array[index].head;
        MeshGraphAdjNode *data_node_aggregation = graph_aggregation->array[index].head;
        for (int index_i = 0; index_i < graph_aggregation->array[index].size; ++index_i)
        {
            double node_coarse_x = data_node_coarse->node->node_x;
            double node_coarse_y = data_node_coarse->node->node_y;
            double node_coarse_z = data_node_coarse->node->node_z;
            double node_aggregation_x = data_node_aggregation->node->node_x;
            double node_aggregation_y = data_node_aggregation->node->node_y;
            double node_aggregation_z = data_node_aggregation->node->node_z;
            if (rbm_order == 1)
            {
                for (int index_tmp = 0; index_tmp < 9; index_tmp = index_tmp + 4)
                {
                    P_operator[index][index_i][index_tmp] = 1.;
                }
                for (int index_tmp = 12; index_tmp < 15; ++index_tmp)
                {
                    P_operator[index][index_i][index_tmp] = 1.;
                }
                for (int index_tmp = 1; index_tmp < 12; index_tmp = index_tmp + 5)
                {
                    P_operator[index][index_i][index_tmp] = 0.;
                }
                P_operator[index][index_i][2] = node_coarse_z - node_aggregation_z;  // z
                P_operator[index][index_i][3] = node_aggregation_y - node_coarse_y;  // -y
                P_operator[index][index_i][5] = node_aggregation_z - node_coarse_z;  // -z
                P_operator[index][index_i][7] = node_coarse_x - node_aggregation_x;  // x
                P_operator[index][index_i][9] = node_coarse_y - node_aggregation_y;  // y
                P_operator[index][index_i][10] = node_aggregation_x - node_coarse_x; // -x
            }
            else if (rbm_order == 2)
            {
                for (int index_tmp = 0; index_tmp < 15; index_tmp = index_tmp + 7)
                {
                    P_operator[index][index_i][index_tmp] = 1.;
                }
                for (int index_tmp = 21; index_tmp < 24; ++index_tmp)
                {
                    P_operator[index][index_i][index_tmp] = 1.;
                }
                for (int index_tmp = 1; index_tmp < 18; index_tmp = index_tmp + 8)
                {
                    P_operator[index][index_i][index_tmp] = 0.;
                }
                P_operator[index][index_i][2] = node_coarse_z - node_aggregation_z; // z
                P_operator[index][index_i][3] = node_aggregation_y - node_coarse_y; // -y
                P_operator[index][index_i][4] = (node_aggregation_z - node_coarse_z) *
                                                (node_coarse_x - node_aggregation_x); // -zx
                P_operator[index][index_i][5] = 0;                                    // 0
                P_operator[index][index_i][6] = (node_aggregation_z - node_coarse_z) *
                                                (node_coarse_y - node_aggregation_y); // -zy
                P_operator[index][index_i][7] = node_aggregation_z - node_coarse_z;   // -z
                P_operator[index][index_i][10] = node_coarse_x - node_aggregation_x;  // x
                P_operator[index][index_i][11] = 0;                                   // 0
                P_operator[index][index_i][12] = (node_aggregation_z - node_coarse_z) *
                                                 (node_coarse_y - node_aggregation_y); // -zy
                P_operator[index][index_i][13] = (node_aggregation_z - node_coarse_z) *
                                                 (node_coarse_x - node_aggregation_x); // -zx
                P_operator[index][index_i][15] = node_coarse_y - node_aggregation_y;   // y
                P_operator[index][index_i][16] = node_aggregation_x - node_coarse_x;   // -x
                P_operator[index][index_i][18] = (node_coarse_x - node_aggregation_x) *
                                                 (node_coarse_x - node_aggregation_x) / 2.; // x * x / 2
                P_operator[index][index_i][19] = (node_coarse_y - node_aggregation_y) *
                                                 (node_coarse_y - node_aggregation_y) / 2.; // y * y / 2
                P_operator[index][index_i][20] = (node_coarse_x - node_aggregation_x) *
                                                 (node_coarse_y - node_aggregation_y); // xy
            }
            data_node_aggregation = data_node_aggregation->next;
        }
    }
}
