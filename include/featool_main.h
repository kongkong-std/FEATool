#ifndef FEATOOL_MAIN_H_
#define FEATOOL_MAIN_H_

// hearder
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// struct
/*
 * 2d mesh element data
 *     element index
 *     node index: 1, 2, 3
 */
typedef struct mesh_element_2d
{
    /* data */
    int ele_index;
    int node_1, node_2, node_3;
} MeshElement2D;

/*
 * 3d mesh element data
 *     element index
 *     node index: 1, 2, 3, 4
 */
typedef struct mesh_element_3d
{
    /* data */
    int ele_index;
    int node_1, node_2, node_3, node_4;
} MeshElement3D;

/*
 * node data
 *     node index
 *     coordinate: x, y, z
 */
typedef struct mesh_node
{
    /* data */
    int node_index;
    double x, y, z;
} MeshNode;

#endif