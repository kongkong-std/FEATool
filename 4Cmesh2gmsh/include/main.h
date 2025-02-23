#ifndef MAIN_H_
#define MAIN_H_

// header
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAX_SIZE 1024

// struct
typedef struct
{
    /* data */
    int idx;
    double x, y, z;
    int omega_flag;
} FourCDatNodeData;

typedef struct
{
    /* data */
    int idx;
    char cell_type[32];
    int *node_list;
    int msh_type; // .msh file element type
    int num_cell;
    int omega_flag;
} FourCDatEleData;

typedef struct
{
    /* data */
    int idx;
    int edge_flag;
    int node_list[2];
    int msh_type;
    int num_cell;
} FourCDatEdgeEleData;

// function prototype
void FourCDatFileProcess(const char * /*path to 4c dat file*/, int * /*number of nodes*/, FourCDatNodeData ** /*node data*/,
                         int * /*number of elements*/, FourCDatEleData ** /*element data*/,
                         int * /*number of edge elements*/, FourCDatEdgeEleData ** /*edge element data*/);
int FourCNumNodeCellType(const char * /*number of nodes of cell type*/, int * /*msh file type*/);
void FourCDat2MSHFile(const char * /*path to msh file*/, int /*number of nodes*/, const FourCDatNodeData * /*node data*/,
                      int /*number of elements*/, const FourCDatEleData * /*element data*/,
                      int /*number of edge elements*/, const FourCDatEdgeEleData * /*edge element data*/);

#endif // main.h