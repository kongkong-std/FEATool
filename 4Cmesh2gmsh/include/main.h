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
} FourCDatNodeData;

typedef struct
{
    /* data */
    int idx;
    char cell_type[32];
    int *node_list;
    int msh_type; // .msh file element type
} FourCDatEleData;

// function prototype
void FourCDatFileProcess(const char * /*path to 4c dat file*/, FourCDatNodeData ** /*node data*/, FourCDatEleData ** /*element data*/);
int FourCNumNodeCellType(const char * /*number of nodes of cell type*/, int * /*msh file type*/);

#endif // main.h