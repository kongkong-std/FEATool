#ifndef MAIN_H_
#define MAIN_H_

// header
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// macro
#define MAX_SIZE 1024

// struct, linked list
typedef enum
{
    TYPE_PHYSICAL_TAG,
    TYPE_NODE,
    TYPE_ELEMENT
} NodeType;

typedef struct list_node
{
    void *data;             // 指向数据的指针
    struct list_node *next; // 指向下一个节点的指针
    NodeType type;
} ListNode;

typedef struct generic_list
{
    ListNode *head; // 链表头
    ListNode *tail; // 链表尾
    int size;       // 链表大小
} GenericList;

typedef struct
{
    /* data */
    int tag_dim; // tag dimensional
    int tag_idx; // tag index
} MeshPhysicalTag;

typedef struct
{
    /* data */
    int node_idx;                  // node index
    double node_x, node_y, node_z; // node coordinate
} MeshNode;

typedef struct
{
    /* data */
    int ele_idx;       // element index
    int ele_type;      // element type
    int ele_num_tag;   // element tag number
    int *ele_tag_idx;  // element tag index, size ele_num_tag
    int *ele_node_lst; // element node list, size map(ele_type)
} MeshElement;

// function prototype
/*
 * output to file of physical tag data
 */
void PrintFileMeshPhysicalTag(const char * /*path to destinal file of physical tag*/, GenericList * /*data*/);

/*
 * output to file of mesh noda data
 */
void PrintFileMeshNode(const char * /*path to destinal file of node*/, GenericList * /*data*/);

/*
 * output to file of mesh element data
 */
void PrintFileMeshElement(const char * /*path to destinal file of element*/, GenericList * /*data*/);

/*
 * linked list to output to file
 */
void PrintListFile(const char * /*path to destinal file*/, GenericList * /*linked list*/);

/*
 * element type to size mapping
 */
int NodeListSizeMapping(int /*element type*/);

/*
 * print physical tag
 */
void PrintMeshPhysicalTag(void * /*data*/);

/*
 * print mesh node
 */
void PrintMeshNode(void * /*data*/);

/*
 * print mesh element
 */
void PrintMeshElement(void * /*data*/);

/*
 * print linked list
 */
void PrintList(GenericList * /*linked list*/);

/*
 * clear linked list
 */
void ClearList(GenericList * /*linked list*/);

/*
 * add node to the list
 */
void AddNodeToList(GenericList * /*linked list*/, void * /*data*/, NodeType /*node type*/);

/*
 * check whether linked list is empty or not
 *     0: not empty
 *     1: empty
 */
int IsListEmpty(GenericList * /*linked list*/);

/*
 * initialize linked list
 */
void InitializeList(GenericList * /*linked list*/);

/*
 * mesh file process
 */
void MeshFileProcess(const char * /*path to mesh file*/,
                     int /*label of bound*/, int /*label of omega*/,
                     GenericList * /*generic data list*/,
                     GenericList * /*mesh node data*/,
                     GenericList * /*mesh bound element data*/,
                     GenericList * /*mesh omega element data*/);

#endif