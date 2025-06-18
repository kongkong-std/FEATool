#include "../include/main.h"

void ClearList(GenericList *list)
{
    // parameter non-null check
    if (list)
    {
        ListNode *current = list->head;
        while (current != NULL)
        {
            ListNode *temp = current;
            current = current->next;
            NodeType type = temp->type;
            if (type == TYPE_ELEMENT)
            {
                MeshElement *element = (MeshElement *)temp->data;
                free(element->ele_tag_idx);
                free(element->ele_node_lst);
                free(element);
            }
            else
            {
                free(temp->data); // 释放节点数据
            }
            free(temp); // 释放节点
        }
        list->head = NULL; // 清空头指针
        list->tail = NULL; // 清空尾指针
        list->size = 0;    // 重置大小
    }
}

void AddNodeToList(GenericList *list, void *data, NodeType type)
{
    ListNode *newNode = (ListNode *)malloc(sizeof(ListNode));
    newNode->data = data; // 设置节点数据
    newNode->next = NULL; // 新节点的下一个指针为空
    newNode->type = type;

    if (list->tail)
    {
        list->tail->next = newNode; // 如果尾指针存在，将新节点添加到尾部
    }
    else
    {
        list->head = newNode; // 如果是第一个节点，设置头指针
    }
    list->tail = newNode; // 更新尾指针
    list->size++;         // 更新链表大小
}

void InitializeList(GenericList *list)
{
    list->size = 0;
    list->head = NULL;
    list->tail = NULL;
}

int IsListEmpty(GenericList *list)
{
    int value = 0;
    if (list->head == NULL &&
        list->tail == NULL &&
        list->size == 0)
    {
        value = 1;
    }

    return value;
}
