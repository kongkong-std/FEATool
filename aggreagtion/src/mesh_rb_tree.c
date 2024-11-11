#if 0
#include "../include/main.h"

// 创建一个新的红黑树节点
MeshNodeRBNode *CreateMeshNodeRBNode(MeshNode data)
{
    MeshNodeRBNode *newNode = (MeshNodeRBNode *)malloc(sizeof(MeshNodeRBNode));
    newNode->data = data;
    newNode->color = RED; // 新插入的节点默认是红色
    newNode->left = newNode->right = newNode->parent = NULL;
    return newNode;
}

// 创建红黑树
MeshNodeRBTree *CreateMeshNodeRBTree(void)
{
    MeshNodeRBTree *tree = (MeshNodeRBTree *)malloc(sizeof(MeshNodeRBTree));
    if (!tree)
    {
        fprintf(stderr, "Memory allocation failed for RBTree\n");
        exit(1);
    }

    // 创建哨兵节点并初始化
    tree->TNULL = CreateMeshNodeRBNode((MeshNode){0, 0, 0, 0}); // 哨兵节点
    tree->TNULL->color = BLACK;
    tree->TNULL->left = tree->TNULL;
    tree->TNULL->right = tree->TNULL;
    tree->TNULL->parent = tree->TNULL;

    tree->root = tree->TNULL; // 树初始化为空
    return tree;
}

// 左旋操作
void LeftRotate(MeshNodeRBTree *tree, MeshNodeRBNode *x)
{
    MeshNodeRBNode *y = x->right;
    x->right = y->left;
    if (y->left != tree->TNULL)
    {
        y->left->parent = x;
    }
    y->parent = x->parent;
    if (x->parent == NULL)
    {
        tree->root = y; // 如果 x 是根节点，y 成为新的根节点
    }
    else if (x == x->parent->left)
    {
        x->parent->left = y;
    }
    else
    {
        x->parent->right = y;
    }
    y->left = x;
    x->parent = y;
}

// 右旋操作
void RightRotate(MeshNodeRBTree *tree, MeshNodeRBNode *x)
{
    MeshNodeRBNode *y = x->left;
    x->left = y->right;
    if (y->right != tree->TNULL)
    {
        y->right->parent = x;
    }
    y->parent = x->parent;
    if (x->parent == NULL)
    {
        tree->root = y;
    }
    else if (x == x->parent->right)
    {
        x->parent->right = y;
    }
    else
    {
        x->parent->left = y;
    }
    y->right = x;
    x->parent = y;
}

// 插入修复
void FixInsert(MeshNodeRBTree *tree, MeshNodeRBNode *k)
{
    MeshNodeRBNode *u;
    while (k->parent->color == RED)
    {
        if (k->parent == k->parent->parent->right)
        {
            u = k->parent->parent->left;
            if (u->color == RED)
            {
                u->color = BLACK;
                k->parent->color = BLACK;
                k->parent->parent->color = RED;
                k = k->parent->parent;
            }
            else
            {
                if (k == k->parent->left)
                {
                    k = k->parent;
                    RightRotate(tree, k);
                }
                k->parent->color = BLACK;
                k->parent->parent->color = RED;
                LeftRotate(tree, k->parent->parent);
            }
        }
        else
        {
            u = k->parent->parent->right;
            if (u->color == RED)
            {
                u->color = BLACK;
                k->parent->color = BLACK;
                k->parent->parent->color = RED;
                k = k->parent->parent;
            }
            else
            {
                if (k == k->parent->right)
                {
                    k = k->parent;
                    LeftRotate(tree, k);
                }
                k->parent->color = BLACK;
                k->parent->parent->color = RED;
                RightRotate(tree, k->parent->parent);
            }
        }
        if (k == tree->root)
        {
            break;
        }
    }
    tree->root->color = BLACK;
}

void InsertRBTree(MeshNodeRBTree *tree, MeshNode data)
{
    MeshNodeRBNode *node = CreateMeshNodeRBNode(data);
    MeshNodeRBNode *yNode = NULL;
    MeshNodeRBNode *xNode = tree->root;

    // 通过循环找到合适的插入位置
    while (xNode != tree->TNULL)
    {
        yNode = xNode;
        if (node->data.node_idx < xNode->data.node_idx)
        {
            xNode = xNode->left;
        }
        else
        {
            xNode = xNode->right;
        }
    }

    node->parent = yNode; // 设置父节点
    if (yNode == NULL)
    {
        tree->root = node; // 如果树为空，新节点为根节点
    }
    else if (node->data.node_idx < yNode->data.node_idx)
    {
        yNode->left = node; // 否则设置为左子节点
    }
    else
    {
        yNode->right = node; // 否则设置为右子节点
    }

    node->left = tree->TNULL;  // 左子节点为哨兵节点
    node->right = tree->TNULL; // 右子节点为哨兵节点
    node->color = RED;         // 新节点默认是红色

    // 插入后修复红黑树
    FixInsert(tree, node);
}

MeshNodeRBNode *SearchRBTree(MeshNodeRBTree *tree, int node_idx)
{
    MeshNodeRBNode *current = tree->root;
    while (current != tree->TNULL && node_idx != current->data.node_idx)
    {
        if (node_idx < current->data.node_idx)
        {
            current = current->left;
        }
        else
        {
            current = current->right;
        }
    }
    return current;
}
#endif
