#ifndef BPLUSTREE_H_INCLUDED
#define BPLUSTREE_H_INCLUDED

#pragma once

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define M (4)
#define LIMIT_M_2 (M % 2 ? (M + 1) / 2 : M / 2)

typedef struct BPlusNode *BPlusTree, *Position;
typedef int KeyType;
struct BPlusNode {
  int KeyNum;
  KeyType Key[M + 1];
  BPlusTree Children[M + 1];
  BPlusTree Next;
};

/* 初始化 */
BPlusTree Initialize();
/* 插入 */
BPlusTree Insert(BPlusTree T, KeyType Key);
/* 删除 */
BPlusTree Remove(BPlusTree T, KeyType Key);
/* 销毁 */
BPlusTree Destroy(BPlusTree T);
/* 遍历节点 */
void Travel(BPlusTree T);
/* 遍历树叶节点的数据 */

#endif