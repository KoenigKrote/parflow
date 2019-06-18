#ifndef _BC_FLAT_H
#define _BC_FLAT_H

typedef struct bc_flat_list {
  int node_type; // Unnecessary?
  int is_single;
  int i, j, k;
  int izl, iyl, ixl;
  int izu, iyu, ixu;
  GrGeomOctree *bc_node;
  struct bc_flat_list *next;
} BCFlatList;


#define NewBCFlatList(ptr) \
  ptr = (BCFlatList*)malloc(sizeof(BCFlatList))

#define FreeBCFlatList(ptr)                    \
  while (ptr != NULL)                         \
  {                                             \
    BCFlatList *temp = ptr->next;             \
    free(ptr);                                  \
    ptr = temp;                                 \
  }

#define NewBCFlatListArray(ptr, ipatch, bc_struct)                      \
  ptr = (BCFlatList**)malloc(sizeof(BCFlatList*) * ((bc_struct)->num_patches)); \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    ptr[(ipatch)] = NULL;                                               \
  }

#define FreeBCFlatListArray(ptr, ipatch, bc_struct)    \
  {                                                    \
    BCFlatList *node;                                  \
    ForBCStructNumPatches((ipatch), bc_struct)         \
    {                                                  \
      node = ptr[(ipatch)];                            \
      FreeBCFlatList(node);                            \
    }                                                  \
    free(ptr);                                         \
  }

#define PatchType(_case, _list)                                         \
  case _case:                                                           \
  {                                                                     \
    NewBCFlatList(_list[ipatch]);                                       \
    CollectBCStructPatchLoop(i, j, k, ival, bc_struct, ipatch, is, _list[ipatch]); \
    break;                                                              \
  }

#define CollectBCFlatLists(ipatch, bc_struct, patches)  \
  ForSubgridI(is, GridSubgrids(grid))                   \
  {                                                     \
    ForBCStructNumPatches(ipatch, bc_struct)            \
    {                                                   \
      switch(BCStructBCType(bc_struct, ipatch))         \
      {                                                 \
        patches;                                        \
      }                                                 \
    }                                                   \
  }


#define CheckBCFlatListNull(ipatch, bc_stuct, list)   \
  ForSubgridI(is, GridSubgrids(grid))                 \
  {                                                   \
    ForBCStructNumPatches(ipatch, bc_struct)          \
    {                                                 \
      BCFlatList *_node_ = bc_list[ipatch];           \
      if (_node_ == NULL)                             \
        fprintf(stderr, "Node null\n");               \
      else {                                          \
        if (_node_->next == NULL)                     \
          fprintf(stderr, "Only head\n");             \
        else                                          \
          while (_node_ != NULL) {                    \
            if (_node_ == NULL) {                     \
              fprintf(stderr, "End of list\n");       \
              break;                                  \
            }                                         \
            if (_node_->next == NULL) {               \
              break;                                  \
            }                                         \
            if (_node_->bc_node == NULL)              \
              fprintf(stderr, "BCNode null\n");       \
            _node_ = _node_->next;                    \
          }                                           \
      }                                               \
    }                                                 \
  }

#endif // _BC_FLAT_H
