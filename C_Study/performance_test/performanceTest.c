#include "performanceTest.h"

/**
 * @brief  Linked list without header
 */
typedef struct _def_LinkedList {
  int element;
  struct _def_LinkedList *next;
} LinkedList;

void test_linkedList_array(int cnt_record) {
  if (cnt_record <= 0) {
    fprintf(stderr, "Warning: empty records.\n");
    return;
  }
  clock_t time_start, time_end;

  int array[cnt_record];
  for (int i = 0; i < cnt_record; i++) {
    array[i] = i;
  }

  // Initialize the linked list without header
  LinkedList *list = NULL;
  LinkedList *tmp_list = (LinkedList *)calloc(1, sizeof(LinkedList));
  list = tmp_list;
  int idx_linked_list = 0;
  for (int i = 0; i < cnt_record; i++) {
    if (i < cnt_record - 1) {
      tmp_list->next = (LinkedList *)calloc(1, sizeof(LinkedList));
    }
    tmp_list->element = i;
    tmp_list = tmp_list->next;
  }

  // Iterate array
  time_start = clock();
  for (int i = 0; i < cnt_record; i++) {
    array[i] = ((array[i] + 10) | 0x212) + 3;
    for (int i = 0; i < array[i]; i++) {
    }
  }
  time_end = clock();
  printf("array time (cnt: %d): %f\n", cnt_record,
         ((float)(time_end - time_start) / CLOCKS_PER_SEC));

  // Iterate linked list
  time_start = clock();
  tmp_list = list;
  while (tmp_list != NULL) {
    tmp_list->element = ((tmp_list->element + 10) | 0x212) + 3;
    for (int i = 0; i < tmp_list->element; i++) {
    }
    tmp_list = tmp_list->next;
  }
  time_end = clock();
  printf("linked list time (cnt: %d): %f\n", cnt_record,
         ((float)(time_end - time_start) / CLOCKS_PER_SEC));

  // Destroy the linked list without header
  tmp_list = list;
  while (tmp_list != NULL) {
    LinkedList *next = tmp_list->next;
    free(tmp_list);
    tmp_list = next;
  }
}
