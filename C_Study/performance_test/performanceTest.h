#ifndef PERFORMANCETEST_H_INCLUDED
#define PERFORMANCETEST_H_INCLUDED

#pragma once

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**
 * @brief  Test the speed of iteration between linked-list and array.
 * @param  cnt_record: count of record to iterate
 */
void test_linkedList_array(int cnt_record);

#endif