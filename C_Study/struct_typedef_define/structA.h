#ifndef STRUCTA_H_INCLUDED
#define STRUCTA_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

static int cnt_strA = 0;

typedef struct structA structA;

structA *init_structA(int data);

void modify_structA(structA *strA, int new_data);

void append_structA(structA *strA, int new_data);

void print_structA(structA *strA);

void traverse_structA(structA *strA);

void destroy_structA(structA *strA);

extern int data_structA(structA *strA);

#endif