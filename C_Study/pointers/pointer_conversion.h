#ifndef POINTER_CONVERSION_H_INCLUDED
#define POINTER_CONVERSION_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct structA structA;

typedef struct structB structB;

typedef struct structC structC;

structA *init_structA(int data1, int data2, int data3);

structB *init_structB(char data1, char data2);

structC *init_structC(void** pointers, int cnt_pointer);

void print_structA(structA *str);

void print_structB(structB *str);

void print_structC(structC *str, char *strTypes);

void destroy_structA(structA *strA);

void destroy_structB(structB *strB);

void destroy_structC(structC *strC);

void destroy_structC_specifiedType(structC *str, char *strTypes);



#endif