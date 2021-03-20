#ifndef GENOMEFAMACROS_H_INCLUDED
#define GENOMEFAMACROS_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>

#include "debug.h"

// BASE_NORECORD is used to pad the coded bases if there is no left bases
#define BASE_NORECORD 0b000

#define BASE_A 0b001
#define BASE_C 0b010
#define BASE_G 0b011
#define BASE_T 0b100
#define BASE_N 0b101

#define BASE_INVALID 0b111

// this BASE_MASK_LEFT and BASE_MASK_RIGHT can be used to extract a coded base
// from a long integer (0b) 111 000 000 ... 000 0
#define BASE_MASK_LEFT 0xE000000000000000
#define BASE_MASK_RIGHT 0x0000000000000007

#define BASE_CODE_LENGTH 3

#define MAX_INFO_LENGTH 2048

#define BP_PER_UINT64 (sizeof(uint64_t) * 8 / 3)

#define BP_PER_LINE 70

typedef uint8_t Base;

char charOfBase(Base bp);

Base baseOfChar(char bp);

inline static void assert_bases_equal(Base x, Base y) {
  // note that the "||" operation can terminate the condition statement
  // instantly when "x == y" succeeds and thus printf will not print. (a great
  // way for debugging using "assert(...)")
  assert(x == y || (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                            x, y) >= 0));
}

#endif