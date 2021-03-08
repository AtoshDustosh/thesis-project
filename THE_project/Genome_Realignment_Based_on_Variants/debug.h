#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#pragma once

// #define NDEBUG 0

/*
 * The following statement can make the program automatically include <assert.h>
 * when debugging and activate asssert(). It can also automatically disable
 * assert() when not debugging, saving space after  compilation.
 */
#ifdef NDEBUG
#define assert(statement) ((void)0)
#else
#include <assert.h>
#endif

#endif