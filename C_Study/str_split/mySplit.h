#ifndef MYSPLIT_H_INCLUDED
#define MYSPLIT_H_INCLUDED

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void _testSet_mySplit();

/**
 * @brief Split the string with delim. This method treat delimiter as a whole
 * string. For example, ">abc - def.gh-i -> j" with delim " - " will be split
 * into {">abs", "def.gh-i -> j"}. And this method will not produce any empty
 * strings. For example, ">abs - def.gh-i-j" with delim ">" will be split into
 * {"abs - def.gh-i -", " j"}. NOTE that "tokens" must be manually freed after
 * usage, otherwise may results in a memory leakage.
 * @todo note that this method has bugs. When handling strings "abcddef" with
 * delim "def", it may not recognize the delim "def". And pleas do not use this
 * method for any frequent executed methods. Planning to fix the problem in the
 * future.
 *
 * @param str string to be split
 * @param delim delimiter (treated as a whole)
 * @param tokens tokens split
 * @return int number of tokens
 */
int mySplit(char *str, char *delim, char ***tokens);

#endif