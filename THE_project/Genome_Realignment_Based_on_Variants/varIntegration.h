#ifndef VARINTEGRATION_H_INCLUDED
#define VARINTEGRATION_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "debug.h"
#include "genomeVcf.h"

/******************
 * Basic Structures
 ******************/

/**
 * @brief  An index about a variant to be integrated into the reference genome.
 */
typedef struct _define_VarTodo {
  uint32_t varIdx;
  struct _define_VarTodo *next;
} VarTodo;

/**
 * @note vts is actually a linked-list with a header. This structure is easier
 * to sort by field values than a linked-list without a header.
 */
typedef struct _define_VarTodoChrom {
  char *name;
  uint32_t varCnt;
  // indexes of variants to be integrated
  VarTodo *vts;
  struct _define_VarTodoChrom *next;
} VarTodoChrom;

typedef struct _define_VarIntegration {
  // note that this gv is not a copy of other GenomeVcf. It's only a reference
  // to another gv object. Thus, you can but must not destroy the gv object
  // using destroy_GenomeVcf() when you want to destroy this VarIntegration
  // Object.
  GenomeVcf *gv;
  // vtcs is a linked-list without header. We don't need to sort the
  // chromosomes, so this is not implemented with a header.
  uint32_t chromCnt;
  VarTodoChrom *vtcs;
} VarIntegration;

static VarTodo *init_VarTodo();

static void destroy_VarTodo(VarTodo *vt);

static VarTodoChrom *init_VarTodoChrom();

static void destroy_VarTodoChrom(VarTodoChrom *vtc);

VarIntegration *init_VarIntegration();

void destroy_VarIntegration(VarIntegration *vi);

/*************************************
 * Debugging Methods for VarIntegration
 ************************************/
void _testSet_varIntegration();

void printVarIntegration(VarIntegration *vi);

/*****************************************
 * Methods for manipulating VarIntegration
 ****************************************/

void addVtToVtc(VarTodo *vt, VarTodoChrom *vtc);

void addVtcToVarInt(VarTodoChrom *vtc, VarIntegration *vi);

VarTodo getVtFromVtc(uint32_t vtIdx, VarTodoChrom *vtc);

VarTodoChrom getVtcFromVarInt(char *chromName, VarIntegration *vi);

void markVar(char *chromName, uint32_t varIdx);

/*******************
 * Methods for users
 *******************/

#endif