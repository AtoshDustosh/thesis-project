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
 * Every field except the vts is just a reference to another object and thus
 * doesn't need to be freed.
 */
typedef struct _define_VarTodoChrom {
  char *name;
  uint32_t varCnt;
  // indexes of variants to be integrated
  VarTodo *vts;
  struct _define_VarTodoChrom *next;
} VarTodoChrom;

typedef struct _define_VarIntegration {
  /*
   * note that this gv is not a copy of other GenomeVcf. It's only a reference
   * to another gv object. Thus, you can but must not destroy the gv object
   * using destroy_GenomeVcf() when you want to destroy this VarIntegration
   * Object.
   */
  GenomeVcf *gv;
  // vtcs is a linked-list without header. We don't need to sort the
  // chromosomes, so this is not implemented with a header.
  uint32_t chromCnt;
  VarTodoChrom *vtcs;
} VarIntegration;

typedef struct _define_VarIntegrationIterator {
  VarIntegration *vi;
  VarTodoChrom *tmpVtc;
  VarTodo *tmpVt;
} VarIntegrationIterator;

static inline uint32_t vtData(VarTodo *vt) { return vt->varIdx; }

static inline char *vtcName(VarTodoChrom *vtc) { return vtc->name; }

VarIntegrationIterator *init_VarIntegrationIterator(VarIntegration *vi);

/**
 * @brief  This iterator return the next chrom to be iterated.
 * @retval pointer to the next VarTodoChrom object to be iterated; NULL if there
 * is no chroms left or the iterator is not initialized with a non-NULL
 * VarIntegration.
 */
VarTodoChrom *viItNextChrom(VarIntegrationIterator *viIt);

/**
 * @brief  This iterator return the next variant to be iterated.
 * @retval pointer to the next VarTodo object to be i terated. NULL if there is
 * no variants left in the temporary chrom, chrom is not selected for iteration
 * (use @viItNextChrom before using this), or the iterator is not intialized
 * with a non-NULL VarIntegration.
 */
VarTodo *viItNextVar(VarIntegrationIterator *viIt);

void destroy_VarIntegrationIterator(VarIntegrationIterator *viIt);

VarTodo *init_VarTodo();

void destroy_VarTodo(VarTodo *vt);

VarTodoChrom *init_VarTodoChrom();

void destroy_VarTodoChrom(VarTodoChrom *vtc);

VarIntegration *init_VarIntegration(GenomeVcf *gv);

void destroy_VarIntegration(VarIntegration *vi);

/*****************************************
 * Methods for manipulating VarIntegration
 ****************************************/

/**
 * @brief  Add a VarTodo object into the VarTodoChrom object.
 * @note   This actually serves as a function for marking which variants should
 * be integrated later.
 */
void addVtToVtc(uint32_t varIdx, VarTodoChrom *vtc);

void addVtcToVarInt(VarTodoChrom *vtc, VarIntegration *vi);

/**
 * @brief  Get the (vtIdx+1)-th VarTodo object in vtc.
 * @param  vtIdx: 0-based index/id for the VarTodo object
 */
VarTodo *getVtFromVtc(uint32_t vtIdx, VarTodoChrom *vtc);

VarTodoChrom *getVtcFromVarInt(char *chromName, VarIntegration *vi);

/*************************************
 * Debugging Methods for VarIntegration
 ************************************/
void _testSet_varIntegration();

void printVarIntegration(VarIntegration *vi);

#endif