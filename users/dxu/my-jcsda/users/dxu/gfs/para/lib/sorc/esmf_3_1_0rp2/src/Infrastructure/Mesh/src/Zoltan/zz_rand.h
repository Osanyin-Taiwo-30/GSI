/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: zz_rand.h,v $
 *    $Author: dneckels $
 *    $Date: 2007/08/08 22:44:00 $
 *    Revision: 1.5 $
 ****************************************************************************/
#ifndef __ZZ_RAND_H
#define __ZZ_RAND_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define ZOLTAN_RAND_MAX 0xffffffff
#define ZOLTAN_RAND_INIT 123456789U

extern unsigned int Zoltan_Seed();
extern unsigned int Zoltan_Rand(unsigned int *);
extern unsigned int Zoltan_Rand_InRange(unsigned int *, unsigned int);
extern void Zoltan_Srand(unsigned int, unsigned int *);
extern void Zoltan_Rand_Perm_Int(int*, int, unsigned int *);
extern void Zoltan_Srand_Sync(unsigned int, unsigned int *, MPI_Comm);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZZ_RAND_H */
