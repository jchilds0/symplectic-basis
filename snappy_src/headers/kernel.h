/**
 *  @file kernel.h
 *  @brief Includes all public and private headers files, for building the kernel.
 *
 *  This file #includes all header files needed for the kernel.
 *  It should be #included in all kernel .c files, but nowhere else.
 */

#ifndef _kernel_
#define _kernel_

#include "SnapPea.h"

#ifdef FORCE_C_LINKAGE
#ifdef __cplusplus
extern "C" {
#endif
#endif

#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "kernel_typedefs.h"
#include "triangulation.h"
#include "positioned_tet.h"
#include "kernel_prototypes.h"
#include "tables.h"

#ifdef FORCE_C_LINKAGE
#ifdef __cplusplus
}
#endif
#endif

#endif
/* Local Variables:                      */
/* mode: c                               */
/* c-basic-offset: 4                     */
/* fill-column: 80                       */
/* comment-column: 0                     */
/* c-file-offsets: ((inextern-lang . 0)) */
/* End:                                  */
