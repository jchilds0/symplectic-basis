//
// Created by joshu on 5/09/2023.
//

#ifndef SYMPLECTIC_BASIS_H
#define SYMPLECTIC_BASIS_H

#include "kernel.h"

extern int** get_symplectic_basis(Triangulation *, int *, int *, int);
/**<
 *  Returns the symplectic basis
 */

extern void free_symplectic_basis(int **, int);
/**<
 *  Returns the symplectic basis
 */

#endif /* SYMPLECTIC_BASIS_H */
