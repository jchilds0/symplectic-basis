/*
 * symplectic_basis.h
 *
 * Exports the functions
 *
 *      int **get_symplectic_basis(Triangulation *, int *, int *, int);
 *
 *      void free_symplectic_basis(int **, int);
 *
 * to the user. The files in symp_src are written with the intention of
 * being collected into one file symplectic_basis.c and included in
 * the SnapPy kernel. The program is split into multiple files currently
 * for development and due to uncertainty around the existence of
 * oscillating curves in general.
 */


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
