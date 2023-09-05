/**
 * @mainpage Documentation
 * \section Overview
 * This space should contain an overview of the kernel code, written in markdown,
 * which will appear on the main page of the html documentation.
 */

/**
 *  @file SnapPea.h
 *
 *  @brief The public interface to the SnapPea kernel.
 *
 *  This file defines the interface between SnapPea's computational kernel
 *  ("the kernel") and the user-interface ("the UI").  Both parts
 *  must \#include this file, and anything shared between the two parts
 *  must be declared in this file.  The only communication between the
 *  two parts is via function calls -- no external variables are shared.
 *
 *  All external symbols in the UI must begin with 'u' followed by a
 *  capital letter.  Nothing in the kernel should begin in this way.
 *
 *  Typedef names use capitals for the first letter of each word,
 *  e.g. Triangulation, CuspIndex.
 *
 *  SnapPea 2.0 was funded by the University of Minnesota's
 *  Geometry Center and the U.S. National Science Foundation.
 *  SnapPea 3.0 is funded by the U.S. National Science Foundation
 *  and the MacArthur Foundation.  SnapPea and its source code may
 *  be used freely for all noncommercial purposes.  Please direct
 *  questions, problems and suggestions to Jeff Weeks (weeks@northnet.org).
 *
 *  Copyright 1999 by Jeff Weeks.  All rights reserved.
 *  Copyright 2008-present by Marc Culler, Nathan Dunfield, Matthias Goerner.
 */

#ifndef _SnapPea_
#define _SnapPea_

#include "real_type.h"

#include "kernel_namespace.h"

/**
 *  Note:  values of the SolutionType enum are stored as integers in
 *  the triangulation.doc file format.  Changing the order of the
 *  entries in the enum would therefore invalidate all previously stored
 *  triangulations.
 */

typedef enum
{
    not_attempted,          /**<  solution not attempted, or user cancelled                   */
    geometric_solution,     /**<  all positively oriented tetrahedra; not flat or degenerate  */
    nongeometric_solution,  /**<  positive volume, but some negatively oriented tetrahedra    */
    flat_solution,          /**<  all tetrahedra flat, but no shapes = {0, 1, infinity}       */
    degenerate_solution,    /**<  at least one tetrahedron has shape = {0, 1, infinity}       */
    other_solution,         /**<  volume <= 0, but not flat or degenerate                     */
    no_solution,            /**<  gluing equations could not be solved                        */
    externally_computed     /**<  tetrahedra shapes were inserted into the triangulation      */
} SolutionType;

/**
 *  The constants complete and filled facilitate reference
 *  to the shape of a Tetrahedron as part of the complete or
 *  Dehn filled hyperbolic structure, respectively.
 */

typedef enum
{
    complete,
    filled
} FillingStatus;

typedef enum
{
    func_OK = 0,
    func_cancelled,
    func_failed,
    func_bad_input
} FuncResult;;

typedef struct
{
    Real  real;
    Real  imag;
} Complex;

typedef char   Boolean;

/**
 *  The values of MatrixParity should not be changed.
 *  (They must correspond to the values in the parity[] table in tables.c.)
 */

typedef enum
{
    orientation_reversing = 0,
    orientation_preserving = 1,
    unknown_parity = -1
} MatrixParity;

/**
 *  SnapPea represents a Moebius transformation as a matrix
 *  in SL(2,C) plus a specification of whether the Moebius
 *  transformation is orientation_preserving or orientation_reversing.
 *
 *  If mt->parity is orientation_preserving, then mt->matrix is
 *  interpreted in the usual way as the Moebius transformation
 *
 *                      az + b
 *              f(z) = --------
 *                      cz + d
 *
 *
 *  If mt->parity is orientation_reversing, then mt->matrix is
 *  interpreted as a function of the complex conjugate z' ("z-bar")
 *
 *                      az' + b
 *              f(z) = ---------
 *                      cz' + d
 */

typedef Complex SL2CMatrix[2][2];

typedef struct
{
    SL2CMatrix      matrix;
    MatrixParity    parity;
} MoebiusTransformation;

/**
 *  Matrices in O(3,1) represent isometries in the Minkowski space
 *  model of hyperbolic 3-space.  The matrices are expressed relative
 *  to a coordinate system in which the metric is
 *
 *                          -1  0  0  0
 *                           0  1  0  0
 *                           0  0  1  0
 *                           0  0  0  1
 *
 *  That is, the first coordinate is timelike, and the remaining
 *  three are spacelike.  O(3,1) matrices represent both
 *  orientation_preserving and orientation_reversing isometries.
 */

typedef Real O31Matrix[4][4];
typedef Real GL4RMatrix[4][4];

/**
 *  An O31Vector is a vector in (3,1)-dimensional Minkowski space.
 *  The 0-th coordinate is the timelike one.
 */

typedef Real O31Vector[4];

/**
 *  MatrixInt22 is a 2 x 2 integer matrix.  A MatrixInt22
 *  may, for example, describe how the peripheral curves of
 *  one Cusp map to those of another.
 */

typedef int MatrixInt22[2][2];

/**
 *  An AbelianGroup is represented as a sequence of torsion coefficients.
 *  A torsion coefficient of 0 represents an infinite cyclic factor.
 *  For example, the group Z + Z + Z/2 + Z/5 is represented as the
 *  sequence (0, 0, 2, 5).  We make the convention that torsion coefficients
 *  are always nonnegative.
 *
 *  The UI may declare pointers to AbelianGroups, but only the kernel
 *  may allocate or deallocate the actual memory used to store an
 *  AbelianGroup.  (This allows the kernel to keep track of memory
 *  allocation/deallocation as a debugging aid.)
 */

typedef struct
{
    int         num_torsion_coefficients;   /**<  number of torsion coefficients              */
    long int    *torsion_coefficients;      /**<  pointer to array of torsion coefficients    */
} AbelianGroup;

/**
 *  A closed geodesic may be topologically a circle or a mirrored interval.
 */

typedef enum
{
    orbifold1_unknown,
    orbifold_s1,    /**<  circle              */
    orbifold_mI     /**<  mirrored interval   */
} Orbifold1;

/**
 *  The following 2-orbifolds may occur as the link of an
 *  edge midpoint in a cell decomposition of a 3-orbifold.
 *
 *  94/10/4.  The UI will see only types orbifold_nn
 *  and orbifold_xnn.  Edges of the other types have 0-cells
 *  of the singular set at their midpoints, and are now
 *  subdivided in Dirichlet_extras.c.  JRW
 */

typedef enum
{
    orbifold_nn,    /**<  (nn)    2-sphere with two cone points (n may be 1)  */
    orbifold_no,    /**<  (n|o)   cross surface with cone point (n may be 1)  */
    orbifold_xnn,   /**<  (*nn)   disk with mirror boundary with two
                     *           corner reflectors
                     */
    orbifold_2xn,   /**<  (2*n)   disk with order two cone point and mirror
                     *           boundary with one corner reflector
                     */
    orbifold_22n    /**<  (22n)   sphere with three cone points               */
} Orbifold2;

/**
 *  A MultiLength records the complex length of a geodesic together with a
 *  parity telling whether it preserves or reverses orientation, a topology
 *  telling whether it's a circle or a mirrored interval, and a multiplicity
 *  telling how many distinct geodesics have that complex length, parity and
 *  topology.  Finally, the matrix of some fundamental group element (and,
 *  optionally, the corresponding word in the original generators) realizing
 *  this geodesic is given.
 */

typedef struct
{
    Complex         length;
    MatrixParity    parity;
    Orbifold1       topology;
    int             multiplicity;
    O31Matrix       matrix;
    int             *word;     /* Added MG 2022-01-24 */
} MultiLength;


/**
 *  A CuspNbhdHoroball records a horoball to be drawn as part of a
 *  picture of a cusp cross section.  Only the kernel should allocate
 *  and free CuspNbhdHoroballs and CuspNbhdHoroballLists.  These
 *  definitions are provided to the UI so it access the data easily.
 */

typedef struct
{
    Complex center;
    Real  radius;
    int     cusp_index;
} CuspNbhdHoroball;

typedef struct
{
    int                 num_horoballs;
    /**
     *  The horoball field points to an array
     *  of num_horoballs CuspNbhdHoroballs.
     */
    CuspNbhdHoroball    *horoball;
} CuspNbhdHoroballList;


/**
 *  A CuspNbhdSegment records a 1-cell to be drawn as part of a
 *  picture of a cusp cross section.  (Typically it's either part of
 *  a triangulation of the cusp cross section, or part of a Ford domain.)
 *  Only the kernel should allocate and free CuspNbhdSegments and
 *  CuspNbhdSegmentLists.  These definitions are provided to the UI
 *  so it can easily access the data.
 *
 *  JRW 99/03/17   When the CuspNbhdSegment describes a triangulation
 *  (as opposed to a Ford domain),
 *
 *      the start_index tells the edge index of the vertical edge
 *          that runs from the given segment's beginning
 *          to the viewer's eye,
 *
 *      the middle_index tells the edge index of the given segment, and
 *
 *      the end_index tells the edge index of the vertical edge
 *          that runs from the given segment's end
 *          to the viewer's eye.
 *
 *  These indices let the viewer see how the horoball picture
 *  "connects up" to form the manifold.
 */

typedef struct
{
    Complex endpoint[2];
    int     start_index,
            middle_index,
            end_index;
} CuspNbhdSegment;

typedef struct
{
    int             num_segments;
    /**
     *  segment is a pointer to an array of num_segments CuspNbhdSegments.
     */
    CuspNbhdSegment *segment;
} CuspNbhdSegmentList;


typedef enum
{
    oriented_manifold,
    nonorientable_manifold,
    unknown_orientability
} Orientability;


typedef enum
{
    torus_cusp,
    Klein_cusp,
    unknown_topology
} CuspTopology;


typedef enum
{
    Dirichlet_interactive,
    Dirichlet_stop_here,
    Dirichlet_keep_going
} DirichletInteractivity;


/**
 *  An LRFactorization specifies the monodromy for a punctured torus
 *  bundle over a circle.  The factorization is_available whenever
 *  (det(monodromy) = +1 and |trace(monodromy)| >= 2) or
 *  (det(monodromy) = -1 and |trace(monodromy)| >= 1).
 *  LR_factors points to an array of L's and R's, interpreted as factors
 *
 *              L = ( 1  0 )            R = ( 1  1 )
 *                  ( 1  1 )                ( 0  1 )
 *
 *  The factors act on a column vector, beginning with the last
 *  (i.e. rightmost) factor.
 *
 *  If negative_determinant is TRUE, the product is left-multiplied by
 *
 *                          ( 0  1 )
 *                          ( 1  0 )
 *
 *  If negative_trace is TRUE, the product is left-multiplied by
 *
 *                          (-1  0 )
 *                          ( 0 -1 )
 *
 *  When the factorization is unavailable, is_available is set to FALSE,
 *  num_LR_factors is set to zero, and LR_factors is set to NULL.
 *  But the negative_determinant and negative_trace flags are still set,
 *  so the UI can display this information correctly.
 */
typedef struct
{
    Boolean is_available;
    Boolean negative_determinant;
    Boolean negative_trace;
    int     num_LR_factors;
    char    *LR_factors;
} LRFactorization;


/**
 *  The full definition of a Shingling appears near the top of shingling.c.
 *  But computationally a Shingling is just a collection of planes in
 *  hyperbolic space (typically viewed as circles on the sphere at infinity).
 *  Each plane has an index (which defines the color of the circle at
 *  infinity).
 */

typedef struct
{
    /**
     *  A plane in hyperbolic 3-space defines a hyperplane through
     *  the origin in the Minkowski space model.  Use the hyperplane's
     *  normal vector to represent the original plane.  [Note:  the
     *  normal is computed once, in the standard coordinate system,
     *  and does not change as the UI rotates the polyhedron.]
     */
    O31Vector   normal;

    /**
     *  A plane in hyperbolic 3-space intersects the sphere at infinity
     *  in a circle.  It's easy to draw the circle if we know its center
     *  and two orthogonal "radials".  (The 0-components of the center
     *  and radials may be ignored.)  [Note:  the center and radials are
     *  rotated in real time according to the polyhedron's current
     *  position, and are scaled according to the window's pixel size.]
     */
    O31Vector   center;
    O31Vector   radialA;
    O31Vector   radialB;

    /**
     *  The face planes of the original Dirichlet domain have index 0,
     *  the face planes of the next layer (cf. shingling.c) have index 1,
     *  and so on.
     */
    int         index;

} Shingle;

/**
 *  A Shingling is just an array of Shingles.
 */
typedef struct
{
    int         num_shingles;
    Shingle     *shingles;

} Shingling;


/*
 *  The following are "opaque typedefs".  They let the UI declare and
 *  pass pointers to Triangulations, IsometryLists, etc. without
 *  knowing what a Triangulation, IsometryList, etc. is.  The definitions
 *  of struct Triangulation, struct IsometryList, etc. are private to the
 *  kernel.  SymmetryLists and IsometryLists are represented by the same
 *  data structure because Symmetries are just Isometries from a manifold
 *  to itself.
 */

typedef struct Triangulation                Triangulation;
typedef struct IsometryList                 IsometryList;
typedef struct SymmetryGroup                SymmetryGroup;
typedef struct SymmetryGroupPresentation    SymmetryGroupPresentation;
typedef struct DualOneSkeletonCurve         DualOneSkeletonCurve;
typedef struct TerseTriangulation           TerseTriangulation;
typedef struct GroupPresentation            GroupPresentation;
typedef struct CuspNeighborhoods            CuspNeighborhoods;
typedef struct NormalSurfaceList            NormalSurfaceList;

#include "end_namespace.h"

/*
 *  When the UI reads a Triangulation from disk, it passes the results
 *  to the kernel using the format described in triangulation_io.h.
 */
#include "triangulation_io.h"

/*  To guarantee thread-safety, it's useful to declare      */
/*  global variables to be "const", for example             */
/*                                                          */
/*      static const Complex    minus_i = {0.0, -1.0};      */
/*                                                          */
/*  Unfortunately the current gcc compiler complains when   */
/*  non-const variables are passed to functions expecting   */
/*  const arguments.  Obviously this is harmless, but gcc   */
/*  complains anyhow.  So for now let's use the following   */
/*  CONST macro, to allow the const declarations to be      */
/*  reactivated if desired.                                 */
/*                                                          */
/*  Note:  In Win32, windef.h also defines CONST = const,   */
/*  so clear its definition before making our own.          */
#undef  CONST
/**
 *   Placeholder for currently unused const declarations.
 */
#define CONST
/* #define CONST const */

#ifdef FORCE_C_LINKAGE
#ifdef __cplusplus
extern "C" {
#endif
#endif

/************************************************************************/

#include "kernel_namespace.h"

/*
 *  The UI provides the following functions for use by the kernel:
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void uAcknowledge(const char *message);
/**<
 *  Presents the string *message to the user and waits for acknowledgment ("OK").
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
int uQuery(const char *message, const int num_responses,
                const char *responses[], const int default_response);
/**<
 *  Presents the string *message to the user and asks the user to choose
 *  one of the responses.  Returns the number of the chosen response
 *  (numbering starts at 0).  In an interactive context, the UI should
 *  present the possible responses evenhandedly -- none should be
 *  presented as a default.  However, in a batch context (when no human
 *  is present), uQuery should return the default_response.
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void uFatalError(const char *function, const char *file);
/**<
 *  Informs the user that a fatal error has occurred in the given
 *  function and file, and then exits.
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void uAbortMemoryFull(void);
/**<
 *  Informs the user that the available memory has been exhausted,
 *  and aborts SnapPea.
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void uPrepareMemFullMessage(void);
/**<
 *  uMemoryFull() is a tricky function, because the system may not find
 *  enough memory to display an error message.  (I tried having it stash
 *  away some memory and then free it to support the desired dialog box,
 *  but at least on the Mac this didn't work for some unknown reason.)
 *  uPrepareMemFullMessage() gives the system a chance to prepare
 *  a (hidden) dialog box.  Call it once when the UI initializes.
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void         uLongComputationBegins(const char *message,
                                    Boolean is_abortable);
/**<
 *  The kernel uses the three functions uLongComputationBegins(),
 *  uLongComputationContinues() and uLongComputationEnds() to inform
 *  the UI of a long computation.  The UI relays this information to
 *  the user in whatever manner it considers appropriate.  For
 *  example, it might wait a second or two after the beginning of a
 *  long computation, and then display a dialog box containing
 *  *message (a typical message might be "finding canonical
 *  triangulation" or "computing hyperbolic structure").  If
 *  is_abortable is TRUE, the dialog box would contain an abort
 *  button.  The reason for waiting a second or two before displaying
 *  the dialog box is to avoid annoying the user with flashing dialog
 *  boxes for computations which turn out not to be so long after all.
 *
 *  The kernel is responsible for calling uLongComputationContinues() at
 *  least every 1/60 second or so during a long computation.
 *  uLongComputationContinues() serves two purposes:
 *
 *  (1) It lets the UI yield time to its window system.  (This is
 *      crucial for smooth background operation in the Mac's
 *      cooperative multitasking environment.  I don't know whether
 *      it is necessary in X or NeXT.)
 *
 *  (2) If the computation is abortable, it checks whether the user
 *      has asked to abort, and returns the result (func_cancelled
 *      to abort, func_OK to continue).
 *
 *  While the kernel is responsible for making sure uLongComputationContinues()
 *  is called often enough, uLongComputationContinues() itself must take
 *  responsibility for not repeating time-consuming operations too often.
 *  For example, it might return immediately from a call if less than
 *  1/60 of a second has elapsed since the last time it carried out
 *  its full duties.
 *
 *  uLongComputationEnds() signals that the long computation is over.
 *  The kernel must call uLongComputationEnds() even after an aborted
 *  computation.  ( uLongComputationContinues() merely informs the kernel
 *  that the user punched the abort button.  The kernel must still call
 *  uLongComputationEnds() to dismiss the dialog box in the usual way.)
 *
 *  If the UI receives a call to uLongComputationEnds() when no long
 *  computation is in progress, or a call to uLongComputationBegins()
 *  when a long computation is already in progress, it should notify
 *  the user of the error and exit.
 *
 *  If the UI receives a call to uLongComputationContinues() when in
 *  fact no long computation is in progress, it should simply take
 *  care of any background responsibilities (see (1) above) and not
 *  complain.  The reason for this provision is that the calls to
 *  uLongComputationBegins() and uLongComputationEnds() occur in high
 *  level functions, while the calls to uLongComputationContinues()
 *  occur at the lowest level, perhaps in a different file.  Someday
 *  those low-level functions (for example, the routines for solving
 *  simultaneous linear equations) might be called as part of some quick,
 *  non-abortable computation.
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
FuncResult   uLongComputationContinues(void);
/**<
 * See uLongComputationBegins().
 */

extern
#ifdef _MSC_VER
 __declspec(dllexport)
#endif
void         uLongComputationEnds(void);
/**<
 * See uLongComputationBegins().
 */


/************************************************************************/
/************************************************************************/


/************************************************************************/
/*                                                                      */
/*                              chern_simons.c                          */
/*                                                                      */
/************************************************************************/

extern void set_CS_value(   Triangulation   *manifold,
                            Real            a_value);
/**<
 *  Set the Chern-Simons invariant of *manifold to a_value.
 */

extern void get_CS_value(   Triangulation   *manifold,
                            Boolean         *value_is_known,
                            Real            *the_value,
                            int             *the_precision,
                            Boolean         *requires_initialization);
/**<
 *  If the Chern-Simons invariant of *manifold is known, sets
 *  *value_is_known to TRUE and writes the current value and its precision
 *  (the number of significant digits to the right of the decimal point)
 *  to *the_value and *the_precision, respectively.
 *
 *  If the Chern-Simons invariant is not known, sets *value_is_known to
 *  FALSE, and then sets *requires_initialization to TRUE if the_value
 *  is unknown because the computation has not been initialized, or
 *  to FALSE if the_value is unknown because the solution contains
 *  negatively oriented Tetrahedra.  The UI might want to convey
 *  these situations to the user in different ways.
 */


/************************************************************************/
/*                                                                      */
/*                              complex.c                               */
/*                                                                      */
/************************************************************************/

/**
 *
 *  Complex arithmetic operator.
 *  Standard complex constants (Zero, One, etc.) are defined in the kernel.
 */
/** @{ */
extern Complex  complex_minus           (Complex z0, Complex z1),
                complex_plus            (Complex z0, Complex z1),
                complex_mult            (Complex z0, Complex z1),
                complex_div             (Complex z0, Complex z1),
                complex_sqrt            (Complex z),
                complex_conjugate       (Complex z),
                complex_negate          (Complex z),
                complex_real_mult       (Real r, Complex z),
                complex_exp             (Complex z),
                complex_log             (Complex z, Real approx_arg);
extern Real     complex_modulus         (Complex z);
extern Boolean  complex_infinite        (Complex z);
extern Real     complex_modulus_squared (Complex z);
extern Boolean  complex_nonzero         (Complex z);
/** @} */


/************************************************************************/
/*                                                                      */
/*                              core_geodesic.c                         */
/*                                                                      */
/************************************************************************/

extern void core_geodesic(  Triangulation   *manifold,
                            int             cusp_index,
                            int             *singularity_index,
                            Complex         *core_length,
                            int             *precision);
/**<
 *  Examines the Cusp of index cusp_index in *manifold.
 *
 *  If the Cusp is unfilled or the Dehn filling coefficients are not
 *  integers, sets *singularity_index to zero and leaves *core_length
 *  undefined.
 *
 *  If the Cusp has relatively prime integer Dehn filling coefficients,
 *  sets *singularity_index to 1 and *core_length to the complex length
 *  of the central geodesic.
 *
 *  If the Cusp has non relatively prime integer Dehn filling coefficients,
 *  sets *singularity_index to the index of the singular locus, and
 *  *core_length to the complex length of the central geodesic in the
 *  smallest manifold cover of a neighborhood of the singular set.
 *
 *  In the latter two cases, if the precision pointer is not NULL,
 *  *precision is set to the number of decimal places of accuracy in
 *  the computed value of *core_length.
 *
 *  core_geodesic() is intended for use by the UI.  Kernel function may
 *  find compute_core_geodesic() (declared in kernel_prototypes.h) more
 *  convenient.
 */


/************************************************************************/
/*                                                                      */
/*                              my_malloc.c                             */
/*                                                                      */
/************************************************************************/

extern void verify_my_malloc_usage(void);
/**<
 *  The UI should call verify_my_malloc_usage() upon exit to verify that
 *  the number of calls to my_malloc() was exactly balanced by the number
 *  of calls to my_free().  In case of error, verify_my_malloc_usage()
 *  passes an appropriate message to uAcknowledge.
 */

/************************************************************************/
/*                                                                      */
/*                          o31_matrices.c                              */
/*                                                                      */
/************************************************************************/

/*
 *  Most of the functions in o31_matrices.c are private to the kernel.
 *  The following have been made available to the UI as well.
 */

extern Real       gl4R_determinant(GL4RMatrix m);
/**<
 * Returns the determininant of a 4 x 4 matrix.
 */

extern Real       o31_trace(O31Matrix m);
/**<
 * Returns the trace of an O(3,1) matrix.
 */

/************************************************************************/
/*                                                                      */
/*                              orient.c                                */
/*                                                                      */
/************************************************************************/

extern void reorient(Triangulation *manifold);
/**<
 *  Reverse a manifold's orientation.
 */


/************************************************************************/
/*                                                                      */
/*                          triangulations.c                            */
/*                                                                      */
/************************************************************************/

extern void data_to_triangulation(  TriangulationData   *data,
                                    Triangulation       **manifold_ptr);
/**<
 *  Uses the TriangulationData (defined in triangulation_io.h) to
 *  construct a Triangulation.  Sets *manifold_ptr to point to the
 *  Triangulation, or to NULL if it fails.
 */

extern void triangulation_to_data(  Triangulation       *manifold,
                                    TriangulationData   **data_ptr);
/**<
 *  Allocates the TriangulationData and writes in the data describing
 *  the manifold.  Sets *data_ptr to point to the result.  The UI
 *  should call free_triangulation_data() when it's done with the
 *  TriangulationData.
 */

extern void free_triangulation_data(TriangulationData *data);
/**<
 *  If the UI lets the kernel allocate a TriangulationData structure
 *      (as in a call to triangulation_to_data()), then the UI should
 *      call free_triangulation_data() to release it.
 *  If the UI allocates its own TriangulationData structure (as in
 *      preparing for a call to data_to_triangulation()), then the UI
 *      should release the structure itself.
 */

extern void free_triangulation(Triangulation *manifold);
/**<
 *  If manifold != NULL, frees up the storage associated with a
 *      triangulation structure.
 *  If manifold == NULL, does nothing.
 */

extern void copy_triangulation(Triangulation *source, Triangulation **destination);
/**<
 *  Makes a copy of the Triangulation *source.
 */


/************************************************************************/
/*                                                                      */
/*                              volume.c                                */
/*                                                                      */
/************************************************************************/

extern Real volume(Triangulation *manifold, int *precision);
/**<
 *  Computes and returns the volume of the manifold.
 *  If the pointer "precision" is not NULL, estimates the number
 *  of decimal places of accuracy, and places the result in the
 *  variable *precision.
 */

#include "end_namespace.h"

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
