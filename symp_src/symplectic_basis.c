/**
 *  Symplectic Basis
 *
 *  Computes a symplectic basis of a triangulated knot or link exterior with
 *  orientable torus cusps. This symplectic matrix extends the Neumann-Zagier
 *  matrix to one which is symplectic up to factors of 2, and which arises
 *  from the triangulation of the manifold.
 *
 *  See - https://arxiv.org/abs/2208.06969
 *
 */

#include "symplectic_kernel.h"

int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     *combinatorial_holonomy(Triangulation *, int);
void                    label_triangulation_edges(Triangulation *);

/*
 * Allocates arrays for symplectic basis and gluing equations.
 * get_gluing_equations find oscillating curves on the manifold.
 * Constructs return array using gluing_equations_for_edge_class
 * and combinatorial_holonomy
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols, int log) {
    int i, j, k;
    debug = log;
    Boolean *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, Boolean);
    Tetrahedron *tet;

    label_triangulation_edges(manifold);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        if (tet->extra != NULL)
            uFatalError("do_oscillating_curves", "symplectic_basis");

        tet->extra = NEW_ARRAY(manifold->num_tetrahedra, Extra);

        for (i = 0; i < manifold->num_tetrahedra; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    tet->extra[i].curve[j][k] = 0;
    }

    // Dual Edge Curves Gamma_i -> symplectic equations
    do_oscillating_curves(manifold, edge_classes);

    // Construct return array
    *num_rows = 2 * (manifold->num_tetrahedra - manifold->num_cusps);
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i]) {
            continue;
        }

        eqns[2 * j]     = gluing_equations_for_edge_class(manifold, i);
        eqns[2 * j + 1] = combinatorial_holonomy(manifold, i);
        j++;
    }

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        my_free(tet->extra);
        tet->extra = NULL;
    }
    my_free(edge_classes);

    *num_cols = 3 * manifold->num_tetrahedra;
    return eqns;
}

/*
 * Copy of get_gluings_equations.c get_gluing_equations() which finds 
 * the edge gluings equations for a given edge index. Used instead 
 * of get_gluing_equations to ensure we have the correct edge index 
 * and simplify memory management since we don't need all the rows of 
 * the gluing equations matrix.
 */

int *gluing_equations_for_edge_class(Triangulation *manifold, int edgeClass) {
    int *eqns, i, T;
    EdgeClass *edge;
    PositionedTet ptet0, ptet;

    T = manifold->num_tetrahedra;
    eqns = NEW_ARRAY(3 * T, int);

    for (i = 0; i < 3 * T; i++)
        eqns[i] = 0;

    /*
     *  Build edge equations.
     */

    for (edge = manifold->edge_list_begin.next; edge != &manifold->edge_list_end; edge = edge->next) {
        if (edge->index == edgeClass)
            break;
    }

    set_left_edge(edge, &ptet0);
    ptet = ptet0;
    do {
        eqns[3 * ptet.tet->index + edge3_between_faces[ptet.near_face][ptet.left_face]]++;
        veer_left(&ptet);
    } while (same_positioned_tet(&ptet, &ptet0) == FALSE);

    return eqns;
}

/*
 * Construct the symplectic equations from the oscillating curves
 */

int *combinatorial_holonomy(Triangulation *manifold, int edge_class) {
    int v, f, ff;
    int *eqns = NEW_ARRAY(3 * manifold->num_tetrahedra, int);
    Tetrahedron *tet;

    for (int i = 0; i < 3 * manifold->num_tetrahedra; i++) {
        eqns[i] = 0;
    }

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which tet vertex
        for (v = 0; v < 4; v++) {
            // which face
            for (f = 0; f < 4; f++) {
                if (f == v)
                    continue;

                ff = (int) remaining_face[v][f];

                eqns[3 * tet->index + edge3_between_faces[f][ff]]
                    += FLOW(tet->extra[edge_class].curve[v][f], tet->extra[edge_class].curve[v][ff]);
            }
        }
    }

    return eqns;
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

/*
 * Give each edge of the triangulation an index to identify the cusp vertices
 */

void label_triangulation_edges(Triangulation *manifold) {
    int i = 0;
    EdgeClass *edge = &manifold->edge_list_begin;

    while ((edge = edge->next)->next != NULL)
        edge->index = i++;

    // incorrect number of edge classes
    if (i != manifold->num_tetrahedra)
        uFatalError("label_triangulation_edges", "symplectic_basis");
}
