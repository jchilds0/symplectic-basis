/*
 * cusp_structure.c
 *
 * Provides the functions
 *
 *      CuspStructure           *init_cusp_structure(Triangulation *, Cusp *, EndMultiGraph *);
 *
 *      void                    free_cusp_structure(CuspStructure *);
 *
 * which are used by oscillating_curves.c to initialise and free the cusp structure.
 * init_cusp_structure call the init functions for cusp triangulation, cusp regions,
 * and train lines. It also provides
 *
 *      void                    construct_cusp_region_dual_graph(CuspStructure *);
 *
 * for updating the dual graph after constructing an oscillating curve.
 */

#include "symplectic_kernel.h"

void                    init_cusp_triangulation(Triangulation *, CuspStructure *);
void                    find_intersection_triangle(Triangulation *, CuspStructure *);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *, CuspTriangle *, int);

/*
 * peripheral_curves.c places a meridian and longitude curve on each cusp. It
 * starts at a base triangle, the intersection point, and searches outwards.
 * Note it does not visit a cusp triangle more than once. So we find a cusp
 * triangle which contains both a meridian and longitude (this should be the
 * same intersection triangle that peripheral_curves sets since it is the same
 * search process) and assert this is the intersection triangle. Currently
 * init_cusp_regions assumes the intersection triangle only contains curves
 * which intersect. This is because we need some information about the curves
 * to construct the cusp regions.
 */

void find_intersection_triangle(Triangulation *manifold, CuspStructure *boundary) {
    FaceIndex   face;
    Cusp *cusp = boundary->cusp;
    int n;

    for (cusp->basepoint_tet = manifold->tet_list_begin.next;
         cusp->basepoint_tet != &manifold->tet_list_end;
         cusp->basepoint_tet = cusp->basepoint_tet->next)

        for (cusp->basepoint_vertex = 0;
             cusp->basepoint_vertex < 4;
             cusp->basepoint_vertex++)
        {
            if (cusp->basepoint_tet->cusp[cusp->basepoint_vertex] != cusp)
                continue;

            for (face = 0; face < 4; face++)
            {
                if (face == cusp->basepoint_vertex)
                    continue;

                for (n = 0; n < 2; n++) {
                    cusp->basepoint_orientation = ORIENTATION(n);

                    if (cusp->basepoint_tet->curve
                        [M]
                        [cusp->basepoint_orientation]
                        [cusp->basepoint_vertex]
                        [face] != 0
                        && cusp->basepoint_tet->curve
                           [L]
                           [cusp->basepoint_orientation]
                           [cusp->basepoint_vertex]
                           [face] != 0) {
                        /*
                         *  We found the basepoint!
                         */

                        boundary->intersect_tet_index  = cusp->basepoint_tet->index;
                        boundary->intersect_tet_vertex = cusp->basepoint_vertex;
                        return;
                    }


                }
            }
        }
}

/*
 * Calculate the number of curves passing around a vertex in the cusp
 * triangulation.
 */

CuspStructure *init_cusp_structure(Triangulation *manifold, Cusp *cusp, EndMultiGraph *multi_graph) {
    CuspStructure *boundary = NEW_STRUCT(CuspStructure);

    // Invalid cusp topology
    if (cusp->topology == Klein_cusp)
        uFatalError("init_cusp_structure", "symplectic_basis");

    boundary->manifold              = manifold;
    boundary->cusp                  = cusp;
    boundary->num_edge_classes      = manifold->num_tetrahedra;
    boundary->num_cusp_triangles    = 0;
    boundary->num_cusp_regions      = 0;

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);
    init_train_line(boundary);

    boundary->dual_graph = NULL;
    construct_cusp_region_dual_graph(boundary);

    return boundary;
}

void free_cusp_structure(CuspStructure *cusp) {
    CuspTriangle *tri;

    while (cusp->cusp_triangle_begin.next != &cusp->cusp_triangle_end) {
        tri = cusp->cusp_triangle_begin.next;
        REMOVE_NODE(tri)
        my_free(tri);
    }

    free_graph(cusp->dual_graph);
    free_cusp_region(cusp);
    free_train_line(cusp);
    my_free(cusp);
}

/*
 * Returns a pointer to the cusp triangle which is the neighbour of tri across
 * face 'face'.
 */

CuspTriangle *find_cusp_triangle(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end,
                                 CuspTriangle *tri, int face) {
    int tet_index, tet_vertex;
    CuspTriangle *pTri;

    tet_index = tri->tet->neighbor[face]->index;
    tet_vertex = EVALUATE(tri->tet->gluing[face], tri->tet_vertex);

    for (pTri = cusp_triangle_begin->next; pTri != cusp_triangle_end; pTri = pTri->next) {
        if (pTri->tet_index == tet_index && pTri->tet_vertex == tet_vertex)
            return pTri;
    }

    // Didn't find a neighbour
    return NULL;
}

/*
 * Construct the cusp triangle doubly linked list which consists of the
 * triangles in the cusp triangulation
 */

void init_cusp_triangulation(Triangulation *manifold, CuspStructure *cusp) {
    int index = 0;
    VertexIndex vertex;
    FaceIndex face;
    Tetrahedron *tet;
    CuspTriangle *tri;

    // Allocate Cusp Triangulation Header and Tail Null nodes
    cusp->cusp_triangle_begin.next      = &cusp->cusp_triangle_end;
    cusp->cusp_triangle_begin.prev      = NULL;
    cusp->cusp_triangle_end.next        = NULL;
    cusp->cusp_triangle_end.prev        = &cusp->cusp_triangle_begin;

    // which tetrahedron are we on
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // while vertex are we on
        for (vertex = 0; vertex < 4; vertex++) {
            // is this vertex on the right cusp
            if (tet->cusp[vertex] != cusp->cusp) {
                continue;
            }

            tri = NEW_STRUCT( CuspTriangle );
            INSERT_BEFORE(tri, &cusp->cusp_triangle_end)
            index++;

            tri->tet = tet;
            tri->cusp = tet->cusp[vertex];
            tri->tet_index = tri->tet->index;
            tri->tet_vertex = vertex;

            for (face = 0; face < 4; face ++) {
                if (tri->tet_vertex == face)
                    continue;

                tri->vertices[face].v1              = tri->tet_vertex;
                tri->vertices[face].v2              = face;
                tri->vertices[face].edge            = tri->tet->edge_class[
                        edge_between_vertices[tri->vertices[face].v1][tri->vertices[face].v2]];
                tri->vertices[face].edge_class      = tri->vertices[face].edge->index;
                tri->vertices[face].edge_index      = -1;
            }
        }
    }

    // which cusp triangle
    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        // which vertex
        for (face = 0; face < 4; face++) {
            if (face == tri->tet_vertex)
                continue;

            tri->neighbours[face] = find_cusp_triangle(&cusp->cusp_triangle_begin, &cusp->cusp_triangle_end, tri, face);
        }
    }

    label_cusp_vertex_indices(&cusp->cusp_triangle_begin, &cusp->cusp_triangle_end, cusp->num_edge_classes);
    cusp->num_cusp_triangles = index;
}

/*
 * Each edge class of the manifold appears as two vertices in the cusp
 * triangulation. We iterate over the cusp triangulation, walking around each
 * vertex to give it the same index.
 */

void label_cusp_vertex_indices(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end, int numEdgeClasses) {
    int i, vertex;
    CuspTriangle *tri;

    int *current_index = NEW_ARRAY(numEdgeClasses, int);

    for (i = 0; i < numEdgeClasses; i++)
        current_index[i] = 0;

    for (tri = cusp_triangle_begin->next; tri != cusp_triangle_end; tri = tri->next) {
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tet_vertex || tri->vertices[vertex].edge_index != -1)
                continue;

            walk_around_cusp_vertex(tri, vertex, current_index[tri->vertices[vertex].edge_class]);
            current_index[tri->vertices[vertex].edge_class]++;
        }
    }

    my_free(current_index);
}

/*
 * Walk around vertex cusp_vertex of triangle *tri and set edge_index to index.
 */

void walk_around_cusp_vertex(CuspTriangle *tri, int cusp_vertex, int index) {
    int gluing_vertex, outside_vertex, old_gluing_vertex, old_cusp_vertex, old_outside_vertex;
    gluing_vertex = (int) remaining_face[cusp_vertex][tri->tet_vertex];
    outside_vertex = (int) remaining_face[tri->tet_vertex][cusp_vertex];

    while (tri->vertices[cusp_vertex].edge_index == -1) {
        tri->vertices[cusp_vertex].edge_index = index;

        // Move to the next cusp triangle
        old_cusp_vertex         = cusp_vertex;
        old_gluing_vertex       = gluing_vertex;
        old_outside_vertex      = outside_vertex;

        cusp_vertex             = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex           = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri                     = tri->neighbours[old_gluing_vertex];
    }
}

/*
 * Construct the graph dual to the cusp regions, using region->index to label
 * each vertex, and adding edges using region->adj_cusp_regions[].
 */

void construct_cusp_region_dual_graph(CuspStructure *cusp) {
    int i, face;
    CuspRegion *region;

    Graph *graph1 = init_graph(cusp->num_cusp_regions, FALSE);
    cusp->dual_graph_regions = NEW_ARRAY(cusp->num_cusp_regions, CuspRegion *);

    int *visited = NEW_ARRAY(graph1->num_vertices, int);

    for (i = 0; i < graph1->num_vertices; i++) {
        visited[i] = FALSE;
        cusp->dual_graph_regions[i] = NULL;
    }

    // Walk around the cusp triangulation inserting edges
    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            if (visited[region->index])
                continue;

            for (face = 0; face < 4; face++) {
                if (!region->adj_cusp_triangle[face])
                    continue;

                // Missing adj region data
                if (region->adj_cusp_regions[face] == NULL)
                    uFatalError("construct_cusp_region_dual_graph", "symplectic_basis");

                insert_edge(graph1, region->index, region->adj_cusp_regions[face]->index, graph1->directed);
                cusp->dual_graph_regions[region->index] = region;
            }

            visited[region->index] = 1;
        }
    }

    free_graph(cusp->dual_graph);
    my_free(visited);

    cusp->dual_graph = graph1;
}
