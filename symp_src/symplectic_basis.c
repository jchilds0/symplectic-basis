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

#include "typedefs.h"
#include <stdio.h>
#include <string.h>
#include "queue.h"
#include "end_multi_graph.h"

#define ATLEAST_TWO(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))
#define TRI_TO_INDEX(tet_index, tet_vertex)     (4 * (tet_index) + (tet_vertex))

#define COPY_PATH_ENDPOINT(new, old)    {                                                       \
                                            (new)->vertex = (old)->vertex;                      \
                                            (new)->face = (old)->face;                          \
                                            (new)->tri = (old)->tri;                            \
                                            (new)->region_index = (old)->region_index;          \
                                            (new)->region = (old)->region;                      \
                                            (new)->node = (old)->node;                          \
                                            (new)->num_adj_curves = (old)->num_adj_curves;      \
                                        }

#define COPY_PATH_NODE(new, old)        {                                                           \
                                            (new)->next = NULL;                                     \
                                            (new)->prev = NULL;                                     \
                                            (new)->next_face = (old)->next_face;                    \
                                            (new)->prev_face = (old)->prev_face;                    \
                                            (new)->inside_vertex = (old)->inside_vertex;            \
                                            (new)->cusp_region_index = (old)->cusp_region_index;    \
                                            (new)->tri = (old)->tri;                                \
                                        }


static int debug = 0;

static int edgesThreeToFour[4][3] = {{1, 2, 3},
                                     {0, 2, 3},
                                     {0, 1, 3},
                                     {0, 1, 2}};


/**
 * Symplectic Basis
 */

int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     *combinatorial_holonomy(Triangulation *, int);
void                    oscillating_curves(Triangulation *manifold, Boolean *edge_classes);

/**
 * Initialisation Functions
 */

CuspStructure           *init_cusp_structure(Triangulation *, Cusp *);
void                    free_cusp_structure(CuspStructure **, int, int);
void                    init_cusp_triangulation(Triangulation *, CuspStructure *);
void                    init_cusp_region(CuspStructure *);
int                     init_intersect_cusp_region(CuspStructure *, CuspTriangle *, int);
int                     init_intersect_vertex_two_zero_flows(CuspStructure *, CuspTriangle *, int);
int                     init_normal_cusp_region(CuspStructure *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspStructure *, CuspTriangle *, const int [4], const Boolean [4], int);
void                    init_train_line(CuspStructure *);
CurveComponent          *init_curve_component(int, int, int);
OscillatingCurves       *init_oscillating_curves(Triangulation *, const Boolean *);
void                    free_oscillating_curves(OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, CuspStructure *);

/**
 * Cusp Functions
 */

int                     net_flow_around_vertex(CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *, CuspTriangle *, int);
void                    update_adj_region_data(CuspStructure *);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    construct_cusp_region_dual_graph(CuspStructure *);
void                    log_structs(Triangulation *, CuspStructure **, OscillatingCurves *, char *);

/**
 * Train lines
 */

void                    do_manifold_train_lines(Triangulation *, CuspStructure **, EndMultiGraph *);
int                     *find_tet_index_for_edge_classes(Triangulation *, const Boolean *);
void                    find_edge_class_edges(Triangulation *, CuspStructure **, Boolean *);
void                    find_edge_class_edges_on_cusp(CuspStructure *, const Boolean *, const int *);
Boolean                 *update_edge_classes_on_cusp(CuspStructure **, Boolean **, int, int, int);

void                    find_primary_train_line(CuspStructure *, Boolean *);
void                    do_initial_train_line_segment_on_cusp(CuspStructure *, PathEndPoint *, PathEndPoint *);
void                    do_train_line_segment_on_cusp(CuspStructure *, PathEndPoint *, PathEndPoint *);
void                    extended_train_line_path(CuspStructure *, PathEndPoint *, PathEndPoint *, EdgeNode *, EdgeNode *);
void                    path_finding_with_loops(CuspStructure *, PathEndPoint *, PathEndPoint *, int, int, EdgeNode *, EdgeNode *);
void                    cycle_path(Graph *, EdgeNode *, EdgeNode *, int, int, int, int, int);
void                    graph_path_to_path_node(CuspStructure *, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_regions_along_train_line_segment(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_region_train_line_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_train_line_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int);

Boolean                 *edge_classes_on_cusp(CuspStructure *, const Boolean *);
PathEndPoint            *next_valid_endpoint_index(CuspStructure *, PathEndPoint *);
void                    tri_endpoint_to_region_endpoint(CuspStructure *, PathEndPoint *);
Boolean                 array_contains_true(const Boolean *, int);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(CuspStructure **, OscillatingCurves *, EndMultiGraph *);
void                    do_one_oscillating_curve(CuspStructure **, OscillatingCurves *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int, int);
CurveComponent          *setup_train_line_component(CuspStructure *, EndMultiGraph *, CurveComponent *, CurveComponent *, CuspEndPoint *, int);
void                    do_curve_component_on_train_line(CuspStructure *, CurveComponent *);
CurveComponent          *setup_first_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
CurveComponent          *setup_last_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
void                    do_curve_component_to_new_edge_class(CuspStructure *, CurveComponent *);
void                    find_single_endpoint(CuspStructure *, PathEndPoint *, int, int);
void                    find_single_matching_endpoint(CuspStructure *, PathEndPoint *, PathEndPoint *, int, int);
void                    find_train_line_endpoint(CuspStructure *, PathEndPoint *, int, int, int, Boolean);

void                    graph_path_to_dual_curve(CuspStructure *, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    endpoint_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *, PathEndPoint *, int);
void                    interior_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *);

void                    split_cusp_regions_along_path(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_path_len_one(CuspStructure *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_region_path_interior(CuspRegion *, CuspRegion *, PathNode *, int);
void                    split_cusp_region_path_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, PathNode *, int);

void                    update_adj_curve_along_path(CuspStructure **, OscillatingCurves *, int, Boolean);
void                    update_adj_curve_at_endpoint(PathEndPoint *, CurveComponent *, int);
void                    update_adj_curve_on_cusp(CuspStructure *);
void                    update_path_holonomy(CurveComponent *, int);

/**
 * Dual Curves
 *
 * Each oscillating curve contributes combinatorial holonomy, we store this in
 * curve[4][4] in a similar way to the curve[4][4] attribute of a Tetrahedron.
 * An array of size num_edge_classes is attached to each Tetrahedron.
 * tet->extra[edge_class]->curve[v][f] is the intersection number of
 * the oscillating curve associated to edge_class with the face 'f' of the
 * cusp triangle at vertex 'v' of tet.
 */

struct extra {
    int                         curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

// Symplectic Basis

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

    peripheral_curves(manifold);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        if (tet->extra != NULL)
            uFatalError("oscillating_curves", "symplectic_basis");

        tet->extra = NEW_ARRAY(manifold->num_tetrahedra, Extra);

        for (i = 0; i < manifold->num_tetrahedra; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    tet->extra[i].curve[j][k] = 0;
    }

    // Dual Edge Curves Gamma_i -> symplectic equations
    oscillating_curves(manifold, edge_classes);

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

/*
 * Initialise cusp structure on each cusp, construct train lines, construct
 * oscillating curves and store the intersection numbers of each curve with the
 * cusp triangles it enters in tet->extra[edge_class]->curve, in the same fashion
 * as the peripheral curves.
 */

void oscillating_curves(Triangulation *manifold, Boolean *edge_classes) {
    int i;
    label_triangulation_edges(manifold);

    CuspStructure **cusps         = NEW_ARRAY(manifold->num_cusps, CuspStructure *);
    EndMultiGraph *multi_graph    = init_end_multi_graph(manifold);
    Cusp *cusp;

    for (i = 0; i < multi_graph->num_edge_classes; i++)
        edge_classes[i] = multi_graph->edge_classes[i] == TRUE ? FALSE : TRUE;

    edge_classes[multi_graph->e0] = FALSE;

    OscillatingCurves *curves   = init_oscillating_curves(manifold, edge_classes);

    for (i = 0; i < manifold->num_cusps; i++) {
        for (cusp = manifold->cusp_list_begin.next; cusp != &manifold->cusp_list_end && cusp->index != i; cusp = cusp->next);

        if (cusp == &manifold->cusp_list_end)
            uFatalError("oscillating_curves", "symplectic_basis");

        cusps[i] = init_cusp_structure(manifold, cusp);
    }

    if (debug) {
        printf("\n");
        printf("Struct Initialisation\n");
        printf("\n");

        log_structs(manifold, cusps, NULL, "gluing");
        log_structs(manifold, cusps, NULL, "homology");
        log_structs(manifold, cusps, NULL, "edge_indices");
        log_structs(manifold, cusps, NULL, "inside_edge");
        log_structs(manifold, cusps, NULL, "cusp_regions");
    }

    do_manifold_train_lines(manifold, cusps, multi_graph);
    do_oscillating_curves(cusps, curves, multi_graph);

    if (debug) {
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("%d, ", cusps[i]->num_cusp_regions);
        }
        printf("\n");
    }

    free_end_multi_graph(multi_graph);
    free_oscillating_curves(curves);
    free_cusp_structure(cusps, manifold->num_cusps, manifold->num_tetrahedra);
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

// ------------------------------------

/*
 * Initialisation Functions
 */

CuspStructure *init_cusp_structure(Triangulation *manifold, Cusp *cusp) {
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

void free_cusp_structure(CuspStructure **cusps, int num_cusps, int num_edge_classes) {
    int cusp_index;
    CuspTriangle *tri;
    CuspRegion *region;
    PathNode *path_node;
    CuspStructure *cusp;

    for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
        cusp = cusps[cusp_index];
        // free graph
        free_graph(cusp->dual_graph);

        // free cusp regions
        for (int i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
            while (cusp->cusp_region_begin[i].next != &cusp->cusp_region_end[i]) {
                region = cusp->cusp_region_begin[i].next;
                REMOVE_NODE(region)
                my_free(region);
            }
        }

        my_free(cusp->cusp_region_begin);
        my_free(cusp->cusp_region_end);

        // free cusp triangle
        while (cusp->cusp_triangle_begin.next != &cusp->cusp_triangle_end) {
            tri = cusp->cusp_triangle_begin.next;
            REMOVE_NODE(tri)
            my_free(tri);
        }

        // free train line path
        while (cusp->train_line_path_begin.next != &cusp->train_line_path_end) {
            path_node = cusp->train_line_path_begin.next;
            REMOVE_NODE(path_node)
            my_free(path_node);
        }

        my_free(cusp->train_line_endpoint[0]);
        my_free(cusp->train_line_endpoint[1]);
        my_free(cusp);
    }

    my_free(cusps);
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

            tri->num_curves = net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][0])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][1])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][2]);

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
 * Initialise the cusp region doubly linked list to cotain the regions bounded 
 * by the meridian and longitude curves.
 */

void init_cusp_region(CuspStructure *cusp) {
    int index;
    CuspTriangle *tri;

    // Header and tailer nodes.
    cusp->cusp_region_begin = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);
    cusp->cusp_region_end   = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);

    for (index = 0; index < 4 * cusp->manifold->num_tetrahedra; index++) {
        cusp->cusp_region_begin[index].next    = &cusp->cusp_region_end[index];
        cusp->cusp_region_begin[index].prev    = NULL;
        cusp->cusp_region_end[index].next      = NULL;
        cusp->cusp_region_end[index].prev      = &cusp->cusp_region_begin[index];
    }

    index = 0;
    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        // Intersection vertex doesn't have a center
        if (tri->tet_index == cusp->intersect_tet_index && tri->tet_vertex == cusp->intersect_tet_vertex) {
            index = init_intersect_cusp_region(cusp, tri, index);
            continue;
        }

        index = init_normal_cusp_region(cusp, tri, index);
    }

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

/*
 * Assume peripheral_curves() has been called, and as a result the only curves 
 * on the intersection triangle are those which intersect, and they give a 
 * valid intersection.
 */

int init_intersect_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, curve_index, vertex, v1, v2, v3;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 1; i < net_flow_around_vertex(tri, vertex); i++) {
            for (curve_index = 0; curve_index < 2; curve_index++) {
                distance[v1]                    = i;
                distance[v2]                    = MIN(distance[v1], 2 * net_flow_around_vertex(tri, vertex) - distance[v1])
                                                + net_flow_around_vertex(tri, v2) + net_flow_around_vertex(tri, v1);
                distance[vertex]                = net_flow_around_vertex(tri, vertex)
                                                - distance[v1] + net_flow_around_vertex(tri, v1);
                distance[tri->tet_vertex]       = -1;

                adj_triangle[v1]                = 1;
                adj_triangle[v2]                = 0;
                adj_triangle[vertex]            = 0;
                adj_triangle[tri->tet_vertex]   = -1;

                set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
                index++;

                // Swap vertices
                v1 = (int) remaining_face[vertex][tri->tet_vertex];
                v2 = (int) remaining_face[tri->tet_vertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (net_flow_around_vertex(tri, v1) && net_flow_around_vertex(tri, v2)) {
            distance[v1]                    = net_flow_around_vertex(tri, v2);
            distance[v2]                    = net_flow_around_vertex(tri, v1);
            distance[vertex]                = MIN(net_flow_around_vertex(tri, v1) + distance[v1],
                                              net_flow_around_vertex(tri, v2) + distance[v2])
                                                      + net_flow_around_vertex(tri, vertex);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[v1]                = 0;
            adj_triangle[v2]                = 0;
            adj_triangle[vertex]            = 1;
            adj_triangle[tri->tet_vertex]   = -1;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }
    }

    // Region of distance 0 to vertex
    v1 = edgesThreeToFour[tri->tet_vertex][0];
    v2 = edgesThreeToFour[tri->tet_vertex][1];
    v3 = edgesThreeToFour[tri->tet_vertex][2];

    // Edge Case: Two vertices with 0 flow
    if (ATLEAST_TWO(!net_flow_around_vertex(tri, v1),
                    !net_flow_around_vertex(tri, v2),
                    !net_flow_around_vertex(tri, v3)))
        return init_intersect_vertex_two_zero_flows(cusp, tri, index);

    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        distance[vertex]               = 0;
        distance[v1]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v1);
        distance[v2]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v2);
        distance[tri->tet_vertex]      = -1;

        adj_triangle[vertex]           = 0;
        adj_triangle[v1]               = 1;
        adj_triangle[v2]               = 1;
        adj_triangle[tri->tet_vertex]  = 0;

        set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
        index++;
    }

    return index;
}

int init_intersect_vertex_two_zero_flows(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int vertex, v1, v2, v3, distance[4];
    Boolean adj_triangle[4];

    v1 = (int) edgesThreeToFour[tri->tet_vertex][0];
    v2 = (int) edgesThreeToFour[tri->tet_vertex][1];
    v3 = (int) edgesThreeToFour[tri->tet_vertex][2];

    distance[v1]                   = net_flow_around_vertex(tri, v1);
    distance[v2]                   = net_flow_around_vertex(tri, v2);
    distance[v3]                   = net_flow_around_vertex(tri, v3);
    distance[tri->tet_vertex]      = -1;

    adj_triangle[v1]               = 1;
    adj_triangle[v2]               = 1;
    adj_triangle[v3]               = 1;
    adj_triangle[tri->tet_vertex]  = -1;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;

    // Find vertex with non-zero flow
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        if (net_flow_around_vertex(tri, vertex)) {
            v1 = vertex;
            v2 = (int) remaining_face[tri->tet_vertex][v1];
            v3 = (int) remaining_face[v1][tri->tet_vertex];
            break;
        }
    }
    distance[v1]                    = 0;
    distance[v2]                    = net_flow_around_vertex(tri, v1);
    distance[v3]                    = net_flow_around_vertex(tri, v1);
    distance[tri->tet_vertex]       = -1;

    adj_triangle[v1]                = 0;
    adj_triangle[v2]                = 1;
    adj_triangle[v3]                = 1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);

    return index + 1;
}

int init_normal_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, vertex, v1, v2;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 0; i < net_flow_around_vertex(tri, vertex); i++) {
            distance[vertex]                = i;
            distance[v1]                    = net_flow_around_vertex(tri, v1)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[v2]                    = net_flow_around_vertex(tri, v2)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[vertex]            = 0;
            adj_triangle[v1]                = 1;
            adj_triangle[v2]                = 1;
            adj_triangle[tri->tet_vertex]   = 0;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }

    }

    // center region
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        distance[vertex]        = net_flow_around_vertex(tri, vertex);
        adj_triangle[vertex]    = 1;
    }

    distance[tri->tet_vertex]       = -1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which allocates the attributes of the 
 * cusp region
 */

void set_cusp_region_data(CuspStructure *cusp, CuspTriangle *tri, const int distance[4],
                          const Boolean adj_cusp_triangle[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *region = NEW_STRUCT( CuspRegion );
    INSERT_BEFORE(region, &cusp->cusp_region_end[TRI_TO_INDEX(tri->tet_index, tri->tet_vertex)])

    region->tri             = tri;
    region->tet_index       = region->tri->tet_index;
    region->tet_vertex      = region->tri->tet_vertex;
    region->index           = index;

    // default values
    for (i = 0; i < 4; i++) {
        region->adj_cusp_triangle[i] = FALSE;
        region->adj_cusp_regions[i]  = NULL;

        for (j = 0; j < 4; j++) {
            region->curve[i][j]             = -1;
            region->dive[i][j]              = 0;
            region->num_adj_curves[i][j]    = 0;
            region->temp_adj_curves[i][j]   = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tet_vertex][i];
        v2 = edgesThreeToFour[tri->tet_vertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tet_vertex][(i + 2) % 3];

        region->curve[v2][v1]   = distance[v1];
        region->curve[v3][v1]   = distance[v1];
        region->dive[v2][v1]    = distance[v1] ? FALSE : TRUE;
        region->dive[v3][v1]    = distance[v1] ? FALSE : TRUE;

        region->adj_cusp_triangle[v1] = adj_cusp_triangle[v1];
    }
}

void init_train_line(CuspStructure *cusp) {
    int edge_class, edge_index;

    cusp->train_line_path_begin.next    = &cusp->train_line_path_end;
    cusp->train_line_path_begin.prev    = NULL;
    cusp->train_line_path_end.next      = NULL;
    cusp->train_line_path_end.prev      = &cusp->train_line_path_begin;

    cusp->train_line_endpoint[0] = NEW_ARRAY(cusp->manifold->num_tetrahedra, PathEndPoint);
    cusp->train_line_endpoint[1] = NEW_ARRAY(cusp->manifold->num_tetrahedra, PathEndPoint);

    for (edge_class = 0; edge_class < cusp->manifold->num_tetrahedra; edge_class++) {
        for (edge_index = 0; edge_index < 2; edge_index++) {
            cusp->train_line_endpoint[edge_index][edge_class].tri               = NULL;
            cusp->train_line_endpoint[edge_index][edge_class].region            = NULL;
            cusp->train_line_endpoint[edge_index][edge_class].node              = NULL;
            cusp->train_line_endpoint[edge_index][edge_class].num_adj_curves    = 0;
        }
    }
}

CurveComponent *init_curve_component(int edge_class_start, int edge_class_finish, int cusp_index) {
    int i;

    CurveComponent *path = NEW_STRUCT(CurveComponent );

    path->path_begin.next     = &path->path_end;
    path->path_begin.prev     = NULL;
    path->path_end.next       = NULL;
    path->path_end.prev       = &path->path_begin;

    path->edge_class[START]     = edge_class_start;
    path->edge_class[FINISH]    = edge_class_finish;
    path->cusp_index            = cusp_index;

    for (i = 0; i < 2; i++) {
        path->endpoints[i].tri              = NULL;
        path->endpoints[i].region           = NULL;
        path->endpoints[i].num_adj_curves   = 0;
    }

    return path;
}

/*
 * Initialise dual curve doubly linked list which stores the oscillating curves
 * on the cusp
 */

OscillatingCurves *init_oscillating_curves(Triangulation *manifold, const Boolean *edge_classes) {
    int i, j;
    OscillatingCurves *curves = NEW_STRUCT(OscillatingCurves );

    curves->num_curves = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++)
        if (edge_classes[i])
            curves->num_curves++;

    curves->curve_begin               = NEW_ARRAY(curves->num_curves, CurveComponent );
    curves->curve_end                 = NEW_ARRAY(curves->num_curves, CurveComponent );
    curves->edge_class                = NEW_ARRAY(curves->num_curves, int);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i])
            continue;

        curves->edge_class[j] = i;
        j++;
    }

    // which curve
    for (i = 0; i < curves->num_curves; i++) {
        curves->curve_begin[i].next    = &curves->curve_end[i];
        curves->curve_begin[i].prev    = NULL;
        curves->curve_end[i].next      = NULL;
        curves->curve_end[i].prev      = &curves->curve_begin[i];
    }

    return curves;
}

void free_oscillating_curves(OscillatingCurves *curves) {
    int i;
    CurveComponent *path;
    PathNode *path_node;

    for (i = 0; i < curves->num_curves; i++) {
        while (curves->curve_begin[i].next != &curves->curve_end[i]) {
            path = curves->curve_begin[i].next;
            REMOVE_NODE(path)

            while (path->path_begin.next != &path->path_end) {
                path_node = path->path_begin.next;
                REMOVE_NODE(path_node)
                my_free(path_node);
            }

            my_free(path);
        }
    }

    my_free(curves->curve_begin);
    my_free(curves->curve_end);
    my_free(curves->edge_class);
    my_free(curves);
}

// ----------------------------------------------------

/*
 * Cusp Functions
 */

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

int net_flow_around_vertex(CuspTriangle *tri, int vertex) {
    int mflow, lflow, retval;

    // Contribution from meridian curves
    mflow = FLOW(tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    retval = ABS(mflow) + ABS(lflow);
    return retval;
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
 * Calculate which regions are located across cusp edges and store the result
 * in the adj_cusp_regions attribute
 */

void update_adj_region_data(CuspStructure *cusp) {
    CuspTriangle *adj_triangle;
    CuspRegion *region;
    FaceIndex f;
    int i, adj_index;

    // Add adjacent region info
    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            for (f = 0; f < 4; f++) {
                if (!region->adj_cusp_triangle[f] || region->tet_vertex == f) {
                    region->adj_cusp_regions[f] = NULL;
                    continue;
                }

                adj_triangle = region->tri->neighbours[f];
                adj_index = TRI_TO_INDEX(adj_triangle->tet_index, adj_triangle->tet_vertex);
                region->adj_cusp_regions[f] = find_adj_region(&cusp->cusp_region_begin[adj_index],
                                                              &cusp->cusp_region_end[adj_index],
                                                              region, f);
            }
        }
    }
}

/*
 * Find the cusp region which is adjacent to x across face.
 */

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
                            CuspRegion *x, int face) {
    int v1, v2, y_vertex1, y_vertex2, y_face, distance_v1, distance_v2, tet_index, tet_vertex;
    Boolean adj_face;
    CuspTriangle *tri = x->tri;
    CuspRegion *region;

    v1 = (int) remaining_face[tri->tet_vertex][face];
    v2 = (int) remaining_face[face][tri->tet_vertex];

    y_vertex1    = EVALUATE(tri->tet->gluing[face], v1);
    y_vertex2    = EVALUATE(tri->tet->gluing[face], v2);
    y_face       = EVALUATE(tri->tet->gluing[face], face);

    // Check current adj region first
    if (x->adj_cusp_regions[face] != NULL) {
        distance_v1      = (x->curve[face][v1] == x->adj_cusp_regions[face]->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == x->adj_cusp_regions[face]->curve[y_face][y_vertex2]);
        adj_face         = x->adj_cusp_regions[face]->adj_cusp_triangle[y_face];

        if (distance_v1 && distance_v2 && adj_face)
            return x->adj_cusp_regions[face];
    }

    /*
     * We search through the regions in reverse as the new regions
     * are added to the end of the doubly linked list
     */
    for (region = cusp_region_end->prev; region != cusp_region_begin; region = region->prev) {
        tet_index    = (tri->neighbours[face]->tet_index == region->tet_index);
        tet_vertex   = (tri->neighbours[face]->tet_vertex == region->tet_vertex);

        if (!tet_index || !tet_vertex)
            continue;

        distance_v1      = (x->curve[face][v1] == region->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == region->curve[y_face][y_vertex2]);
        adj_face         = region->adj_cusp_triangle[y_face];

        // missing distance
        if (region->curve[y_face][y_vertex1] == -1 || region->curve[y_face][y_vertex2] == -1)
            uFatalError("find_adj_region", "symplectic_basis");

        if (distance_v1 && distance_v2 && adj_face)
            return region;
    }

    // We didn't find a cusp region
    //uFatalError("find_cusp_region", "symplectic_basis");
    return NULL;
}

/*
 * region1 splits into region1 and region2, set them up to be split
 */

void copy_region(CuspRegion *region1, CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri            = region1->tri;
    region2->tet_index      = region1->tet_index;
    region2->tet_vertex     = region1->tet_vertex;

    for (i = 0; i < 4; i++) {
        region2->adj_cusp_triangle[i]   = region1->adj_cusp_triangle[i];
        region2->adj_cusp_regions[i]    = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]            = region1->curve[i][j];
            region2->dive[i][j]             = FALSE;
            region2->num_adj_curves[i][j]   = region1->num_adj_curves[i][j];
            region2->temp_adj_curves[i][j]  = region1->temp_adj_curves[i][j];
        }
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

/*
 * Types: gluing, train_lines, cusp_regions, homology, edge_indices,
 * dual_curves, inside_edge, graph, endpoints
 */

void log_structs(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves, char *type) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    CuspTriangle *tri;
    CuspRegion *region;
    EdgeNode *edge_node;
    PathNode *path_node;
    CurveComponent *path;
    Graph *g;
    CuspStructure *cusp;
    PathEndPoint *endpoint;

    if (strcmp(type, "gluing") == 0) {
        printf("Triangle gluing info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);
            cusp = cusps[i];

            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                for (j = 0; j < 4; j++) {
                    if (j == tri->tet_vertex)
                        continue;

                    x_vertex1 = (int) remaining_face[tri->tet_vertex][j];
                    x_vertex2 = (int) remaining_face[j][tri->tet_vertex];
                    y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                    y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                    printf("    (Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                           "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                           tri->tet_index,               // Tet Index
                           tri->tet_vertex,                // Tet Vertex
                           j,      // Cusp Edge
                           tri->tet->neighbor[j]->index,                              // Tet Index
                           EVALUATE(tri->tet->gluing[j], tri->tet_vertex),             // Tet Vertex
                           EVALUATE(tri->tet->gluing[j], j),   // Cusp Edge
                           x_vertex1, y_vertex1,
                           x_vertex2, y_vertex2
                    );
                }
            }
        }
    } else if (strcmp(type, "train_lines") == 0) {
        printf("Train Lines\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            printf("    Train Line Path: \n");

            for (path_node = cusp->train_line_path_begin.next; path_node != &cusp->train_line_path_end; path_node = path_node->next) {
                printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                       path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                       path_node->next_face, path_node->prev_face, path_node->inside_vertex
                );
            }

            printf("    Train Line Endpoints\n");
            for (j = 0; j < cusp->num_edge_classes; j++) {
                for (k = 0; k < 2; k++) {
                    if (cusp->train_line_endpoint[k][j].tri == NULL)
                        continue;

                    endpoint = &cusp->train_line_endpoint[k][j];
                    printf("        Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d)\n",
                           endpoint->region_index, endpoint->tri->tet_index,
                           endpoint->tri->tet_vertex, endpoint->face, endpoint->vertex,
                           endpoint->tri->vertices[endpoint->vertex].edge_class,
                           endpoint->tri->vertices[endpoint->vertex].edge_index);
                }
            }
        }
    } else if (strcmp(type, "cusp_regions") == 0) {
        printf("Cusp Region info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (j = 0; j < 4 * cusp->manifold->num_tetrahedra; j++) {
                printf("    Cusp Triangle (Tet Index %d Tet Vertex %d)\n", j / 4, j % 4);
                for (region = cusp->cusp_region_begin[j].next;
                     region != &cusp->cusp_region_end[j]; region = region->next) {
                    v1 = edgesThreeToFour[region->tet_vertex][0];
                    v2 = edgesThreeToFour[region->tet_vertex][1];
                    v3 = edgesThreeToFour[region->tet_vertex][2];

                    printf("    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                           " (Curves: [%d %d] [%d %d] [%d %d]) (Adj Curves: [%d %d] [%d %d] [%d %d]) (Dive: [%d %d] [%d %d] [%d %d])\n",
                           region->index, region->tet_index, region->tet_vertex,
                           region->adj_cusp_triangle[v1], region->adj_cusp_triangle[v2], region->adj_cusp_triangle[v3],
                           region->adj_cusp_regions[v1] == NULL ? -1 : region->adj_cusp_regions[v1]->index,
                           region->adj_cusp_regions[v2] == NULL ? -1 : region->adj_cusp_regions[v2]->index,
                           region->adj_cusp_regions[v3] == NULL ? -1 : region->adj_cusp_regions[v3]->index,
                           region->curve[v2][v1], region->curve[v3][v1],
                           region->curve[v1][v2], region->curve[v3][v2],
                           region->curve[v1][v3], region->curve[v2][v3],
                           region->num_adj_curves[v2][v1], region->num_adj_curves[v3][v1],
                           region->num_adj_curves[v1][v2], region->num_adj_curves[v3][v2],
                           region->num_adj_curves[v1][v3], region->num_adj_curves[v2][v3],
                           region->dive[v2][v1], region->dive[v3][v1],
                           region->dive[v1][v2], region->dive[v3][v2],
                           region->dive[v1][v3], region->dive[v2][v3]
                    );
                }
            }
        }

    } else if (strcmp(type, "homology") == 0) {
        printf("Homology info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            cusp = cusps[i];

            printf("Boundary %d\n", i);
            printf("Intersect Tet Index %d, Intersect Tet Vertex %d\n", cusp->intersect_tet_index, cusp->intersect_tet_vertex);
            printf("    Meridian\n");

            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[M][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][3]
                );
            }
            printf("    Longitude\n");
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[L][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][3]
                );
            }
        }

    } else if (strcmp(type, "edge_indices") == 0) {
        printf("Edge classes\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                v1 = edgesThreeToFour[tri->tet_vertex][0];
                v2 = edgesThreeToFour[tri->tet_vertex][1];
                v3 = edgesThreeToFour[tri->tet_vertex][2];

                printf("    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                       "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                       tri->tet_index, tri->tet_vertex,
                       v1, tri->vertices[v1].edge_class, tri->vertices[v1].edge_index,
                       v2, tri->vertices[v2].edge_class, tri->vertices[v2].edge_index,
                       v3, tri->vertices[v3].edge_class, tri->vertices[v3].edge_index
                );
            }
        }
    } else if (strcmp(type, "dual_curves") == 0) {
        printf("Oscillating curve paths\n");

        // which dual curve
        for (i = 0; i < curves->num_curves; i++) {
            j = 0;

            printf("Dual Curve %d\n", i);
            // which curve component
            for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
                printf("    Part %d: \n", j);

                for (path_node = path->path_begin.next;
                     path_node != &path->path_end;
                     path_node = path_node->next)
                    printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                           path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                           path_node->next_face, path_node->prev_face, path_node->inside_vertex
                    );
                j++;
            }
        }
    } else if (strcmp(type, "inside_edge") == 0) {
        printf("Inside edge info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("    (Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                       tri->tet_index,               // Tet Index
                       tri->tet_vertex,                // Tet Vertex
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][1]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][1]]
                );
            }
        }
    } else if (strcmp(type, "graph") == 0) {
        printf("Graph info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            cusp = cusps[i];

            printf("Boundary %d\n", i);
            g = cusp->dual_graph;
            for (j = 0; j < g->num_vertices; j++) {
                if (cusp->dual_graph_regions[j] == NULL)
                    continue;

                printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ", j,
                       cusp->dual_graph_regions[j]->tet_index,
                       cusp->dual_graph_regions[j]->tet_vertex
                );
                for (edge_node = g->edge_list_begin[j].next;
                     edge_node != &g->edge_list_end[j];
                     edge_node = edge_node->next)
                    printf("%d ", edge_node->y);

                printf("\n");
            }
        }
    } else if (strcmp(type, "endpoints") == 0) {
        printf("EndPoint Info\n");

        // which curve
        for (i = 0; i < curves->num_curves; i++) {
            printf("Dual Curve %d\n", i);

            j = 0;
            // which component
            for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
                printf("    Part %d Cusp %d\n", j, path->endpoints[0].tri->tet->cusp[path->endpoints[0].tri->tet_vertex]->index);
                for (k = 0; k < 2; k++) {
                    if (k == 0)
                        printf("        Start: ");
                    else
                        printf("        End:   ");

                    x_vertex1 = (int) remaining_face[path->endpoints[k].tri->tet_vertex][path->endpoints[k].vertex];
                    x_vertex2 = (int) remaining_face[path->endpoints[k].vertex][path->endpoints[k].tri->tet_vertex];

                    printf("Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d) Adj Curves %d\n",
                           path->endpoints[k].region_index, path->endpoints[k].tri->tet_index,
                           path->endpoints[k].tri->tet_vertex, path->endpoints[k].face, path->endpoints[k].vertex,
                           path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_class,
                           path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_index,
                           path->endpoints[k].num_adj_curves);
                }

                j++;
            }
        }
    } else {
        printf("Unknown type: %s\n", type);
    }
    printf("-------------------------------\n");
}

// ------------------------------------

void do_manifold_train_lines(Triangulation *manifold, CuspStructure **cusps, EndMultiGraph *multi_graph) {
    int cusp_index;
    EdgeClass *edge;
    Boolean *edge_class_on_cusp, *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, Boolean);

    // pick edge classes for train lines
    for (edge = manifold->edge_list_begin.next; edge != &manifold->edge_list_end; edge = edge->next) {
        if (multi_graph->edge_classes[edge->index] || multi_graph->e0 == edge->index) {
            edge_classes[edge->index] = TRUE;
        } else {
            edge_classes[edge->index] = FALSE;
        }
    }

    find_edge_class_edges(manifold, cusps, edge_classes);

    for (cusp_index = 0; cusp_index < manifold->num_cusps; cusp_index++) {
        edge_class_on_cusp = edge_classes_on_cusp(cusps[cusp_index], edge_classes);

        find_primary_train_line(cusps[cusp_index], edge_class_on_cusp);
        update_adj_curve_on_cusp(cusps[cusp_index]);
    }

    if (debug) {
        printf("\n");
        printf("Manifold Train Lines\n");
        printf("\n");
        printf("-------------------------------\n");

        log_structs(manifold, cusps, NULL, "train_lines");
        log_structs(manifold, cusps, NULL, "cusp_regions");
        log_structs(manifold, cusps, NULL, "graph");
    }
    my_free(edge_classes);
}

/*
 * Find a bipartite matching from the graph g which has a vertex for
 * each target edge class, a vertex for each tetrahedron and an
 * edge (tet, edge_class) iff edge_class corresponds to the edge index
 * of an edge of tet.
 *
 * edge_classes: array of booleans which are true for target edge classes
 */

int *find_tet_index_for_edge_classes(Triangulation *manifold, const Boolean *edge_classes) {
    int i, j, num_edge_classes = manifold->num_tetrahedra;
    int edge_source = 2 * num_edge_classes, tet_sink = 2 * num_edge_classes + 1;
    int *edge_class_to_tet_index = NEW_ARRAY(num_edge_classes, int);
    int **residual_network;
    Graph *g = init_graph(2 * num_edge_classes + 2, TRUE);
    Tetrahedron *tet;

    for (i = 0; i < num_edge_classes; i++)
        edge_class_to_tet_index[i] = -1;

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        for (i = 0; i < 6; i++) {
            if (edge_classes[tet->edge_class[i]->index])
                insert_edge(g, tet->edge_class[i]->index, num_edge_classes + tet->index, g->directed);
        }
    }

    /*
     * Convert the graph to a maximum flow problem
     */
    for (i = 0; i < num_edge_classes; i++)
        insert_edge(g, edge_source, i, g->directed);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next)
        insert_edge(g, num_edge_classes + tet->index, tet_sink, g->directed);

    residual_network = ford_fulkerson(g, edge_source, tet_sink);

    for (i = 0; i < num_edge_classes; i++) {
        for (j = num_edge_classes; j < 2 * num_edge_classes; j++)
            if (residual_network[j][i] == 1)
                edge_class_to_tet_index[i] = j - num_edge_classes;

        if (edge_classes[i] && edge_class_to_tet_index[i] == -1) {
            uFatalError("find_tet_index_for_edge_classes", "symplectic_basis");
        }
    }

    for (i = 0; i < 2 * num_edge_classes + 2; i++)
        my_free(residual_network[i]);

    free_graph(g);
    my_free(residual_network);
    return edge_class_to_tet_index;
}

/*
 * Assign a cusp triangle, face and vertex to each PathEndPoint of the train 
 * line. This is done in a breadth first search fashion from the first cusp, 
 * adding cusps to the search queue after diving through them.
 */

void find_edge_class_edges(Triangulation *manifold, CuspStructure **cusps, Boolean *edge_classes) {
    int edge_class, cusp_index, other_cusp_index;
    int *edge_class_to_tet_index = find_tet_index_for_edge_classes(manifold, edge_classes);
    Boolean found_edge_class;
    Boolean *visited_cusps, **edge_class_on_cusp = NEW_ARRAY(manifold->num_cusps, Boolean *);
    Queue *queue = init_queue(manifold->num_cusps);
    CuspStructure *cusp;

    for (cusp_index = 0; cusp_index < manifold->num_cusps; cusp_index++) {
        edge_class_on_cusp[cusp_index] = edge_classes_on_cusp(cusps[cusp_index], edge_classes);
    }

    enqueue(queue, 0);

    while (!empty_queue(queue)) {
        cusp_index = dequeue(queue);
        cusp = cusps[cusp_index];

        found_edge_class = FALSE;
        for (edge_class = 0; edge_class < cusp->num_edge_classes; edge_class++) {
            if (edge_class_on_cusp[cusp_index][edge_class])
                found_edge_class = TRUE;
        }

        if (!found_edge_class)
            continue;

        // assign edges to edge classes
        find_edge_class_edges_on_cusp(cusps[cusp_index], edge_class_on_cusp[cusp_index], edge_class_to_tet_index);

        // update dive edges classes
        visited_cusps = update_edge_classes_on_cusp(cusps, edge_class_on_cusp, manifold->num_cusps,
                                                    cusp->num_edge_classes,cusp_index);

        for (other_cusp_index = 0; other_cusp_index < manifold->num_cusps; other_cusp_index++) {
            if (!visited_cusps[other_cusp_index])
                continue;

            enqueue(queue, other_cusp_index);
        }

        my_free(visited_cusps);
    }

    for (cusp_index = 0; cusp_index < manifold->num_cusps; cusp_index++) {
        my_free(edge_class_on_cusp[cusp_index]);
    }

    my_free(edge_class_on_cusp);
    my_free(edge_class_to_tet_index);
    free_queue(queue);
}

/*
 * Find a cusp triangle, face and vertex for each edge class which is true
 * in edge_classes, using edge_class_to_tet_index to pick the tet for each
 * edge_class.
 */

void find_edge_class_edges_on_cusp(CuspStructure *cusp, const Boolean *edge_classes, const int *edge_class_to_tet_index) {
    int edge_class;
    VertexIndex v1, v2;
    FaceIndex face;
    CuspTriangle *tri;
    CuspVertex *vertex1, *vertex2;
    Boolean found;

    for (edge_class = 0; edge_class < cusp->num_edge_classes; edge_class++) {
        if (!edge_classes[edge_class])
            continue;

        if (edge_class_to_tet_index[edge_class] == -1)
            uFatalError("find_edge_class_edges_on_cusp", "symplectic_basis");

        found = FALSE;

        // find a cusp edge incident to the edge class
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            if (found || edge_class_to_tet_index[edge_class] != tri->tet_index)
                continue;

            for (face = 0; face < 4; face++) {
                if (face == tri->tet_vertex || found)
                    continue;

                v1 = remaining_face[tri->tet_vertex][face];
                v2 = remaining_face[face][tri->tet_vertex];

                vertex1 = &tri->vertices[v1];
                vertex2 = &tri->vertices[v2];

                if (vertex1->edge_class == edge_class && vertex1->edge_index == 0) {
                    cusp->train_line_endpoint[0][edge_class].tri = tri;
                    cusp->train_line_endpoint[0][edge_class].face = face;
                    cusp->train_line_endpoint[0][edge_class].vertex = v1;
                    found = TRUE;
                } else if (vertex2->edge_class == edge_class && vertex2->edge_index == 0) {
                    cusp->train_line_endpoint[0][edge_class].tri = tri;
                    cusp->train_line_endpoint[0][edge_class].face = face;
                    cusp->train_line_endpoint[0][edge_class].vertex = v2;
                    found = TRUE;
                }
            }
        }

        if (!found)
            uFatalError("find_edge_class_edges_on_cusp", "symplectic_basis");
    }
}

/*
 * Each edge we choose to add to the list of edges in find_edge_class_edges_on_cusp
 * has a corresponding edge on another cusp, which represents diving through the
 * manifold along that edge. Find these corresponding edges and set them in the
 * edge_begin and edges_end arrays, so the final train lines are consistent.
 */

Boolean *update_edge_classes_on_cusp(CuspStructure **cusps, Boolean **edge_classes,
                                     int num_cusps, int num_edge_classes, int current_cusp_index) {
    int cusp_index, other_cusp_index, edge_class, edge_index;
    VertexIndex v1, v2, vertex;
    CuspVertex *vertex1, *vertex2;
    Boolean *visited_cusp = NEW_ARRAY(num_edge_classes, Boolean);
    PathEndPoint *endpoint;
    CuspTriangle *tri, *other_tri;

    for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
        visited_cusp[cusp_index] = FALSE;
    }

    for (edge_class = 0; edge_class < num_edge_classes; edge_class++) {
        if (!edge_classes[current_cusp_index][edge_class])
            continue;

        endpoint = &cusps[current_cusp_index]->train_line_endpoint[0][edge_class];
        other_cusp_index = endpoint->tri->tet->cusp[endpoint->vertex]->index;

        if (other_cusp_index == current_cusp_index) {
            edge_index = 1;
        } else {
            edge_index = 0;
        }

        v1 = remaining_face[endpoint->tri->tet_vertex][endpoint->face];
        v2 = remaining_face[endpoint->face][endpoint->tri->tet_vertex];

        vertex1 = &endpoint->tri->vertices[v1];
        vertex2 = &endpoint->tri->vertices[v2];

        if (vertex1->edge_class == edge_class && vertex1->edge_index == 0)
            vertex = v1;
        else if (vertex2->edge_class == edge_class && vertex2->edge_index == 0)
            vertex = v2;
        else
            continue;

        other_tri = NULL;

        for (tri = cusps[other_cusp_index]->cusp_triangle_begin.next;
             tri != &cusps[other_cusp_index]->cusp_triangle_end;
             tri = tri->next) {
            if (tri->tet_vertex != vertex || tri->tet_index != endpoint->tri->tet_index)
                continue;

            other_tri = tri;
        }

        if (other_tri == NULL)
            uFatalError("update_edge_classes_on_cusp", "symplectic_basis");

        edge_classes[current_cusp_index][edge_class] = FALSE;
        edge_classes[other_cusp_index][edge_index * num_edge_classes + edge_class] = FALSE;
        visited_cusp[other_cusp_index] = TRUE;

        cusps[other_cusp_index]->train_line_endpoint[edge_index][edge_class].tri = other_tri;
        cusps[other_cusp_index]->train_line_endpoint[edge_index][edge_class].vertex = endpoint->tri->tet_vertex;
        cusps[other_cusp_index]->train_line_endpoint[edge_index][edge_class].face = endpoint->face;
    }

    return visited_cusp;
}

/*
 * edge_classes is a collection 'C' of edge classes indicated by TRUE in the array.
 * edge_classes_on_cusp returns C intersect { edge classes on 'cusp'} (the edge classes
 * which have an end lyine at 'cusp'.
 */

Boolean *edge_classes_on_cusp(CuspStructure *cusp, const Boolean *edge_classes) {
    CuspTriangle *tri;
    VertexIndex v;
    Boolean *edge_class_on_cusp = NEW_ARRAY(2 * cusp->manifold->num_tetrahedra, Boolean);
    int edge_class, edge_index;

    for (int i = 0; i < 2 * cusp->manifold->num_tetrahedra; i++) {
        edge_class_on_cusp[i] = FALSE;
    }

    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        for (v = 0; v < 4; v++) {
            if (v == tri->tet_vertex)
                continue;

            edge_class = tri->vertices[v].edge_class;
            edge_index = tri->vertices[v].edge_index;
            edge_class_on_cusp[cusp->manifold->num_tetrahedra * edge_index + edge_class] = edge_classes[edge_class];
        }
    }

    return edge_class_on_cusp;
}

// ------------------------------------

/*
 * Train lines
 */

/*
 * Use the regions on either side of the target edges to find a curve
 * through a cusp which passes along each target edge.
 */

void find_primary_train_line(CuspStructure *cusp, Boolean *edge_classes) {
    int start_index, start_class, finish_index, finish_class;
    PathEndPoint *start, *finish;
    Triangulation *manifold = cusp->manifold;

    start = next_valid_endpoint_index(cusp, NULL);
    tri_endpoint_to_region_endpoint(cusp, start);
    start_index = start->tri->vertices[start->vertex].edge_index;
    start_class = start->tri->vertices[start->vertex].edge_class;
    edge_classes[start_index * manifold->num_tetrahedra + start_class] = FALSE;

    if (!array_contains_true(edge_classes, 2 * manifold->num_tetrahedra)) {
        return;
    }

    finish = next_valid_endpoint_index(cusp, start);
    tri_endpoint_to_region_endpoint(cusp, finish);
    finish_index = finish->tri->vertices[finish->vertex].edge_index;
    finish_class = finish->tri->vertices[finish->vertex].edge_class;
    edge_classes[finish_index * manifold->num_tetrahedra + finish_class] = FALSE;
    do_initial_train_line_segment_on_cusp(cusp, start, finish);

    while (array_contains_true(edge_classes, 2 * manifold->num_tetrahedra)) {
        start = finish;
        finish = next_valid_endpoint_index(cusp, start);
        tri_endpoint_to_region_endpoint(cusp, finish);
        finish_index = finish->tri->vertices[finish->vertex].edge_index;
        finish_class = finish->tri->vertices[finish->vertex].edge_class;
        edge_classes[finish_index * manifold->num_tetrahedra + finish_class] = FALSE;

        do_train_line_segment_on_cusp(cusp, start, finish);
    }

    my_free(edge_classes);
}


/*
 * Construct the first segment of a train line. Essentially the same process
 * as do_curve_component_to_new_edge_class but stores the result in the cusp train
 * line.
 */

void do_initial_train_line_segment_on_cusp(CuspStructure *cusp, PathEndPoint *start_endpoint,
                                           PathEndPoint *finish_endpoint) {
    EdgeNode node_begin, node_end;

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    path_finding_with_loops(cusp, start_endpoint, finish_endpoint, 0, 0, &node_begin, &node_end);

    if (finish_endpoint == NULL)
        uFatalError("do_initial_train_line_segment_on_cusp", "symplectic_basis");

    // split along curve
    graph_path_to_dual_curve(cusp, &node_begin, &node_end,
                             &cusp->train_line_path_begin, &cusp->train_line_path_end,
                             start_endpoint, finish_endpoint);
    split_cusp_regions_along_path(cusp, &cusp->train_line_path_begin,
                                  &cusp->train_line_path_end, start_endpoint, finish_endpoint);

    start_endpoint->node = cusp->train_line_path_begin.next;
    finish_endpoint->node = cusp->train_line_path_end.prev;

    free_edge_node(&node_begin, &node_end);
}

/*
 * Construct the next train line segment after the first. The start endpoint
 * is already set, so we set the region the start endpoint is in to visited
 * before starting the breadth first search and instead start from the region
 * adjacent across the face of the cusp triangle we dive along.
 */

void do_train_line_segment_on_cusp(CuspStructure *cusp, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    EdgeNode node_begin, node_end;
    PathNode *start_node;
    CuspRegion *region;
    int start_index;

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    // update start_endpoint region
    start_index = TRI_TO_INDEX(start_endpoint->tri->tet_index, start_endpoint->tri->tet_vertex);
    for (region = cusp->cusp_region_begin[start_index].next;
         region != &cusp->cusp_region_end[start_index];
         region = region->next) {
        if (region->tet_index != start_endpoint->tri->tet_index ||
            region->tet_vertex != start_endpoint->tri->tet_vertex)
            continue;

        if (!region->adj_cusp_triangle[start_endpoint->face] ||
            !region->dive[start_endpoint->face][start_endpoint->vertex])
            continue;

        if (start_endpoint->face == cusp->train_line_path_end.prev->prev_face
            && region->curve[start_endpoint->face][start_endpoint->vertex] != 1)
            continue;

        start_endpoint->region_index = region->index;
        start_endpoint->region = region;
    }

    if (start_endpoint->region == NULL)
        uFatalError("do_train_line_segment_on_cusp", "symplectic_basis");

    /*
     * We require curves run between distinct sides of each cusp triangle
     * it enters. Hence, we need to remove the edge of the dual graph
     * corresponding to the last curve segment we drew. This edge will be
     * added back when the dual graph is reconstructed.
     */

    extended_train_line_path(cusp, start_endpoint, finish_endpoint, &node_begin, &node_end);

    if (finish_endpoint == NULL)
        uFatalError("do_initial_train_line_segment_on_cusp", "symplectic_basis");

    // split along curve
    start_node = cusp->train_line_path_end.prev;
    graph_path_to_path_node(cusp, &node_begin, &node_end,
                            &cusp->train_line_path_begin, &cusp->train_line_path_end, 
                            start_endpoint, finish_endpoint);
    split_cusp_regions_along_train_line_segment(cusp, start_node, &cusp->train_line_path_end,
                                                start_endpoint, finish_endpoint);

    finish_endpoint->node = cusp->train_line_path_end.prev;

    free_edge_node(&node_begin, &node_end);
}

PathEndPoint *next_valid_endpoint_index(CuspStructure *cusp, PathEndPoint *current_endpoint) {
    int start_index, start_class, edge_class;

    if (current_endpoint == NULL) {
        start_index = 0;
        start_class = -1;
    } else {
        start_index = current_endpoint->tri->vertices[current_endpoint->vertex].edge_index;
        start_class = current_endpoint->tri->vertices[current_endpoint->vertex].edge_class;
    }

    if (start_index == 0) {
        for (edge_class = start_class + 1; edge_class < cusp->num_edge_classes; edge_class++) {
            if (cusp->train_line_endpoint[START][edge_class].tri == NULL)
                continue;

            return &cusp->train_line_endpoint[START][edge_class];
        }

        start_class = -1;
    }

    for (edge_class = start_class + 1; edge_class < cusp->num_edge_classes; edge_class++) {
        if (cusp->train_line_endpoint[FINISH][edge_class].tri == NULL)
            continue;

        return &cusp->train_line_endpoint[FINISH][edge_class];
    }

    return NULL;
}

/*
 * Find a path from start_endpoint to finish_endpoint, which
 * goes around a cycle so the center is on the same side as the face
 * the finish endpoint dives along.
 */

void path_finding_with_loops(CuspStructure *cusp, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint,
                             int loop_edge_class, int loop_edge_index, EdgeNode *node_begin, EdgeNode *node_end) {
    int *parent;
    Boolean *discovered, *processed;

    construct_cusp_region_dual_graph(cusp);
    processed = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    init_search(cusp->dual_graph, processed, discovered, parent);
    bfs(cusp->dual_graph, start_endpoint->region_index, processed, discovered, parent);
    find_path(start_endpoint->region_index, finish_endpoint->region_index, parent, node_begin, node_end);

    my_free(processed);
    my_free(discovered);
    my_free(parent);
}

/*
 * Find a path for the train line from start to finish endpoint
 * and store in doubly linked list node_begin -> node_end.
 *
 * Uses the cycle ensured by path finding with loops to find
 * a path if the finish endpoint is not in the subgraph.
 */

void extended_train_line_path(CuspStructure *cusp, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint,
                              EdgeNode *node_begin, EdgeNode *node_end) {
    int cycle_start, cycle_end, start, finish, visited;
    Boolean cycle;
    Boolean *discovered, *processed;
    int *parent;

    construct_cusp_region_dual_graph(cusp);
    processed = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    if (start_endpoint->face == cusp->train_line_path_end.prev->prev_face) {
        // curve dives along the face it passes through

        start = start_endpoint->region_index;
        visited = start_endpoint->region->adj_cusp_regions[start_endpoint->face]->index;
    } else {
        // curve dives through the vertex opposite the face it passes through or
        // curve travells around the vertex it dives through

        start = start_endpoint->region->adj_cusp_regions[start_endpoint->face]->index;
        visited = start_endpoint->region_index;
    }

    finish = finish_endpoint->region_index;
    init_search(cusp->dual_graph, processed, discovered, parent);
    delete_edge(cusp->dual_graph, visited, start, cusp->dual_graph->directed);
    bfs(cusp->dual_graph, start, processed, discovered, parent);

    if (parent[finish] == -1 && start != finish) {
        /*
         * The finish endpoint is not in the subgraph we created by removing the edge
         * (visited, start). Assume there exists a cycle in this subgraph, we use this to
         * 'turn the curve around' and use the edge (visited, start).
         */

        init_search(cusp->dual_graph, processed, discovered, parent);
        cycle = cycle_exists(cusp->dual_graph, start, processed, discovered, parent, &cycle_start, &cycle_end);

        if (cycle == FALSE)
            // nothing we can do, train line does not work
            uFatalError("do_train_line_segment_on_cusp", "symplectic_basis");

        // reset parent array
        init_search(cusp->dual_graph, processed, discovered, parent);
        bfs(cusp->dual_graph, start, processed, discovered, parent);

        cycle_path(cusp->dual_graph, node_begin, node_end, start, visited,
                   finish, cycle_start, cycle_end);
    } else {
        find_path(start, finish, parent, node_begin, node_end);
    }

    my_free(processed);
    my_free(discovered);
    my_free(parent);
}

void cycle_path(Graph *g, EdgeNode *node_begin, EdgeNode *node_end, int start, int prev, int finish,
                int cycle_start, int cycle_end) {
    EdgeNode *node, *temp_node, temp_begin, temp_end;
    Boolean *discovered, *processed;
    int *parent;

    temp_begin.next = &temp_end;
    temp_begin.prev = NULL;
    temp_end.next = NULL;
    temp_end.prev = &temp_begin;

    processed = NEW_ARRAY(g->num_vertices, Boolean);
    discovered = NEW_ARRAY(g->num_vertices, Boolean);
    parent = NEW_ARRAY(g->num_vertices, int);

    // find a path from start -> cycle_end
    find_path(start, cycle_end, parent, node_begin, node_end);

    // duplicate the path start -> cycle_start, and reverse it
    find_path(start, cycle_start, parent, &temp_begin, &temp_end);
    for (node = temp_end.prev; node != &temp_begin; node = node->prev) {
        temp_node = NEW_STRUCT( EdgeNode );
        temp_node->y = node->y;
        INSERT_BEFORE(temp_node, node_end)
    }

    // find a path from visited -> target
    init_search(g, processed, discovered, parent);
    bfs(g, prev, processed, discovered, parent);
    find_path(prev, finish, parent, node_end->prev, node_end);

    free_edge_node(&temp_begin, &temp_end);
    my_free(processed);
    my_free(discovered);
    my_free(parent);
}

/*
 * Find a valid region for a path endpoint
 */

void tri_endpoint_to_region_endpoint(CuspStructure *cusp, PathEndPoint *endpoint) {
    CuspRegion *region;
    int index;

    if (endpoint == NULL || endpoint->tri == NULL)
        uFatalError("tri_endpoint_to_region_endpoint", "symplectic_basis");

    index = TRI_TO_INDEX(endpoint->tri->tet_index, endpoint->tri->tet_vertex);
    for (region = cusp->cusp_region_begin[index].next; region != &cusp->cusp_region_end[index]; region = region->next) {
        if (region->tet_index != endpoint->tri->tet_index || region->tet_vertex != endpoint->tri->tet_vertex)
            continue;

        if (!region->adj_cusp_triangle[endpoint->face] || !region->dive[endpoint->face][endpoint->vertex])
            continue;

        endpoint->region = region;
        endpoint->region_index = region->index;
    }

    if (endpoint->region == NULL)
        uFatalError("tri_endpoint_to_region_endpoint", "symplectic_basis");
}

Boolean array_contains_true(const Boolean *array, int len) {
    Boolean found = FALSE;

    for (int i = 0; i < len; i++) {
        if (array[i])
            found = TRUE;
    }

    return found;
}

/*
 * Convert the path found by BFS for a train line segment which is stored as
 * EdgeNode's in the node_begin -> node_end doubly linked list, to a path of
 * PathNode's int the path_begin -> path_end doubly linked list.
 */

void graph_path_to_path_node(CuspStructure *cusp, EdgeNode *node_begin, EdgeNode *node_end, PathNode *path_begin,
                             PathNode *path_end, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    FaceIndex face;
    EdgeNode *edge_node, *node;
    PathNode *path_node;
    CuspRegion *region;
    VertexIndex v1, v2;

    if (node_begin->next == node_end) {
        // path len 0
        return;
    } else if (node_begin->next->next == node_end) {
        // path len 1
        region = cusp->dual_graph_regions[node_begin->next->y];
        path_end->prev->next_face = start_endpoint->face;

        path_node = NEW_STRUCT( PathNode );
        INSERT_BEFORE(path_node, path_end)
        path_node->next_face = finish_endpoint->face;
        path_node->prev_face = EVALUATE(start_endpoint->tri->tet->gluing[start_endpoint->face], start_endpoint->face);
        path_node->cusp_region_index = node_begin->next->y;
        path_node->tri = region->tri;

        for (face = 0; face < 4; face++)
            if (region->tet_vertex != face &&
                path_node->next_face != face &&
                path_node->prev_face != face)
                break;

        path_node->inside_vertex = face;
        return;
    }

    // Set Header node
    path_end->prev->next_face = -1;

    // Add in a node for the start pos when the start endpoint is not in the same cusp tri as the first node.
    region = cusp->dual_graph_regions[node_begin->next->y];

    if (region->tet_index != start_endpoint->tri->tet_index || region->tet_vertex != start_endpoint->tri->tet_vertex) {
        node = NEW_STRUCT( EdgeNode );
        INSERT_AFTER(node, node_begin)
        node->y = start_endpoint->region_index;
    }

    for (face = 0; face < 4; face++) {
        if (!start_endpoint->region->adj_cusp_triangle[face])
            continue;

        if (start_endpoint->region->adj_cusp_regions[face]->index != cusp->dual_graph_regions[node_begin->next->next->y]->index)
            continue;

        path_end->prev->next_face = face;
    }

    if (path_end->prev->next_face == -1)
        uFatalError("graph_path_to_path_node", "symplectic_basis");

    v1 = remaining_face[start_endpoint->region->tet_vertex][path_end->prev->prev_face];
    v2 = remaining_face[path_end->prev->prev_face][start_endpoint->region->tet_vertex];

    if (path_end->prev->next_face == v1)
        path_end->prev->inside_vertex = v2;
    else
        path_end->prev->inside_vertex = v1;

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next)
        interior_edge_node_to_path_node(cusp->dual_graph_regions[edge_node->y], path_end, edge_node);

    // Set Tail node
    endpoint_edge_node_to_path_node(cusp->dual_graph_regions[edge_node->y], path_end, edge_node,
                                    finish_endpoint, FINISH);
}

/*
 * Split the cusp regions along the path path_begin -> path_end.
 * Handles the first node differently to split_cusp_regions_along_path,
 * due to the linking with the previous train line segment.
 */

void split_cusp_regions_along_train_line_segment(CuspStructure *cusp, PathNode *path_begin, PathNode *path_end, 
                                                 PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, split_type, region_index;
    PathNode *node;
    CuspRegion *p_region;
    Graph *g = cusp->dual_graph;

    if (path_begin->tri->tet_index == start_endpoint->tri->tet_index 
        && path_begin->tri->tet_vertex == start_endpoint->tri->tet_vertex) {
        node = path_begin;
    } else if (path_begin->next->tri->tet_index == start_endpoint->tri->tet_index 
               && path_begin->next->tri->tet_vertex == start_endpoint->tri->tet_vertex) {
        node = path_begin->prev;
    } else {
        uFatalError("split_cusp_regions_along_train_line_segment", "symplectic_basis");
        return;
    }

    if (node->next == path_end) {
        // empty path
        return ;
    }

    if (start_endpoint->face == node->prev_face) {
        // curve dives along the face it passes through
        split_type = 0;
    } else if (start_endpoint->vertex == node->prev_face) {
        // curve dives through the vertex opposite the face it passes through
        split_type = 1;
    } else {
        // curve travells around the vertex it dives through
        split_type = 2;
    }

    p_region = cusp->dual_graph_regions[start_endpoint->region_index];
    region_index = TRI_TO_INDEX(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_train_line_endpoints(&cusp->cusp_region_begin[region_index],
                                              &cusp->cusp_region_end[region_index],
                                              p_region, node, start_endpoint, START);
    split_cusp_region_train_line_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                          node, start_endpoint, index, split_type);
    index++;

    // interior edges
    for (node = node->next; node->next != path_end; node = node->next) {
        p_region = cusp->dual_graph_regions[node->cusp_region_index];
        region_index = TRI_TO_INDEX(p_region->tet_index, p_region->tet_vertex);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin[region_index],
                                           &cusp->cusp_region_end[region_index],
                                           p_region, node);
        split_cusp_region_path_interior(&cusp->cusp_region_end[region_index], p_region, node, index);
        index++;
    }

    // update last region
    p_region = cusp->dual_graph_regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   p_region, finish_endpoint, node, FINISH);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                    node, finish_endpoint, index, FINISH);
    index++;

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

/*
 * Updates the cusp regions for the cusp triangle where the train line segments
 * join.
 */

void update_cusp_triangle_train_line_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
                                               CuspRegion *region, PathNode *node, PathEndPoint *path_endpoint, int pos) {
    VertexIndex vertex1, vertex2;
    CuspRegion *current_region;

    vertex1 = remaining_face[region->tet_vertex][node->next_face];
    vertex2 = remaining_face[node->next_face][region->tet_vertex];

    for (current_region = cusp_region_start->next;
         current_region != cusp_region_end;
         current_region = current_region->next) {

        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (!current_region->adj_cusp_triangle[node->next_face])
            continue;

        // Curve goes around the vertex or passes through the face that carries it
        if (current_region->curve[node->next_face][vertex1] > region->curve[node->next_face][vertex1]) {
            current_region->curve[node->next_face][vertex1]++;
            current_region->dive[node->next_face][vertex1] = FALSE;

        } else if (current_region->curve[node->next_face][vertex1] < region->curve[node->next_face][vertex1]) {
            current_region->curve[node->next_face][vertex2]++;
            current_region->dive[node->next_face][vertex2] = FALSE;
        }
    }
}

/*
 * Split the cusp region where the train line segments join.
 */

void split_cusp_region_train_line_endpoint(CuspRegion *region_end, CuspRegion *region, PathNode *node,
                                           PathEndPoint *path_endpoint, int index, int split_type) {
    VertexIndex vertex1, vertex2, other_vertex;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    copy_region(region, new_region);
    new_region->index = index;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */

    if (split_type == 0) {
        if (node->next_face == path_endpoint->vertex) {
            // curve dives through the face opposite the next face
            other_vertex = (VertexIndex) (path_endpoint->face == vertex1 ? vertex2 : vertex1);

            new_region->curve[node->next_face][path_endpoint->face]++;
            new_region->dive[path_endpoint->face][other_vertex]   = region->dive[path_endpoint->face][other_vertex];
            new_region->dive[path_endpoint->vertex][other_vertex] = region->dive[path_endpoint->vertex][other_vertex];
            new_region->adj_cusp_triangle[other_vertex]           = FALSE;

            region->curve[node->next_face][other_vertex]++;
            region->dive[path_endpoint->face][other_vertex]       = FALSE;
            region->dive[path_endpoint->vertex][other_vertex]     = FALSE;
            region->adj_cusp_triangle[path_endpoint->face]        = FALSE;
        } else {
            new_region->curve[node->next_face][path_endpoint->face]++;
            new_region->dive[vertex1][path_endpoint->vertex]        = region->dive[vertex1][path_endpoint->vertex];
            new_region->dive[vertex2][path_endpoint->vertex]        = region->dive[vertex2][path_endpoint->vertex];
            new_region->adj_cusp_triangle[path_endpoint->face]      = FALSE;
            new_region->adj_cusp_triangle[path_endpoint->vertex]    = FALSE;

            region->curve[node->next_face][path_endpoint->vertex]++;
            region->dive[vertex1][path_endpoint->vertex]            = FALSE;
            region->dive[vertex2][path_endpoint->vertex]            = FALSE;
        }
    } else if (split_type == 1 || split_type == 2) {
        other_vertex = (VertexIndex) (path_endpoint->face == vertex1 ? vertex2 : vertex1);
        new_region->curve[path_endpoint->face][other_vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex]
                = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex]    = FALSE;
        new_region->adj_cusp_triangle[other_vertex]             = FALSE;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face][path_endpoint->vertex] = FALSE;
    } else
        uFatalError("split_cusp_region_train_line_endpoint", "symplectic_basis");

    INSERT_BEFORE(new_region, region_end)
}

// ------------------------------------

/*
 * Find oscillating curves. Each curve is made up of an even number of
 * components, with each component contained in a cusp, and connecting
 * two cusp vertices. Each oscillating curve is associated to an edge
 * of the triangulation, the rest of the edges come from the end multi
 * graph.
 *
 * The result is stored in tet->extra[edge_class].curve[f][v] array
 * on each tetrahedron.
 */

void do_oscillating_curves(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    CuspEndPoint cusp_path_begin, cusp_path_end, *temp_cusp;
    int i;

    cusp_path_begin.next = &cusp_path_end;
    cusp_path_begin.prev = NULL;
    cusp_path_end.next   = NULL;
    cusp_path_end.prev   = &cusp_path_begin;

    for (i = 0; i < curves->num_curves; i++) {
        find_multi_graph_path(cusps[0]->manifold, multi_graph,
                              &cusp_path_begin, &cusp_path_end, curves->edge_class[i]);
        do_one_oscillating_curve(cusps, curves, multi_graph, &cusp_path_begin, &cusp_path_end,
                                 curves->edge_class[i], i);

        while (cusp_path_begin.next != &cusp_path_end) {
            temp_cusp = cusp_path_begin.next;
            REMOVE_NODE(temp_cusp)
            my_free(temp_cusp);
        }

        if (debug) {
            printf("\n");
            printf("Oscillating Curve %d\n", i);
            printf("\n");
            printf("-------------------------------\n");

            log_structs(cusps[0]->manifold, cusps, curves, "dual_curves");
            log_structs(cusps[0]->manifold, cusps, curves, "endpoints");
            log_structs(cusps[0]->manifold, cusps, curves, "cusp_regions");
            log_structs(cusps[0]->manifold, cusps, curves, "graph");
        }
    }
}

/*
 * Construct a curve dual to the edge class 'edge_class'. The first and last 
 * components connect to edge_class which is not in the end multi graph so 
 * we need to find a new curve. Any intermediate components, if they exist, will
 * make use of the train lines, as they consist of curves between edge classes 
 * in the end multi graph and thus is a segment of the train line.
 */

void do_one_oscillating_curve(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph,
                              CuspEndPoint *cusp_path_begin, CuspEndPoint *cusp_path_end,
                              int edge_class, int curve_index) {
    int orientation = START;
    CuspEndPoint *endpoint = cusp_path_begin->next;
    CurveComponent *path,
                   *curve_begin = &curves->curve_begin[curve_index],
                   *curve_end = &curves->curve_end[curve_index];

    curve_begin->edge_class[FINISH] = edge_class;
    curve_end->edge_class[START]    = edge_class;

    if (cusps[endpoint->cusp_index]->train_line_endpoint[0][edge_class].tri == NULL) {
        path = setup_first_curve_component(cusps[endpoint->cusp_index], multi_graph, endpoint,
                                           curve_begin, curve_end);
        do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    } else {
        path = setup_train_line_component(cusps[endpoint->cusp_index], multi_graph, curve_begin, curve_end,
                                          endpoint, orientation);
        do_curve_component_on_train_line(cusps[path->cusp_index], path);
    }
    update_path_holonomy(path, edge_class);

    // interior curve components, coming from train lines
    for (endpoint = endpoint->next; endpoint->next != cusp_path_end; endpoint = endpoint->next) {
        orientation = (orientation == START ? FINISH : START);

        path = setup_train_line_component(cusps[endpoint->cusp_index], multi_graph, curve_begin, curve_end,
                                          endpoint, orientation);
        do_curve_component_on_train_line(cusps[path->cusp_index], path);
        update_path_holonomy(path, edge_class);
    }
    orientation = (orientation == START ? FINISH : START);

    if (cusps[endpoint->cusp_index]->train_line_endpoint[0][edge_class].tri == NULL) {
        path = setup_last_curve_component(cusps[endpoint->cusp_index], multi_graph, endpoint,
                                          curve_begin, curve_end);
        do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    } else {
        path = setup_train_line_component(cusps[endpoint->cusp_index], multi_graph, curve_begin, curve_end,
                                          endpoint, orientation);
        do_curve_component_on_train_line(cusps[path->cusp_index], path);
    }
    update_path_holonomy(path, edge_class);

    update_adj_curve_along_path(cusps, curves, curve_index,
                                (Boolean) (cusp_path_begin->next->next->next != cusp_path_end));
}

/*
 * Initialise the curve component for a path which lies on a train line.
 * Set the edge classes and copy path endpoints from the train line
 * endpoints.
 */

CurveComponent *setup_train_line_component(CuspStructure *cusp, EndMultiGraph *multi_graph, 
                                CurveComponent *curve_begin, CurveComponent *curve_end, 
                                CuspEndPoint *endpoint, int orientation) {
    CurveComponent *path; 

    path = init_curve_component(endpoint->edge_class[orientation],
                                endpoint->edge_class[orientation == START ? FINISH : START],
                                endpoint->cusp_index);

    INSERT_BEFORE(path, curve_end)

    if (cusp->train_line_endpoint[FINISH][path->edge_class[START]].tri != NULL && orientation == START) {
        COPY_PATH_ENDPOINT(&path->endpoints[START], &cusp->train_line_endpoint[FINISH][path->edge_class[START]])
    } else {
        COPY_PATH_ENDPOINT(&path->endpoints[START], &cusp->train_line_endpoint[START][path->edge_class[START]])
    }

    if (cusp->train_line_endpoint[FINISH][path->edge_class[FINISH]].tri != NULL && orientation == START) {
        COPY_PATH_ENDPOINT(&path->endpoints[FINISH], &cusp->train_line_endpoint[FINISH][path->edge_class[FINISH]])
    } else {
        COPY_PATH_ENDPOINT(&path->endpoints[FINISH], &cusp->train_line_endpoint[START][path->edge_class[FINISH]])
    }

    return path;
}

/*
 * Find a curve along the train line and copy it to 'curve'
 */

void do_curve_component_on_train_line(CuspStructure *cusp, CurveComponent *curve) {
    int orientation = 0;
    PathNode *node, *new_node, *start_node, *finish_node;
    FaceIndex temp_face;

    start_node = curve->endpoints[START].node;
    finish_node = curve->endpoints[FINISH].node;

    for (node = cusp->train_line_path_begin.next; node != &cusp->train_line_path_end; node = node->next) {
        if (node == start_node) {
            orientation = 1;
            break;
        } else if (node == finish_node) {
            orientation = -1;
            break;
        }
    }

    if (orientation == 1) {
        // copy the train line into the curve
        for (node = start_node; node != finish_node->next; node = node->next) {
            new_node = NEW_STRUCT( PathNode );
            COPY_PATH_NODE(new_node, node)
            INSERT_BEFORE(new_node, &curve->path_end)
        }
    } else if (orientation == -1) {
        // copy the train line into the curve
        for (node = start_node; node != finish_node->prev; node = node->prev) {
            new_node = NEW_STRUCT( PathNode );
            COPY_PATH_NODE(new_node, node)
            INSERT_BEFORE(new_node, &curve->path_end)

            // reverse direction of faces
            temp_face = new_node->next_face;
            new_node->next_face = new_node->prev_face;
            new_node->prev_face = temp_face;
        }
    } else
        uFatalError("do_curve_component_on_train_line", "symplectic_basis");

    // correct endpoint inside vertices
    curve->path_begin.next->prev_face = curve->endpoints[START].face;
    curve->path_end.prev->next_face = curve->endpoints[FINISH].face;

    if (curve->path_begin.next->prev_face == curve->path_begin.next->next_face) {
        curve->path_begin.next->inside_vertex = -1;
    }

    if (curve->path_end.prev->prev_face == curve->path_end.prev->next_face) {
        curve->path_end.prev->inside_vertex = -1;
    }
}

/*
 * Initalise the first curve component of an oscillating curve.
 * Set edge classes and find path endpoints.
 */

CurveComponent *setup_first_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                            CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START],
                                endpoint->edge_class[FINISH],
                                endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end)

    construct_cusp_region_dual_graph(cusp);
    find_single_endpoint(cusp, &path->endpoints[START],
                         path->edge_class[START], START);
    find_train_line_endpoint(cusp, &path->endpoints[FINISH], path->edge_class[FINISH],
                             START, multi_graph->e0, (Boolean) (endpoint->next->next != NULL));
    return path;
}

/*
 * Initalise the last curve component of an oscillating curve.
 * Set edge classes and find path endpoints.
 */

CurveComponent *setup_last_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                           CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START],
                                endpoint->edge_class[FINISH],
                                endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end)

    construct_cusp_region_dual_graph(cusp);
    find_single_matching_endpoint(cusp,
                                  &curves_begin->next->endpoints[START],
                                  &path->endpoints[START],
                                  path->edge_class[START], FINISH);

    find_single_matching_endpoint(cusp,
                                  &path->prev->endpoints[FINISH],
                                  &path->endpoints[FINISH],
                                  path->edge_class[FINISH], FINISH);

    return path;
}

/*
 * Construct an oscillating curve component, which is either the
 * first or last component of an oscillating curve.
 */

void do_curve_component_to_new_edge_class(CuspStructure *cusp, CurveComponent *curve) {
    int *parent;
    Boolean *processed, *discovered;
    EdgeNode node_begin, node_end;

    processed   = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered  = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent      = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    // Find curve using bfs
    init_search(cusp->dual_graph, processed, discovered, parent);
    bfs(cusp->dual_graph, curve->endpoints[START].region_index, processed, discovered, parent);

    find_path(curve->endpoints[START].region_index, curve->endpoints[FINISH].region_index,
              parent, &node_begin, &node_end);
    graph_path_to_dual_curve(cusp, &node_begin, &node_end,
                             &curve->path_begin, &curve->path_end,
                             &curve->endpoints[START], &curve->endpoints[FINISH]);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the curve
    split_cusp_regions_along_path(cusp, &curve->path_begin, &curve->path_end,
                                  &curve->endpoints[START], &curve->endpoints[FINISH]);

    free_edge_node(&node_begin, &node_end);
}


/*
 * Find a cusp region which can dive along a face into a vertex of
 * the cusp triangle which corresponds to 'edge_class' and 'edge_index',
 * and store the result in path_endpoint.
 */

void find_single_endpoint(CuspStructure *cusp, PathEndPoint *path_endpoint, int edge_class, int edge_index) {
    int i;
    VertexIndex vertex;
    FaceIndex face1, face2, face;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < cusp->dual_graph->num_vertices; i++) {
        if (cusp->dual_graph_regions[i] == NULL) {
            continue;
        }

        region = cusp->dual_graph_regions[i];
        // which vertex to dive through
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == region->tet_vertex)
                continue;

            if (region->tri->vertices[vertex].edge_class != edge_class)
                continue;

            if (region->tri->vertices[vertex].edge_index != edge_index)
                continue;

            face1 = remaining_face[region->tet_vertex][vertex];
            face2 = remaining_face[vertex][region->tet_vertex];

            if (region->dive[face1][vertex])
                face = face1;
            else if (region->dive[face2][vertex])
                face = face2;
            else
                continue;

            path_endpoint->region           = region;
            path_endpoint->tri              = region->tri;
            path_endpoint->vertex           = vertex;
            path_endpoint->face             = face;
            path_endpoint->region_index     = i;
            path_endpoint->num_adj_curves   = region->num_adj_curves[path_endpoint->face][path_endpoint->vertex];

            return ;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_endpoints", "symplectic_basis");
}

/*
 * Find a cusp region which can dive into a vertex of the cusp triangle
 * corresponding 'edge_class' and 'edge_index', while matching path_endpoint1.
 *
 * See 'region_index', 'region_vertex', 'region_dive', 'region_curve' for the
 * conditions for a matching endpoint.
 */

void find_single_matching_endpoint(CuspStructure *cusp, PathEndPoint *path_endpoint1, PathEndPoint *path_endpoint2,
                                   int edge_class, int edge_index) {
    int i;
    Boolean region_index, region_vertex, region_dive, region_curve;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < cusp->dual_graph->num_vertices; i++) {
        if (cusp->dual_graph_regions[i] == NULL)
            continue;

        region = cusp->dual_graph_regions[i];

        // are we in the matching endpoint
        region_index    = (Boolean) (region->tet_index != path_endpoint1->tri->tet_index);
        region_vertex   = (Boolean) (region->tet_vertex != path_endpoint1->vertex);
        region_dive     = (Boolean) !region->dive[path_endpoint1->face][path_endpoint1->tri->tet_vertex];
        region_curve    = (Boolean) (region->num_adj_curves[path_endpoint1->face][path_endpoint1->tri->tet_vertex]
                != path_endpoint1->num_adj_curves);

        if (region_index || region_vertex || region_dive || region_curve)
            continue;

        path_endpoint2->region          = region;
        path_endpoint2->tri             = region->tri;
        path_endpoint2->vertex          = path_endpoint1->tri->tet_vertex;
        path_endpoint2->face            = path_endpoint1->face;
        path_endpoint2->region_index    = i;
        path_endpoint2->num_adj_curves  = region->num_adj_curves[path_endpoint2->face][path_endpoint2->vertex];

        return ;
    }

    // didn't find valid path endpoints
    uFatalError("find_single_matching_endpoints", "symplectic_basis");
}

/*
 * find a path endpoint which matches the train line endpoint found during
 * do_manifold_train_lines().
 */

void find_train_line_endpoint(CuspStructure *cusp, PathEndPoint *endpoint, int edge_class, int edge_index,
                                       int e0, Boolean is_train_line) {
    int i;
    Boolean region_index, region_vertex, region_dive, region_curve;
    CuspRegion *region;
    PathEndPoint *train_line_endpoint = &cusp->train_line_endpoint[edge_index][edge_class];

    // which cusp region
    for (i = 0; i < cusp->dual_graph->num_vertices; i++) {
        if (cusp->dual_graph_regions[i] == NULL)
            continue;

        region = cusp->dual_graph_regions[i];
        region_index    = (Boolean) (region->tet_index != train_line_endpoint->tri->tet_index);
        region_vertex   = (Boolean) (region->tet_vertex != train_line_endpoint->tri->tet_vertex);
        region_dive     = (Boolean) !region->dive[train_line_endpoint->face][train_line_endpoint->vertex];

        if (is_train_line) {
            region_curve = (Boolean) (region->num_adj_curves[train_line_endpoint->face][train_line_endpoint->vertex] !=
                    train_line_endpoint->num_adj_curves);
        } else {
            region_curve = (Boolean) (region->num_adj_curves[train_line_endpoint->face][train_line_endpoint->vertex] != 0);
        }

        if (region_index || region_vertex || region_dive || region_curve)
            continue;

        endpoint->region          = region;
        endpoint->tri             = region->tri;
        endpoint->vertex          = train_line_endpoint->vertex;
        endpoint->face            = train_line_endpoint->face;
        endpoint->region_index    = region->index;
        endpoint->num_adj_curves  = region->num_adj_curves[endpoint->face][endpoint->vertex];

        return ;
    }

    uFatalError("find_train_line_endpoint", "symplectic_basis");
}

/*
 * After finding a path, each node contains the index of the region it lies in. 
 * Update path info calculates the face the path crosses to get to the next node 
 * and the vertex it cuts off to simplify combinatorial holonomy calculation.
 */

void graph_path_to_dual_curve(CuspStructure *cusp, EdgeNode *node_begin, EdgeNode *node_end, PathNode *path_begin,
                              PathNode *path_end, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    FaceIndex face;
    EdgeNode *edge_node;
    PathNode *path_node;
    CuspRegion *region;

    // path len 0
    if (node_begin->next == node_end)
        return;

    edge_node = node_begin->next;
    // path len 1
    if (edge_node->next == node_end) {
        for (face = 0; face < 4; face++)
            if (cusp->dual_graph_regions[edge_node->y]->tet_vertex != face &&
                start_endpoint->vertex != face &&
                finish_endpoint->vertex != face)
                break;

        region = cusp->dual_graph_regions[edge_node->y];

        path_node = NEW_STRUCT( PathNode );
        INSERT_BEFORE(path_node, path_end)
        path_node->next_face = finish_endpoint->face;
        path_node->prev_face = start_endpoint->face;
        path_node->cusp_region_index = edge_node->y;
        path_node->tri = region->tri;
        path_node->inside_vertex = face;
        return;
    }

    // Set Header node
    endpoint_edge_node_to_path_node(cusp->dual_graph_regions[edge_node->y], path_end, edge_node,
                                    start_endpoint, START);

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next)
        interior_edge_node_to_path_node(cusp->dual_graph_regions[edge_node->y], path_end, edge_node);

    // Set Tail node
    endpoint_edge_node_to_path_node(cusp->dual_graph_regions[edge_node->y], path_end, edge_node,
                                    finish_endpoint, FINISH);
}

void endpoint_edge_node_to_path_node(CuspRegion *region, PathNode *path_end, EdgeNode *edge_node,
                                     PathEndPoint *path_endpoint, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        path_node->next_face = -1;
        for (face = 0; face < 4; face++) {
            if (face == region->tet_vertex || !region->adj_cusp_triangle[face] || path_node->next_face != -1)
                continue;

        if (region->adj_cusp_regions[face]->index == edge_node->next->y)
            path_node->next_face = face;
        }

        // next node isn't in an adjacent region
        if (path_node->next_face == -1)
            uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

        path_node->prev_face = path_endpoint->face;

        if (path_node->next_face == path_endpoint->vertex) {
            if (path_endpoint->face == vertex1)
                path_node->inside_vertex = vertex2;
            else
                path_node->inside_vertex = vertex1;
        } else if (path_node->next_face == path_endpoint->face) {
            path_node->inside_vertex = -1;
        } else {
            path_node->inside_vertex = path_endpoint->vertex;
        }
    } else {
        path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face],
                                        path_end->prev->next_face);
        path_node->next_face = path_endpoint->face;

        if (path_node->prev_face == path_endpoint->vertex) {
            if (path_endpoint->face == vertex1)
                path_node->inside_vertex = vertex2;
            else
                path_node->inside_vertex = vertex1;
        } else if (path_node->prev_face == path_endpoint->face) {
            path_node->inside_vertex = -1;
        } else {
            path_node->inside_vertex = path_endpoint->vertex;
        }
    }

    INSERT_BEFORE(path_node, path_end)
}

/*
 * node lies in 'region', find the vertex which the subpath 
 * node->prev->y --> node->y --> node->next->y cuts off of the cusp triangle 
 * >tri.
 */

void interior_edge_node_to_path_node(CuspRegion *region, PathNode *path_end, EdgeNode *edge_node) {
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face],
                                    path_end->prev->next_face);

    vertex1 = remaining_face[path_node->tri->tet_vertex][path_node->prev_face];
    vertex2 = remaining_face[path_node->prev_face][path_node->tri->tet_vertex];

    if (region->adj_cusp_triangle[vertex1] && region->adj_cusp_regions[vertex1]->index == edge_node->next->y) {
        path_node->next_face = vertex1;
        path_node->inside_vertex = vertex2;
    } else if (region->adj_cusp_triangle[vertex2] && region->adj_cusp_regions[vertex2]->index == edge_node->next->y) {
        path_node->next_face = vertex2;
        path_node->inside_vertex = vertex1;
    } else
        uFatalError("interior_edge_node_to_path_node", "symplectic_basis");

    INSERT_BEFORE(path_node, path_end)
}

/*
 * The oscillating curve splits the region it passes through into two regions. 
 * Split each region in two and update attributes
 */

void split_cusp_regions_along_path(CuspStructure *cusp, PathNode *path_begin, PathNode *path_end,
                                   PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, region_index;
    PathNode *node;
    CuspRegion *region;

    // empty path
    if (path_begin->next == path_end)
        return ;

    // path of len 1
    if (path_begin->next->next == path_end) {
        split_path_len_one(cusp, path_begin->next, start_endpoint, finish_endpoint);
        return;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and region to the left of the curve.
     */
    node = path_begin->next;
    region = cusp->dual_graph_regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   region, start_endpoint, node, START);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], region,
                                    node, start_endpoint, index, START);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        region = cusp->dual_graph_regions[node->cusp_region_index];
        region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin[region_index],
                                           &cusp->cusp_region_end[region_index], region, node);
        split_cusp_region_path_interior(&cusp->cusp_region_end[region_index], region, node, index);
        index++;
    }

    // update last region
    region = cusp->dual_graph_regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   region, finish_endpoint, node, FINISH);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], region,
                                    node, finish_endpoint, index, FINISH);
    index++;

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

void split_path_len_one(CuspStructure *cusp, PathNode *node, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, region_index;
    FaceIndex face;
    CuspRegion *new_region, *old_region, *region;

    new_region = NEW_STRUCT(CuspRegion);
    old_region = cusp->dual_graph_regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(old_region->tet_index, old_region->tet_vertex);
    INSERT_BEFORE(new_region, &cusp->cusp_region_end[region_index])
    copy_region(old_region, new_region);

    face = node->inside_vertex;

    new_region->index = index;
    new_region->adj_cusp_triangle[start_endpoint->vertex]                   = FALSE;
    new_region->adj_cusp_triangle[finish_endpoint->vertex]                  = FALSE;
    new_region->dive[face][start_endpoint->vertex]                          = TRUE;
    new_region->dive[face][finish_endpoint->vertex]                         = TRUE;
    new_region->dive[start_endpoint->vertex][finish_endpoint->vertex]       = (Boolean) (face != finish_endpoint->face);
    new_region->dive[finish_endpoint->vertex][start_endpoint->vertex]       = (Boolean) (face != start_endpoint->face);
    new_region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
    new_region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;

    old_region->adj_cusp_triangle[face]             = FALSE;
    old_region->dive[face][start_endpoint->vertex]  = (Boolean) (face == start_endpoint->face);
    old_region->dive[face][finish_endpoint->vertex] = (Boolean) (face == finish_endpoint->face);
    old_region->temp_adj_curves[face][start_endpoint->vertex]++;
    old_region->temp_adj_curves[face][finish_endpoint->vertex]++;

    // update other cusp regions
    for (region = cusp->cusp_region_begin[region_index].next;
         region != &cusp->cusp_region_end[region_index];
         region = region->next) {

        if (new_region->tet_index != region->tet_index || new_region->tet_vertex != region->tet_vertex)
            continue;

        if (region == new_region || region == old_region)
            continue;

        if (region->adj_cusp_triangle[start_endpoint->vertex] || region->adj_cusp_triangle[finish_endpoint->vertex]) {
            region->temp_adj_curves[face][finish_endpoint->vertex]++;
            region->temp_adj_curves[face][start_endpoint->vertex]++;

        } else {
            region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
            region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;
        }
    }

    update_adj_region_data(cusp);
    cusp->num_cusp_regions++;
}

/*
 * Set the new and old region data. Draw a picture to see how the attributes 
 * change in each case
 */

void split_cusp_region_path_interior(CuspRegion *region_end, CuspRegion *region, PathNode *node, int index) {
    int v1, v2;
    CuspRegion *new_region = NEW_STRUCT( CuspRegion );

    v1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    v2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    /*
     * new_region becomes the cusp region closest to the inside vertex and
     * region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    // Update new region
    new_region->curve[v1][v2]++;
    new_region->curve[v2][v1]++;
    new_region->dive[v1][node->inside_vertex]           = region->dive[v1][node->inside_vertex];
    new_region->dive[v2][node->inside_vertex]           = region->dive[v2][node->inside_vertex];
    new_region->adj_cusp_triangle[node->inside_vertex]  = FALSE;

    // Update region
    region->curve[v1][node->inside_vertex]++;
    region->curve[v2][node->inside_vertex]++;
    region->dive[v1][node->inside_vertex]           = FALSE;
    region->dive[v2][node->inside_vertex]           = FALSE;

    INSERT_BEFORE(new_region, region_end)
}

void split_cusp_region_path_endpoint(CuspRegion *region_end, CuspRegion *region, PathNode *path_node, 
                                     PathEndPoint *path_endpoint, int index, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;
    path_endpoint->region = NULL;

    if (pos == START) {
        face = path_node->next_face;
    } else {
        face = path_node->prev_face;
    }

    if (face == path_endpoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        new_region->curve[path_endpoint->vertex][vertex2]++;
        new_region->temp_adj_curves[vertex1][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]      = (Boolean) (path_endpoint->face == vertex1);
        new_region->dive[vertex2][path_endpoint->vertex]      = region->dive[vertex2][path_endpoint->vertex];
        new_region->dive[vertex2][vertex1]                    = region->dive[vertex2][vertex1];
        new_region->dive[path_endpoint->vertex][vertex1]      = region->dive[path_endpoint->vertex][vertex1];
        new_region->adj_cusp_triangle[vertex1]                = FALSE;

        region->curve[path_endpoint->vertex][vertex1]++;
        region->temp_adj_curves[vertex2][path_endpoint->vertex]++;
        region->dive[vertex2][path_endpoint->vertex]         = (Boolean) (path_endpoint->face == vertex2);
        region->dive[vertex2][vertex1]                       = FALSE;
        region->dive[path_endpoint->vertex][vertex1]         = FALSE;
        region->adj_cusp_triangle[vertex2]                   = FALSE;
    } else if (face == path_endpoint->face) {
        // curve passes through the face that carries it
        new_region->curve[path_endpoint->face][path_endpoint->face == vertex1 ? vertex2 : vertex1]++;
        new_region->temp_adj_curves[face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex]
                        = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex]                                 = FALSE;
        new_region->adj_cusp_triangle[path_endpoint->face == vertex1 ? vertex2 : vertex1]    = FALSE;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
    } else {
        // Curve goes around the vertex
        new_region->curve[face][path_endpoint->face]++;
        new_region->temp_adj_curves[path_endpoint->face][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]              = region->dive[vertex1][path_endpoint->vertex];
        new_region->dive[vertex2][path_endpoint->vertex]              = region->dive[vertex2][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->face]            = FALSE;
        new_region->adj_cusp_triangle[path_endpoint->vertex]          = FALSE;

        region->curve[face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex] = FALSE;
    }

    INSERT_BEFORE(new_region, region_end)
}

/*
 * After splitting each region the path travels through, the attributes for 
 * other regions in the same cusp triangle is now out of date. Update cusp 
 * triangles for nodes in the interior of the path.
 */

void update_cusp_triangle_path_interior(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
                                        CuspRegion *region, PathNode *node) {
    int face1, face2;
    CuspRegion *current_region;

    face1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    face2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    for (current_region = cusp_region_start->next;
         current_region != cusp_region_end;
         current_region = current_region->next) {

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (current_region->curve[face1][node->inside_vertex] > region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][node->inside_vertex]++;
        }
        else if (current_region->curve[face1][node->inside_vertex] < region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][face2]++;
        }

        if (current_region->curve[face2][node->inside_vertex] > region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][node->inside_vertex]++;
        }
        else if (current_region->curve[face2][node->inside_vertex] < region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][face1]++;
        }
    }
}

/*
 * After splitting each curveRegion the path travels through, the attributes 
 * for other regions in the same cusp triangle is now out of date. Update cusp 
 * triangles for nodes at the end of the path.
 */

void update_cusp_triangle_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end, CuspRegion *region,
                                    PathEndPoint *path_endpoint, PathNode *node, int pos) {
    FaceIndex face, face1, face2;
    CuspRegion *current_region;

    face1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    face2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        face = node->next_face;
    } else {
        face = node->prev_face;
    }

    for (current_region = cusp_region_start->next;
         current_region != cusp_region_end;
         current_region = current_region->next) {
        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (face == path_endpoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (!current_region->adj_cusp_triangle[face]) {
                if (!current_region->adj_cusp_triangle[face1]) {
                    current_region->temp_adj_curves[face1][path_endpoint->vertex]++;
                } else if (!current_region->adj_cusp_triangle[face2]) {
                    current_region->temp_adj_curves[face2][path_endpoint->vertex]++;
                } else {
                    uFatalError("update_cusp_triangle_endpoints", "symplectic_basis");
                }
            } else if (current_region->curve[path_endpoint->vertex][face1] > region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face1]++;
                current_region->temp_adj_curves[face2][path_endpoint->vertex]++;
            } else if (current_region->curve[path_endpoint->vertex][face1] < region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face2]++;
                current_region->temp_adj_curves[face1][path_endpoint->vertex]++;
            }

            continue;
        }

        if (!current_region->adj_cusp_triangle[face]) {
            current_region->temp_adj_curves[face][path_endpoint->vertex]++;
            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (current_region->curve[face][path_endpoint->vertex] > region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][path_endpoint->vertex]++;
            current_region->temp_adj_curves[face][path_endpoint->vertex]++;

        } else if (current_region->curve[face][path_endpoint->vertex] < region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][face == face1 ? face2 : face1]++;
            current_region->temp_adj_curves[face == face1 ? face2 : face1][path_endpoint->vertex]++;
        }
    }
}

void update_adj_curve_along_path(CuspStructure **cusps, OscillatingCurves *curves, int curve_index, Boolean train_line) {
    int cusp_index, edge_class, edge_index;
    CurveComponent *curve,
               *dual_curve_begin = &curves->curve_begin[curve_index],
               *dual_curve_end   = &curves->curve_end[curve_index];
    CuspStructure *cusp;
    Triangulation *manifold = cusps[0]->manifold;

    // Update regions curve data
    for (curve = dual_curve_begin->next; curve != dual_curve_end; curve = curve->next)
        update_adj_curve_on_cusp(cusps[curve->cusp_index]);

    // update train line endpoints
    for (cusp_index = 0; cusp_index < manifold->num_cusps; cusp_index++) {
        cusp = cusps[cusp_index];

        for (edge_class = 0; edge_class < manifold->num_tetrahedra; edge_class++) {
            for (edge_index = 0; edge_index < 2; edge_index++) {
                if (cusp->train_line_endpoint[edge_index][edge_class].tri == NULL)
                    continue;

                update_adj_curve_at_endpoint(&cusp->train_line_endpoint[edge_index][edge_class],
                                             dual_curve_begin->next, START);
                if (train_line == FALSE)
                    update_adj_curve_at_endpoint(&cusp->train_line_endpoint[edge_index][edge_class],
                                                 dual_curve_end->prev, FINISH);
            }
        }
    }
}

/*
 * curve_begin and curve_end are header and tailer nodes of a doubly linked list of path
 * components for a new path. Update the path_endpoint->num_adj_curves attribute to account for this
 * new curve.
 */

void update_adj_curve_at_endpoint(PathEndPoint *path_endpoint, CurveComponent *path, int pos) {
    PathEndPoint *curve_end_point;

    curve_end_point = &path->endpoints[pos];

    if (curve_end_point->tri->tet_index != path_endpoint->tri->tet_index ||
        curve_end_point->tri->tet_vertex != path_endpoint->tri->tet_vertex ||
        curve_end_point->face != path_endpoint->face ||
        curve_end_point->vertex != path_endpoint->vertex)
        return;

    path_endpoint->num_adj_curves++;
}

/*
 * Move the temp adj curves into the current num of adj curves.
 */

void update_adj_curve_on_cusp(CuspStructure *cusp) {
    int i, j, k;
    CuspRegion *region;

    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            // which cusp region
            for (j = 0; j < 4; j++) {
                for (k = 0; k < 4; k++) {
                    region->num_adj_curves[j][k] += region->temp_adj_curves[j][k];
                    region->temp_adj_curves[j][k] = 0;
                }
            }
        }
    }
}

void update_path_holonomy(CurveComponent *path, int edge_class) {
    PathNode *path_node;

    for (path_node = path->path_begin.next; path_node != &path->path_end; path_node = path_node->next) {
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->next_face]++;
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->prev_face]--;
    }
}

