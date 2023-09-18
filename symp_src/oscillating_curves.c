/*
 * oscillating_curves.c
 *
 * Provides the function
 *
 *      void do_oscillating_curves(Triangulation *manifold, Boolean *edge_classes);
 *
 * which finds oscillating curves and stores the intersection nums. in the
 * extra attribute of manifold. It also exposes a number of other functions used to
 * split cusp regions and update cusp regions for use in cusp_train_lines.c
 */

#include "symplectic_kernel.h"

OscillatingCurves       *init_oscillating_curves(Triangulation *, const Boolean *);
void                    free_oscillating_curves(OscillatingCurves *);

void                    find_oscillating_curves(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph);
void                    do_one_oscillating_curve(CuspStructure **, OscillatingCurves *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int, int);
CurveComponent          *setup_first_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
CurveComponent          *setup_last_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
void                    do_curve_component_to_new_edge_class(CuspStructure *, CurveComponent *);

void                    find_single_endpoint(CuspStructure *, PathEndPoint *, int, int);
void                    find_single_matching_endpoint(CuspStructure *, PathEndPoint *, PathEndPoint *, int, int);
void                    find_train_line_endpoint(CuspStructure *, PathEndPoint *, int, int, int, Boolean);

void                    split_path_len_one(CuspStructure *, PathNode *, PathEndPoint *, PathEndPoint *);

void                    update_adj_curve_along_path(CuspStructure **, OscillatingCurves *, int, Boolean);
void                    update_adj_curve_at_endpoint(PathEndPoint *, CurveComponent *, int);
void                    update_path_holonomy(CurveComponent *, int);

CurveComponent *init_curve_component(int edge_class_start, int edge_class_finish, int cusp_index) {
    CurveComponent *path = NEW_STRUCT(CurveComponent );

    *path = (CurveComponent) {
        {edge_class_start, edge_class_finish},
        cusp_index,
        (PathNode) {
            -1, -1, -1, -1,
            NULL, &path->path_end, NULL
        },
        (PathNode) {
            -1, -1, -1, -1,
            NULL, NULL, &path->path_begin
        },
        {
            (PathEndPoint) {
                -1, -1, -1, 0,
                NULL, NULL, NULL
            },
            (PathEndPoint) {
                -1, -1, -1, 0,
                NULL, NULL, NULL
            }
        },
        NULL,
        NULL
    };

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

/*
 * Initialise cusp structure on each cusp, construct train lines, construct
 * oscillating curves and store the intersection numbers of each curve with the
 * cusp triangles it enters in tet->extra[edge_class]->curve, in the same fashion
 * as the peripheral curves.
 */

void do_oscillating_curves(Triangulation *manifold, Boolean *edge_classes) {
    int i;
    char buf[200];
    CuspStructure **cusps         = NEW_ARRAY(manifold->num_cusps, CuspStructure *);
    EndMultiGraph *multi_graph    = init_end_multi_graph(manifold);
    Cusp *cusp;

    for (i = 0; i < multi_graph->num_edge_classes; i++)
        edge_classes[i] = multi_graph->edge_classes[i] == TRUE ? FALSE : TRUE;

    edge_classes[multi_graph->e0] = FALSE;

    OscillatingCurves *curves   = init_oscillating_curves(manifold, edge_classes);

    log_structs(NULL, NULL, NULL, "Generating Cusp Structures");

    for (i = 0; i < manifold->num_cusps; i++) {
        for (cusp = manifold->cusp_list_begin.next; cusp != &manifold->cusp_list_end && cusp->index != i; cusp = cusp->next);

        if (cusp == &manifold->cusp_list_end)
            uFatalError("do_oscillating_curves", "symplectic_basis");

        cusps[i] = init_cusp_structure(manifold, cusp, multi_graph);
    }

    log_structs(manifold, cusps, NULL, "cusp structure");

    do_manifold_train_lines(manifold, cusps, multi_graph);
    find_oscillating_curves(cusps, curves, multi_graph);

    log_structs(NULL, NULL, NULL, "Oscillating Curves Complete");

    for (i = 0; i < manifold->num_cusps; i++) {
        memset(buf, '\0', sizeof buf);
        sprintf(buf, "Cusp %d: %d Cusp Regions", i, cusps[i]->num_cusp_regions);
        log_structs(NULL, NULL, NULL, buf);
    }

    if (debug) {
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("%d, ", cusps[i]->num_cusp_regions);
        }
        printf("\n");
    }

    for (i = 0; i < manifold->num_cusps; i++)
        free_cusp_structure(cusps[i]);

    free_end_multi_graph(multi_graph);
    free_oscillating_curves(curves);
    my_free(cusps);
}

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

void find_oscillating_curves(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    char buf[100];
    CuspEndPoint cusp_path_begin, cusp_path_end, *temp_cusp;
    int i;

    cusp_path_begin.next = &cusp_path_end;
    cusp_path_begin.prev = NULL;
    cusp_path_end.next   = NULL;
    cusp_path_end.prev   = &cusp_path_begin;

    for (i = 0; i < curves->num_curves; i++) {
        memset(buf, '\0', sizeof buf);
        sprintf(buf, "Constructing Oscillating Curve %d", i);
        log_structs(NULL, NULL, NULL, buf);

        find_multi_graph_path(cusps[0]->manifold, multi_graph,
                              &cusp_path_begin, &cusp_path_end, curves->edge_class[i]);
        do_one_oscillating_curve(cusps, curves, multi_graph, &cusp_path_begin, &cusp_path_end,
                                 curves->edge_class[i], i);

        while (cusp_path_begin.next != &cusp_path_end) {
            temp_cusp = cusp_path_begin.next;
            REMOVE_NODE(temp_cusp)
            my_free(temp_cusp);
        }

        log_structs(cusps[0]->manifold, cusps, curves, "oscillating curve");
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

            *path_endpoint = (PathEndPoint) {
                face,
                vertex,
                i,
                region->num_adj_curves[face][vertex],
                NULL,
                region,
                region->tri
            };

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

        *path_endpoint2 = (PathEndPoint) {
            path_endpoint1->face,
            path_endpoint1->tri->tet_vertex,
            i,
            region->num_adj_curves[path_endpoint2->face][path_endpoint2->vertex],
            NULL,
            region,
            region->tri
        };

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

        *endpoint = (PathEndPoint) {
            train_line_endpoint->face,
            train_line_endpoint->vertex,
            region->index,
            region->num_adj_curves[train_line_endpoint->face][train_line_endpoint->vertex],
            NULL,
            region,
            region->tri
        };

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

        *path_node = (PathNode) {
            edge_node->y,
            finish_endpoint->face,
            start_endpoint->face,
            face,
            region->tri,
            path_node->next,
            path_node->prev
        };

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

