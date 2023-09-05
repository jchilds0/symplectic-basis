/*
 * cusp_train_lines.c
 *
 * Provides the functions
 *
 *      void                    init_train_line(CuspStructure *);
 *
 *      void                    free_train_line(CuspStructure *);
 *
 *      void                    do_manifold_train_lines(Triangulation *, CuspStructure **, EndMultiGraph *);
 *
 * which are used by cusp_structure.c when initialising the cusp structure.
 * It also provides
 *
 *      CurveComponent          *setup_train_line_component(CuspStructure *, EndMultiGraph *, CurveComponent *, CurveComponent *, CuspEndPoint *, int);
 *
 *      void                    do_curve_component_on_train_line(CuspStructure *, CurveComponent *);
 *
 * which are used by oscillating_curves.c to construct oscillating
 * curves on the train lines.
 */

#include "symplectic_kernel.h"

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

void free_train_line(CuspStructure *cusp) {
    PathNode *path_node;

    while (cusp->train_line_path_begin.next != &cusp->train_line_path_end) {
        path_node = cusp->train_line_path_begin.next;
        REMOVE_NODE(path_node)
        my_free(path_node);
    }

    my_free(cusp->train_line_endpoint[0]);
    my_free(cusp->train_line_endpoint[1]);

}

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

//  -----------------

// Oscillating Curves

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
