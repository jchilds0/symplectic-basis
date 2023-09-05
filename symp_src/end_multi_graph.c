/*
 * end_multi_graph.c
 *
 * Provides the functions,
 *
 *      EndMultiGraph           *init_end_multi_graph(Triangulation *);
 *
 *      void                    free_end_multi_graph(EndMultiGraph *);
 *
 * which are used by oscillating_curve.c to initialise the end multi graph.
 * It also provides the function
 *
 *      void find_multi_graph_path(Triangulation *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int);
 *
 * which is used by oscillating_curves.c to find a path through the end
 * multi graph.
 */


#include "symplectic_kernel.h"

int                     **find_end_multi_graph_edge_classes(EndMultiGraph *, Triangulation *);
int                     find_edge_class(Triangulation *, int, int);
void                    cusp_graph(Triangulation *, Graph *);
void                    color_graph(Graph *);
int                     find_same_color_edge(Triangulation *, EndMultiGraph *, Graph *);
int                     find_path_len(int, int, int *, int);
void                    graph_path_to_cusp_path(EndMultiGraph *, EdgeNode *, EdgeNode *, CuspEndPoint *, CuspEndPoint *, int);
void                    find_edge_ends(Graph *, Triangulation *, int, int *, int *);

/*
 * End Multi Graph
 *
 * The end multi graph is a graph with vertices for the cusps of M and
 * edges for each edge of the triangulation. We also refer to the spanning
 * tree of this graph as the end multi graph. The end multi graph structure also
 * keeps track of a special E_0 edge, which is used to construct a path of even
 * length through the graph.
 */

EndMultiGraph *init_end_multi_graph(Triangulation *manifold) {
    int i, j;
    int *parent;
    EndMultiGraph *multi_graph = NEW_STRUCT( EndMultiGraph );

    multi_graph->num_cusps = manifold->num_cusps;
    multi_graph->num_edge_classes = manifold->num_tetrahedra;

    Graph *g = init_graph(multi_graph->num_cusps, FALSE);
    cusp_graph(manifold, g);

    parent = NEW_ARRAY(g->num_vertices, int);

    multi_graph->multi_graph = spanning_tree(g, 0, parent);
    color_graph(multi_graph->multi_graph);

    multi_graph->edges = find_end_multi_graph_edge_classes(multi_graph, manifold);
    multi_graph->e0 = find_same_color_edge(manifold, multi_graph, g);

    multi_graph->edge_classes = NEW_ARRAY(multi_graph->num_edge_classes, Boolean);
    for (i = 0; i < multi_graph->num_edge_classes; i++) {
        multi_graph->edge_classes[i] = FALSE;
    }

    for (i = 0; i < multi_graph->num_cusps; i++) {
        for (j = 0; j < multi_graph->num_cusps; j++) {
            if (multi_graph->edges[i][j] == -1)
                continue;

            multi_graph->edge_classes[multi_graph->edges[i][j]] = TRUE;
        }
    }

    free_graph(g);
    my_free(parent);
    return multi_graph;
}

void free_end_multi_graph(EndMultiGraph *multi_graph) {
    int i;

    free_graph(multi_graph->multi_graph);

    for (i = 0; i < multi_graph->num_cusps; i++)
        my_free(multi_graph->edges[i]);

    my_free(multi_graph->edge_classes);
    my_free(multi_graph->edges);
    my_free(multi_graph);
}

void cusp_graph(Triangulation *manifold, Graph *g) {
    int vertex1, vertex2;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (vertex1 = 0; vertex1 < 4; vertex1++) {
            // which vertex of the cusp triangle at vertex1
            for (vertex2 = 0; vertex2 < 4; vertex2++) {
                if (vertex1 == vertex2)
                    continue;

                insert_edge(g, tet->cusp[vertex1]->index, tet->cusp[vertex2]->index, g->directed);
            }
        }
    }
}

/*
 * Assign an edge class to each edge of the graph g and return an array of
 * Booleans indicating if an edge class is in the graph.
 */

int **find_end_multi_graph_edge_classes(EndMultiGraph *multi_graph, Triangulation *manifold) {
    int i, j, edge_class, **cusps;
    EdgeNode *edge_node;
    Graph *g = multi_graph->multi_graph;

    cusps = NEW_ARRAY(multi_graph->num_cusps, int *);

    for (i = 0; i < multi_graph->num_cusps; i++) {
        cusps[i] = NEW_ARRAY(multi_graph->num_cusps, int);

        for (j = 0; j < multi_graph->num_cusps; j++)
            cusps[i][j] = -1;
    }

    for (i = 0; i < g->num_vertices; i++) {
        for (edge_node = g->edge_list_begin[i].next; edge_node != &g->edge_list_end[i]; edge_node = edge_node->next) {
            edge_class = find_edge_class(manifold, i, edge_node->y);
            cusps[i][edge_node->y] = edge_class;
            cusps[edge_node->y][i] = edge_class;
        }
    }

    return cusps;
}

/*
 * Find an edge class whose edge connects cusp1 and cusp2
 */

int find_edge_class(Triangulation *manifold, int cusp1, int cusp2) {
    int v1, v2;
    EdgeClass *edge;
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        for (v1 = 0; v1 < 4; v1++) {
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                if (tet->cusp[v1]->index != cusp1 || tet->cusp[v2]->index != cusp2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                return edge->index;
            }
        }
    }

    uFatalError("find_edge_class", "symplectic_basis");
    return 0;
}

void color_graph(Graph *g) {
    int color = 0, v;
    Queue *q = init_queue(g->num_vertices);
    EdgeNode *node;

    g->color[0] = color;
    q = enqueue(q, 0);

    while (!empty_queue(q)) {
        v = dequeue(q);
        color = g->color[v];

        for (node = g->edge_list_begin[v].next; node != &g->edge_list_end[v]; node = node->next) {
            // graph is not bipartite
            if (g->color[node->y] == color)
                uFatalError("color_graph", "symplectic_basis");

            if (g->color[node->y] != -1)
                continue;

            g->color[node->y] = !color;
            q = enqueue(q, node->y);
        }
    }

    free_queue(q);
}

/*
 * g1 is the colored spanning tree of g2, return the edge class of the edge in
 * g2 which connects vertices in g1 of the same color
 */

int find_same_color_edge(Triangulation *manifold, EndMultiGraph *multi_graph, Graph *g2) {
    int cusp;
    EdgeNode *node;
    Graph *g1 = multi_graph->multi_graph;

    for (cusp = 0; cusp < g2->num_vertices; cusp++) {
        for (node = g2->edge_list_begin[cusp].next; node != &g2->edge_list_end[cusp]; node = node->next) {
            if (g1->color[cusp] == g1->color[node->y] && multi_graph->edges[cusp][node->y] == -1)
                // we found an edge
                return find_edge_class(manifold, cusp, node->y);
        }
    }

    // we didn't find an edge connecting vertices of the same color
    uFatalError("find_same_color_edge", "symplectic_basis");
    return -1;
}

/*
 * Find the length of a path between start and end
 */

int find_path_len(int start, int end, int *parents, int path_length) {
    if ((start == end) || (end == -1)) {
        return path_length;
    } else {
        return find_path_len(start, parents[end], parents, path_length + 1);
    }
}

/*
 * Find a path through the end multi graph, starting at 'edge_class' and which
 * an odd length, since this corresponds to an even number of oscillating curve
 * components. The path is stored in the doubly linked list cusp_path_begin ->
 * cusp_path_end.
 */

void find_multi_graph_path(Triangulation *manifold, EndMultiGraph *multi_graph, CuspEndPoint *cusp_path_begin,
                           CuspEndPoint *cusp_path_end, int edge_class) {
    Graph *g = multi_graph->multi_graph;
    Boolean *processed     = NEW_ARRAY(g->num_vertices, Boolean);
    Boolean *discovered    = NEW_ARRAY(g->num_vertices, Boolean);
    int *parent         = NEW_ARRAY(g->num_vertices, int);
    int start, end, startE0, endE0, path_len = 0;
    EdgeNode node_begin, node_end;

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    find_edge_ends(g, manifold, edge_class, &start, &end);
    find_edge_ends(g, manifold, multi_graph->e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    path_len = find_path_len(start, end, parent, path_len);

    if (path_len % 2 == 1) {
        find_path(start, end, parent, &node_begin, &node_end);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, &node_begin, &node_end);

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end.prev, &node_end);
    }

    graph_path_to_cusp_path(multi_graph, &node_begin, &node_end, cusp_path_begin, cusp_path_end, edge_class);

    free_edge_node(&node_begin, &node_end);
    my_free(parent);
    my_free(discovered);
    my_free(processed);
}

/*
 * Converts the EdgeNode path through the cusps in the end multigraph to
 * a CuspEndPoint path which contains the edge classes on each cusp.
 * A CuspEndPoint corresponding to one section of an oscillating curve, and
 * constructing such a section for all CuspEndPoints gives the whole curve.
 *
 * node_begin -> node_end is a doubly linked list through the end multi graph
 * and the result path is stored in the doubly linked list cusp_path_begin ->
 * cusp_path_end.
 */

void graph_path_to_cusp_path(EndMultiGraph *multi_graph, EdgeNode *node_begin, EdgeNode *node_end,
                             CuspEndPoint *cusp_path_begin, CuspEndPoint *cusp_path_end, int edge_class) {
    int cusp, prev_edge_class;
    EdgeNode *node;
    CuspEndPoint *endpoint;

    prev_edge_class = edge_class;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        endpoint = NEW_STRUCT( CuspEndPoint );
        INSERT_BEFORE(endpoint, cusp_path_end);

        endpoint->cusp_index = cusp;
        endpoint->edge_class[START] = prev_edge_class;
        endpoint->edge_class[FINISH] = multi_graph->edges[node->y][node->next->y];

        if (endpoint->edge_class[FINISH] == -1)
            endpoint->edge_class[FINISH] = multi_graph->e0;

        prev_edge_class = endpoint->edge_class[FINISH];
    }

    endpoint = NEW_STRUCT( CuspEndPoint );
    INSERT_BEFORE(endpoint, cusp_path_end);

    endpoint->cusp_index = node->y;
    endpoint->edge_class[START] = prev_edge_class;
    endpoint->edge_class[FINISH] = edge_class;
}

void find_edge_ends(Graph *g, Triangulation *manifold, int edge_class, int *start, int *end) {
    int v1, v2;
    Tetrahedron *tet;
    EdgeClass *edge;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (v1 = 0; v1 < 4; v1++) {
            // which vertex of the cusp triangle at v1
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                if (edge->index != edge_class)
                    continue;

                *start = tet->cusp[v1]->index;
                *end   = tet->cusp[v2]->index;
                return;
            }
        }
    }


    // didn't find the edge class in the graph
    uFatalError("find_edge_ends", "symplectic_basis");
}
