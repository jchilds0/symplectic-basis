/*
 * graph.c
 *
 * Basic graph implementation using adjacency list,
 * and breadth first search for oscillating_curves.c and end_multi_graph.c
 */


#include "symplectic_kernel.h"

int         augment_path(Graph *, int **, Boolean *, int, int, int);
Boolean     contains(int *, int, int);

/*
 * Initialise the arrays of the graph 'g' to their default values
 */

Graph *init_graph(int max_vertices, Boolean directed) {
    int i;
    Graph *g = NEW_STRUCT(Graph);

    g->num_vertices = max_vertices;
    g->directed = directed;

    g->edge_list_begin      = NEW_ARRAY(max_vertices, EdgeNode);
    g->edge_list_end        = NEW_ARRAY(max_vertices, EdgeNode);
    g->degree               = NEW_ARRAY(max_vertices, int);
    g->color                = NEW_ARRAY(max_vertices, int);

    for (i = 0; i < max_vertices; i++) {
        g->degree[i] = 0;
        g->color[i] = -1;

        g->edge_list_begin[i].next     = &g->edge_list_end[i];
        g->edge_list_begin[i].prev     = NULL;
        g->edge_list_end[i].next       = NULL;
        g->edge_list_end[i].prev       = &g->edge_list_begin[i];
    }

    return g;
}

void free_graph(Graph *g) {
    if (g == NULL)
        return;

    for (int i = 0; i < g->num_vertices; i++) {
        free_edge_node(&g->edge_list_begin[i], &g->edge_list_end[i]);
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->color);
    my_free(g);
}

/*
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(Graph *g, int x, int y, Boolean directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    EdgeNode *p = NEW_STRUCT( EdgeNode);
    INSERT_AFTER(p, &g->edge_list_begin[x]);
    p->y = y;
    g->degree[x]++;

    if (!directed) {
        insert_edge(g, y, x, TRUE);
    }

    return x;
}

/*
 * Remove the edge from vertex x to vertex y
 */

void delete_edge(Graph *g, int vertex_x, int vertex_y, Boolean directed) {
    EdgeNode *node;

    for (node = g->edge_list_begin[vertex_x].next;
         node != &g->edge_list_end[vertex_x] && node->y != vertex_y;
         node = node->next);

    if (node == &g->edge_list_end[vertex_x])
        return;

    REMOVE_NODE(node)
    my_free(node);

    if (!directed) {
        delete_edge(g, vertex_y, vertex_x, TRUE);
    }
}

/*
 * Check if an edge already exists in the graph
 */

Boolean edge_exists(Graph *g, int v1, int v2) {
    EdgeNode *node = &g->edge_list_begin[v1];

    while ((node = node->next)->next != NULL) {
        if (node->y == v2) {
            return TRUE;
        }
    }

    return FALSE;
}

/*
 * Initialise default values for bfs arrays
 */

void init_search(Graph *g, Boolean *processed, Boolean *discovered, int *parent) {
    int i;

    for (i = 0; i < g->num_vertices; i ++) {
        processed[i] = FALSE;
        discovered[i] = FALSE;
        parent[i] = -1;
    }
}

/*
 * Graph search algorithm starting at vertex 'start'.
 */

void bfs(Graph *g, int start, Boolean *processed, Boolean *discovered, int *parent) {
    Queue *q = init_queue(10);
    int v, y;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = TRUE;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = TRUE;
        p = &g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
}

/*
 * Recover the path through the graph from the parents array and store
 * in the doubly linked list node_begin -> ... -> node_end.
 */

void find_path(int start, int end, int *parents, EdgeNode *node_begin, EdgeNode *node_end) {
    int u;

    if (start != end && parents[end] == -1) {
        uFatalError("find_path", "symplectic_basis");
    }

    u = end;
    while (u != start) {
        EdgeNode *new_node = NEW_STRUCT(EdgeNode);
        new_node->y = u;
        INSERT_AFTER(new_node, node_begin);

        u = parents[u];
    };

    EdgeNode *new_node = NEW_STRUCT(EdgeNode);
    new_node->y = start;
    INSERT_AFTER(new_node, node_begin);
}

/*
 * Find a cycle
 */

Boolean cycle_exists(Graph *g, int start, Boolean *processed, Boolean *discovered,
                     int *parent, int *cycle_start, int *cycle_end) {
    Queue *q = init_queue(10);
    int v, y;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = TRUE;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = TRUE;
        p = &g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (processed[y] && y != parent[v]) {
                free_queue(q);
                *cycle_start = y;
                *cycle_end = v;
                return TRUE;
            }

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
    *cycle_start = 0;
    *cycle_end = 0;
    return FALSE;
}

int **ford_fulkerson(Graph *g, int source, int sink) {
    int i, j, augment;
    int **residual_network = NEW_ARRAY(g->num_vertices, int *);
    EdgeNode *node;
    Boolean *visited = NEW_ARRAY(g->num_vertices, Boolean);

    for (i = 0; i < g->num_vertices; i++) {
        residual_network[i] = NEW_ARRAY(g->num_vertices, int);

        for (j = 0; j < g->num_vertices; j++)
            residual_network[i][j] = -1;
    }


    for (i = 0; i < g->num_vertices; i++) {
        for (node = g->edge_list_begin[i].next; node != &g->edge_list_end[i]; node = node->next) {
            residual_network[i][node->y] = 1;
            residual_network[node->y][i] = 0;
        }
    }

    augment = 1;
    while (augment > 0) {
        for (i = 0; i < g->num_vertices; i++)
            visited[i] = FALSE;

        augment = augment_path(g, residual_network, visited, source, sink, INT_MAX);
    }

    my_free(visited);
    return residual_network;
}

int augment_path(Graph *g, int **residual_network, Boolean *visited, int u, int t, int bottleneck) {
    int residual, augment, v;
    EdgeNode *node;

    if (u == t)
        return bottleneck;

    visited[u] = TRUE;

    for (node = g->edge_list_begin[u].next; node != &g->edge_list_end[u]; node = node->next) {
        v = node->y;
        residual = residual_network[u][v];

        if (residual > 0 && !visited[v]) {
            augment = augment_path(g, residual_network, visited, v, t, MIN(bottleneck, residual));

            if (augment > 0) {
                residual_network[u][v] -= augment;
                residual_network[v][u] += augment;
                return augment;
            }
        }
    }

    return 0;
}

/*
 * Modified breadth first search which starts at vertex 'start' and
 * searches until we find an element of 'targets'.
 */

int bfs_target_list(Graph *g, int start, int *targets, int num_targets, Boolean *processed,
                    Boolean *discovered, int *parent) {
    Queue *q = init_queue(10);
    Boolean found = FALSE;
    int v, y, retval = -1;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = TRUE;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = TRUE;
        p = &g->edge_list_begin[v];

        if (!found && contains(targets, num_targets, v)) {
            found = TRUE;
            retval = v;
        }

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
    return retval;
}

Boolean contains(int *array, int len, int target) {
    Boolean found = FALSE;

    for (int i = 0; i < len; i++)
        if (array[i] == target)
            found = TRUE;

    return found;
}

void free_edge_node(EdgeNode *node_begin, EdgeNode *node_end) {
    EdgeNode *node;

    while (node_begin->next != node_end) {
        node = node_begin->next;
        REMOVE_NODE(node);
        my_free(node);
    }
}

/*
 * Find a spanning tree of graph1
 */

Graph *spanning_tree(Graph *graph1, int start, int *parent) {
    int i;

    Boolean *processed = NEW_ARRAY(graph1->num_vertices, Boolean);
    Boolean *discovered = NEW_ARRAY(graph1->num_vertices, Boolean);

    Graph *graph2 = init_graph(graph1->num_vertices, graph1->directed);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->num_vertices; i++) {
        if (parent[i] == -1)
            continue;

        insert_edge(graph2, i, parent[i], graph2->directed);
    }

    my_free(processed);
    my_free(discovered);

    return graph2;
}
