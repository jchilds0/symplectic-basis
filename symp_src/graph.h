//
// Created by joshu on 5/09/2023.
//

#ifndef GRAPH_H
#define GRAPH_H

#include "kernel.h"

typedef struct EdgeNode {
    int                         y;
    struct EdgeNode             *next;
    struct EdgeNode             *prev;
} EdgeNode;

typedef struct Graph {
    EdgeNode                    *edge_list_begin;        /** header node of doubly linked list */
    EdgeNode                    *edge_list_end;          /** tailer node ... */
    int                         *degree;                 /** degree of each vertex */
    int                         *color;                  /** color a tree bipartite */
    int                         num_vertices;            /** number of vertices in the graph */
    Boolean                     directed;                /** is the graph directed */
} Graph;

Graph                   *init_graph(int, Boolean);
void                    free_graph(Graph *);
int                     insert_edge(Graph *, int, int, Boolean);
void                    delete_edge(Graph *, int, int, Boolean);
Boolean                 edge_exists(Graph *, int, int);

void                    init_search(Graph *, Boolean *, Boolean *, int *);
void                    bfs(Graph *, int, Boolean *, Boolean *, int *);
void                    find_path(int, int, int *, EdgeNode *, EdgeNode *);
Boolean                 cycle_exists(Graph *, int, Boolean *, Boolean *, int *, int *, int *);
int                     **ford_fulkerson(Graph *, int, int);
int                     augment_path(Graph *, int **, Boolean *, int, int, int);
int                     bfs_target_list(Graph *, int, int *, int, Boolean *, Boolean *, int *);
Boolean                 contains(int *, int, int);
void                    free_edge_node(EdgeNode *, EdgeNode *);

#endif /* GRAPH_H */
